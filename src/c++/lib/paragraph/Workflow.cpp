// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2017 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 * \brief Workflow implementation
 *
 * \file Workflow.cpp
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>

#include <boost/filesystem.hpp>

#include "common/Error.hh"
#include "common/Threads.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/Workflow.hh"

namespace paragraph
{

Workflow::Workflow(
    bool jointInputs, const InputPaths& inputPaths, const InputPaths& inputIndexPaths,
    const std::vector<std::string>& graph_spec_paths, const std::string& output_file_path,
    const std::string& output_folder_path, const Parameters& parameters, const std::string& reference_path,
    const std::string& target_regions)
    : graphSpecPaths_(graph_spec_paths)
    , outputFilePath_(output_file_path)
    , outputFolderPath_(output_folder_path)
    , parameters_(parameters)
    , referencePath_(reference_path)
    , targetRegions_(target_regions)
{
    if (jointInputs)
    {
        unprocessedInputs_.push_back(Input(inputPaths, inputIndexPaths, std::begin(graphSpecPaths_)));
    }
    else
    {
        for (std::size_t i = 0; inputPaths.size() > i; ++i)
        {
            const auto& inputPath = inputPaths[i];
            const auto& inputIndexPath = inputIndexPaths.at(i);
            unprocessedInputs_.push_back(
                Input(InputPaths(1, inputPath), InputPaths(1, inputIndexPath), std::begin(graphSpecPaths_)));
        }
    }
}

static void dumpOutput(const std::string& output, std::ostream& os, const std::string& file)
{
    os << output;
    if (!os)
    {
        error("ERROR: Failed to write output to '%s' error: '%s'", file.c_str(), std::strerror(errno));
    }
}

std::string Workflow::processGraph(
    const std::string& graphSpecPath, const Parameters& parameters, const InputPaths& inputPaths,
    std::vector<common::BamReader>& readers)
{
    common::ReadBuffer allReads;
    for (common::BamReader& reader : readers)
    {
        common::extractReads(reader, parameters.target_regions(), (int)(parameters.max_reads()), allReads);
    }
    Json::Value outputJson = alignAndDisambiguate(parameters, allReads);
    if (inputPaths.size() == 1)
    {
        outputJson["bam"] = inputPaths.front();
    }
    else
    {
        outputJson["bam"] = Json::arrayValue;
        for (const auto& inputPath : inputPaths)
        {
            outputJson["bam"].append(inputPath);
        }
    }
    return Json::FastWriter().write(outputJson);
}

void Workflow::makeOutputFile(const std::string& output, const std::string& graphSpecPath)
{
    boost::filesystem::path inputPath(graphSpecPath);
    boost::filesystem::path outputPath = boost::filesystem::path(outputFolderPath_) / inputPath.filename();
    std::ofstream ofs(outputPath.string());
    if (!ofs)
    {
        error("ERROR: Failed to open output file '%s'. Error: '%s'", outputPath.string().c_str(), std::strerror(errno));
    }
    dumpOutput(output, ofs, outputPath.string());
}

void Workflow::processGraphs(std::ofstream& outputFileStream)
{
    for (Input& input : unprocessedInputs_)
    {
        std::vector<common::BamReader> readers;
        //        for (const std::string& inputPath : input.inputPaths_)
        for (size_t i = 0; i != input.inputPaths_.size(); ++i)
        {
            //            LOG()->info("Opening {} with {}", inputPath, referencePath_);
            //            readers.push_back(common::BamReader(inputPath, referencePath_));
            const auto& bamPath = input.inputPaths_[i];
            const auto& bamIndexPath = input.inputIndexPaths_[i];
            LOG()->info("Opening {}/{} with {}", bamPath, bamIndexPath, referencePath_);
            readers.emplace_back(bamPath, bamIndexPath, referencePath_);
        }

        std::lock_guard<std::mutex> lock(mutex_);
        while (graphSpecPaths_.end() != input.unprocessedGraphs_)
        {
            const std::string& graphSpecPath = *(input.unprocessedGraphs_++);
            std::string output;
            ASYNC_BLOCK_WITH_CLEANUP([this](bool failure) { terminate_ |= failure; })
            {
                if (terminate_)
                {
                    LOG()->warn("terminating");
                    break;
                }
                common::unlock_guard<std::mutex> unlock(mutex_);
                Parameters parameters = parameters_;
                LOG()->info("Loading parameters {}", graphSpecPath);
                parameters.load(graphSpecPath, referencePath_, targetRegions_);
                LOG()->info("Done loading parameters");

                output = processGraph(graphSpecPath, parameters, input.inputPaths_, readers);

                if (!outputFolderPath_.empty())
                {
                    makeOutputFile(output, graphSpecPath);
                }
            }

            if ("-" == outputFilePath_)
            {
                dumpOutput(output, std::cout, "standard output");
            }
            else if (!outputFilePath_.empty())
            {
                dumpOutput(output, outputFileStream, outputFilePath_);
            }
        }
    }
}

void Workflow::run()
{
    std::ofstream outputFileStream;
    if (!outputFilePath_.empty())
    {
        if ("-" != outputFilePath_)
        {
            LOG()->info("Output file path: {}", outputFilePath_);
            outputFileStream.open(outputFilePath_);
            if (!outputFileStream)
            {
                error(
                    "ERROR: Failed to open output file '%s'. Error: '%s'", outputFilePath_.c_str(),
                    std::strerror(errno));
            }
        }
        else
        {
            LOG()->info("Output to stdout");
        }
    }

    common::CPU_THREADS(parameters_.threads()).execute([this, &outputFileStream]() {
        processGraphs(outputFileStream);
    });
}

} /* namespace workflow */
