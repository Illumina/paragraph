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
 * \brief Workflow implementation for grmpy
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
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Error.hh"
#include "common/JsonHelpers.hh"
#include "common/Threads.hh"
#include "grmpy/AlignSamples.hh"
#include "grmpy/CountAndGenotype.hh"
#include "grmpy/Workflow.hh"
#include "paragraph/Disambiguation.hh"

namespace grmpy
{

Workflow::Workflow(
    const std::vector<std::string>& graphSpecPaths, const std::string& genotypingParameterPath,
    const genotyping::Samples& mainfest, const std::string& outputFilePath, const std::string& outputFolderPath,
    bool gzipOutput, const Parameters& parameters, const std::string& referencePath, bool progress)
    : graphSpecPaths_(graphSpecPaths)
    , genotypingParameterPath_(genotypingParameterPath)
    , manifest_(mainfest)
    , outputFilePath_(outputFilePath)
    , outputFolderPath_(outputFolderPath)
    , gzipOutput_(gzipOutput)
    , parameters_(parameters)
    , referencePath_(referencePath)
    , progress_(progress)
{
    alignedSamples_.resize(std::max<std::size_t>(1, graphSpecPaths_.size()));
    for (const genotyping::SampleInfo& sample : manifest_)
    {
        unalignedSamples_.push_back(UnalignedSample(
            sample, sample.get_alignment_data().isNull() ? std::begin(graphSpecPaths_) : std::end(graphSpecPaths_)));
        if (graphSpecPaths_.size())
        {
            for (std::size_t i = 0; i < graphSpecPaths.size(); ++i)
            {
                alignedSamples_[i].push_back(sample);
            }
        }
        else
        {
            // no graphs given. Assume all samples are pre-aligned
            assert(!sample.get_alignment_data().isNull());
            alignedSamples_[0].push_back(sample);
        }
    }
}

void Workflow::makeOutputFile(const Json::Value& output, const std::string& graphSpecPath) const
{
    const boost::filesystem::path inputPath(graphSpecPath);
    boost::filesystem::path outputPath = boost::filesystem::path(outputFolderPath_) / inputPath.filename();
    if (gzipOutput_)
    {
        outputPath += ".gz";
    }
    boost::iostreams::basic_file_sink<char> of(outputPath.string());
    if (!of.is_open())
    {
        error("ERROR: Failed to open output file '%s'. Error: '%s'", outputPath.string().c_str(), std::strerror(errno));
    }

    boost::iostreams::filtering_ostream fos;
    if (gzipOutput_)
    {
        fos.push(boost::iostreams::gzip_compressor());
    }
    fos.push(of);
    fos << common::writeJson(output);
}

void Workflow::alignSamples()
{
    std::lock_guard<std::mutex> lock(mutex_);
    for (std::size_t i = 0; i < unalignedSamples_.size(); ++i)
    {
        UnalignedSample& input = unalignedSamples_[i];
        if (progress_)
        {
            LOG()->critical(
                "Starting alignment for sample {} ({}/{})", input.sample_.sample_name(), i + 1,
                unalignedSamples_.size());
        }
        common::BamReader reader(input.sample_.filename(), input.sample_.index_filename(), referencePath_);
        while (graphSpecPaths_.end() != input.unprocessedGraphs_)
        {
            const GraphSpecPaths::const_iterator ourGraph = input.unprocessedGraphs_++;
            ASYNC_BLOCK_WITH_CLEANUP([this](bool failure) { terminate_ |= failure; })
            {
                if (terminate_)
                {
                    LOG()->warn("terminating");
                    break;
                }
                common::unlock_guard<std::mutex> unlock(mutex_);

                auto ourGraphIndex = static_cast<unsigned long>(std::distance(graphSpecPaths_.begin(), ourGraph));
                alignSingleSample(
                    parameters_, *ourGraph, referencePath_, reader, alignedSamples_.at(ourGraphIndex).at(i));

                if (progress_)
                {
                    LOG()->critical(
                        "Sample {}: Alignment {} / {} finished", input.sample_.sample_name(), ourGraphIndex + 1,
                        alignedSamples_.size());
                }
            }
        }
    }
}

void Workflow::genotypeGraphs(
    std::ostream& outputFileStream, std::vector<genotyping::Samples>::const_iterator& ungenotypedSamples)
{
    std::lock_guard<std::mutex> lock(mutex_);
    while (alignedSamples_.end() != ungenotypedSamples)
    {
        const std::vector<genotyping::Samples>::const_iterator ourGraphSamples = ungenotypedSamples++;
        const auto graphIndex = static_cast<unsigned long>(std::distance(alignedSamples_.cbegin(), ourGraphSamples));
        Json::Value output;

        ASYNC_BLOCK_WITH_CLEANUP([this](bool failure) { terminate_ |= failure; })
        {
            if (terminate_)
            {
                LOG()->warn("terminating");
                break;
            }
            common::unlock_guard<std::mutex> unlock(mutex_);

            const std::string& graphSpecPath = graphSpecPaths_.empty() ? std::string() : graphSpecPaths_.at(graphIndex);
            output = countAndGenotype(graphSpecPath, referencePath_, genotypingParameterPath_, *ourGraphSamples);
            if (!outputFolderPath_.empty())
            {
                makeOutputFile(output, graphSpecPath);
            }
            if (progress_)
            {
                LOG()->critical("Genotyping finished for graph {} / {}", graphIndex + 1, alignedSamples_.size());
            }
        }

        if (!outputFilePath_.empty())
        {
            if (firstPrinted_)
            {
                outputFileStream << ',';
            }
            outputFileStream << common::writeJson(output);
            firstPrinted_ = true;
        }
    }
}

void Workflow::run()
{
    boost::iostreams::filtering_ostream fos;
    if (gzipOutput_)
    {
        fos.push(boost::iostreams::gzip_compressor());
    }
    if (!outputFilePath_.empty())
    {
        if ("-" != outputFilePath_)
        {
            LOG()->info("Output file path: {}", outputFilePath_);
            boost::iostreams::basic_file_sink<char> of(outputFilePath_);
            if (!of.is_open())
            {
                error(
                    "ERROR: Failed to open output file '%s'. Error: '%s'", outputFilePath_.c_str(),
                    std::strerror(errno));
            }
            fos.push(of);
        }
        else
        {
            LOG()->info("Output to stdout");
            fos.push(std::cout);
        }
    }

    if (!outputFilePath_.empty() && 1 < graphSpecPaths_.size())
    {
        fos << "[";
    }

    LOG()->info("Aligning for {} graphs", graphSpecPaths_.size());
    common::CPU_THREADS(parameters_.threads()).execute([this]() { alignSamples(); });

    LOG()->info("Genotyping {} samples", alignedSamples_.size());
    std::vector<genotyping::Samples>::const_iterator ungenotypedSamples = alignedSamples_.begin();
    common::CPU_THREADS(parameters_.threads()).execute([this, &fos, &ungenotypedSamples]() {
        genotypeGraphs(fos, ungenotypedSamples);
    });

    if (!outputFilePath_.empty() && 1 < graphSpecPaths_.size())
    {
        fos << "]\n";
    }
}

} /* namespace grmpy */
