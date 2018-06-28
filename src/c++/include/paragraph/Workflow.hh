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
 * \brief Workflow declaration
 *
 * \file Workflow.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include "common/ReadExtraction.hh"
#include "paragraph/Parameters.hh"

namespace paragraph
{

class Workflow
{
    typedef std::vector<std::string> GraphSpecPaths;
    typedef std::vector<std::string> InputPaths;
    struct Input
    {
        Input(
            const InputPaths& inputPaths, const InputPaths& inputIndexPaths,
            GraphSpecPaths::const_iterator unprocessedGraphs)
            : inputPaths_(inputPaths)
            , inputIndexPaths_(inputIndexPaths)
            , unprocessedGraphs_(unprocessedGraphs)
        {
        }
        const InputPaths inputPaths_;
        const InputPaths inputIndexPaths_;
        GraphSpecPaths::const_iterator unprocessedGraphs_;
    };
    std::vector<Input> unprocessedInputs_;
    const GraphSpecPaths& graphSpecPaths_;
    const std::string& outputFilePath_;
    const std::string& outputFolderPath_;
    const bool gzipOutput_;
    const Parameters& parameters_;
    const std::string& referencePath_;
    const std::string& targetRegions_;

    mutable std::mutex mutex_;
    bool terminate_ = false;

    bool firstPrinted_ = false;

    std::string processGraph(
        const std::string& graphSpecPath, const Parameters& parameters, const InputPaths& inputPaths,
        std::vector<common::BamReader>& readers);
    void processGraphs(std::ostream& outputFileStream);
    void makeOutputFile(const std::string& output, const std::string& graphSpecPath);

public:
    Workflow(
        bool jointInputs, const std::vector<std::string>& inpuPaths, const InputPaths& inputIndexPaths,
        const std::vector<std::string>& graphSpecPaths, const std::string& outputFilePath,
        const std::string& outputFolderPath, bool gzipOutput, const Parameters& parameters,
        const std::string& referencePath, const std::string& targetRegions);
    void run();
};

} /* namespace paragraph */
