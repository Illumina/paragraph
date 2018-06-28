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
 * \brief Workflow declaration for grmpy
 *
 * \file Workflow.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <mutex>

#include "common/ReadExtraction.hh"
#include "grmpy/Parameters.hh"

namespace grmpy
{

class Workflow
{
    typedef std::vector<std::string> GraphSpecPaths;
    typedef std::vector<std::string> InputPaths;
    struct UnalignedSample
    {
        UnalignedSample(const genotyping::SampleInfo& sample, GraphSpecPaths::const_iterator unprocessedGraphs)
            : sample_(sample)
            , unprocessedGraphs_(unprocessedGraphs)
        {
        }
        const genotyping::SampleInfo& sample_;
        GraphSpecPaths::const_iterator unprocessedGraphs_;
    };
    const GraphSpecPaths& graphSpecPaths_;
    const std::string genotypingParameterPath_;
    const genotyping::Samples& manifest_;
    const std::string outputFilePath_;
    const std::string outputFolderPath_;
    const bool gzipOutput_;
    const Parameters& parameters_;
    const std::string referencePath_;
    std::vector<UnalignedSample> unalignedSamples_;
    // [graphs][samples]
    std::vector<genotyping::Samples> alignedSamples_;

    mutable std::mutex mutex_;
    bool terminate_ = false;

    bool firstPrinted_ = false;

    bool progress_ = true;

    void genotypeGraphs(
        std::ostream& outputFileStream, std::vector<genotyping::Samples>::const_iterator& ungenotypedSamples);
    void alignSamples();
    void makeOutputFile(const Json::Value& output, const std::string& graphSpecPath) const;

public:
    Workflow(
        const std::vector<std::string>& graphSpecPaths, const std::string& genotypingParameterPath,
        const genotyping::Samples& mainfest, const std::string& outputFilePath, const std::string& outputFolderPath,
        bool gzipOutput, const Parameters& parameters, const std::string& referencePath, bool progress);
    void run();
};

} /* namespace grmpy */
