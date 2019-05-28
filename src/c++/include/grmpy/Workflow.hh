// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Paragraph
// Copyright (c) 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//		http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
// See the License for the specific language governing permissions and limitations
//
//

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
