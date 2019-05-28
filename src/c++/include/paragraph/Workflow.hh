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
