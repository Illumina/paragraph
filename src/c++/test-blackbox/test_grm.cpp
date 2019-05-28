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

#include "gtest/gtest.h"

#include <iostream>
#include <sstream>
#include <string>

#include "common.hh"
#include "common/BamReader.hh"
#include "grmpy/AlignSamples.hh"
#include "grmpy/CountAndGenotype.hh"
#include "grmpy/Parameters.hh"
#include "json/json.h"

TEST(Grmpy, GenotypesSingleSwap)
{
    std::string graph_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chrX_graph_typing.2sample.json";
    std::string reference_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chrX_graph_typing.fa";
    std::string manifest_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chrX_graph_typing.manifest";
    std::string genotype_parameter_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/param.json";
    grmpy::Parameters parameters;

    genotyping::Samples samples = genotyping::loadManifest(manifest_path);

    for (auto& sample : samples)
    {
        common::BamReader reader(sample.filename(), sample.index_filename(), reference_path);

        alignSingleSample(parameters, graph_path, reference_path, reader, sample);
    }

    const Json::Value genotype = grmpy::countAndGenotype(graph_path, reference_path, genotype_parameter_path, samples);

    std::stringstream out_stream;
    out_stream << genotype;

    std::istream* json_stream = &out_stream;
    Json::Value result;
    (*json_stream) >> result;

    EXPECT_EQ("REF", result["samples"]["SAMPLE1"]["gt"]["GT"].asString());
    EXPECT_EQ("REF/REF", result["samples"]["SAMPLE2"]["gt"]["GT"].asString());
    std::cout << out_stream.str() << std::endl;
}
