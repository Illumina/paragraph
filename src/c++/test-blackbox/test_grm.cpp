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
