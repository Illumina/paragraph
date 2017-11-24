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

#include <sstream>
#include <string>

#include "common.hh"
#include "grm/CountAndGenotype.hh"
#include "grm/Parameters.hh"

TEST(Grmpy, GenotypesSingleSwap)
{
    std::string input_path
        = g_testenv->getBasePath() + "/../share/test-data/genotyping_test/chr4_graph_typing.2sample.json";
    std::string reference_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4_graph_typing.fa";
    std::string manifest_path
        = g_testenv->getBasePath() + "/../share/test-data/genotyping_test/chr4_graph_typing.manifest";
    std::string output_path = std::string();
    grm::Parameters parameters(true, 0.01, 16, 40, false);
    parameters.load(input_path, reference_path, manifest_path, output_path);
    std::ostringstream out_stream;
    grm::countAndGenotype(parameters, &out_stream);
    std::string expected_output = "#Sample\tGenotype\nSAMPLE1,0/0\nSAMPLE2,0/0\n";
    EXPECT_EQ(expected_output, out_stream.str());
}