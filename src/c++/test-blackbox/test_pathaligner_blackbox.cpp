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

#include "grm/PathAligner.hh"

#include "grm/GraphInput.hh"

#include <list>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "common.hh"
#include "common/Error.hh"
#include "common/JsonHelpers.hh"

using common::Read;
using graphtools::Graph;
using graphtools::Path;
using grm::PathAligner;
using grm::graphFromJson;

using std::list;
using std::string;

using namespace testing;
using namespace common;

TEST(PathAligner, Aligns_PGHetInsReads)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/pg-het-ins/pg-het-ins.json";

    const Json::Value input_graph = common::getJSON(graph_spec_path);
    const Graph graph = graphFromJson(input_graph, g_testenv->getHG38Path());
    const auto paths = grm::pathsFromJson(&graph, input_graph["paths"]);

    PathAligner aligner;
    aligner.setGraph(&graph, paths);

    {
        Read r;
        // clang-format off
        r.setCoreInfo(
            "HSX250:46:HFCCJCCXX:2:2203:14742:58321",
            "GCCGGCTTGTGTCAGCACTGAGCGAGGCCAGCACCTTTGAGGACCC"
            "TCAGCGCCTCTACCACCTGGGCCTCCCCAGCCACGGTGAGGACCCA"
            "CCCTGGCATGATCCCCCTCATCACCTCCCCAGCCACGGTGAGGACC"
            "CACCCTGGCATGA",
            "AAFFFKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
            "KKKKKKFKFKKKKKKKKKKFFKKKKKKKKFFKKKKKFKKKKKKKKK"
            "KKFKFKKAKKKKKKKKKKK7AK<FKKKKKKKKKKKKFK<K,A<AAK"
            "FFKF<<FFAA<A,");
        // clang-format on
        aligner.alignRead(r);
    }
}
