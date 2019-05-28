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
