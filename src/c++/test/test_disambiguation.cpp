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

#include "common/Threads.hh"
#include "grm/Align.hh"
#include "grm/GraphAligner.hh"
#include "grm/GraphInput.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/GraphVariants.hh"

#include "gtest/gtest.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;

using namespace testing;
using namespace common;
using namespace grm;
using namespace graphtools;

class DisambiguationTest : public Test
{
public:
    ReadBuffer reads;
    Graph graph;

    virtual void SetUp()
    {
        reads.emplace_back(
            new Read("f0", "AAAAAAAAAATTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "AAAAAAAAAATTTTTTTTTTTTTTTTTTTTAAAAAAAAAA"));
        reads.emplace_back(new Read("f1", "AAAAAAAAAATTTTTTTTTTT", "AAAAAAAAAATTTTTTTTTTT"));
        reads.emplace_back(
            new Read("f2", "AAAAAAAAAATTTTTTTTTTGGGGGGGGGGAAAAAAAAAA", "AAAAAAAAAATTTTTTTTTTGGGGGGGGGGAAAAAAAAAA"));
        reads.emplace_back(new Read("f3", "AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA"));

        graph = Graph(5);
        graph.setNodeName(0, "LF");
        graph.setNodeSeq(0, "AAAAAAAAAA");
        graph.setNodeName(1, "R1");
        graph.setNodeSeq(1, "TTTTTTTTTT");
        graph.setNodeName(2, "R2");
        graph.setNodeSeq(2, "TTTTTTTTTT");
        graph.setNodeName(3, "A1");
        graph.setNodeSeq(3, "GGGGGGGGGG");
        graph.setNodeName(4, "RF");
        graph.setNodeSeq(4, "AAAAAAAAAA");
        /**
         * LF--R1--R2-->RF
         *  |  |        |
         *  |  >--A1----^
         *  |           |
         *  >-----------^
         */
        graph.addEdge(0, 1);
        graph.addEdge(0, 4);
        graph.addEdge(1, 2);
        graph.addEdge(1, 3);
        graph.addEdge(2, 4);
        graph.addEdge(3, 4);
        graph.addLabelToEdge(0, 1, "R");
        graph.addLabelToEdge(1, 2, "R");
        graph.addLabelToEdge(2, 4, "R");
        graph.addLabelToEdge(0, 4, "D");

        LOG()->set_level(spdlog::level::err);
        // this ensures results will be in predictable order
        common::CPU_THREADS().reset(1);
        std::list<Path> paths;
        grm::alignReads(&graph, paths, reads, nullptr, false, true, false, false, false);
        paragraph::disambiguateReads(&graph, reads);
    }
};

TEST_F(DisambiguationTest, DisambiguatesReads)
{
    ASSERT_EQ(1ull, reads[0]->graph_sequences_supported().size());
    ASSERT_EQ("R", reads[0]->graph_sequences_supported(0));
    ASSERT_EQ(1ull, reads[1]->graph_sequences_supported().size());
    ASSERT_EQ("R", reads[1]->graph_sequences_supported(0));
    ASSERT_EQ(0ull, reads[2]->graph_sequences_supported().size());
    ASSERT_EQ(1ull, reads[3]->graph_sequences_supported().size());
    ASSERT_EQ("D", reads[3]->graph_sequences_supported(0));
}