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

#include "graphs/GraphCoordinates.hh"

#include "gtest/gtest.h"

#include "common/Error.hh"

using namespace testing;
using namespace graphs;

class GraphCoordinatesTest : public Test
{
public:
    WalkableGraph graph;

    void SetUp() override
    {
        Graph init_graph;
        init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());

        init_graph.nodes[0] = std::make_shared<Node>();
        init_graph.nodes[0]->set_id(0);
        init_graph.nodes[0]->set_name("LF");
        init_graph.nodes[0]->set_sequence("AAAAAAAAAAA");

        init_graph.nodes[1] = std::make_shared<Node>();
        init_graph.nodes[1]->set_id(1);
        init_graph.nodes[1]->set_name("P1");
        init_graph.nodes[1]->set_sequence("TTTTTT");

        init_graph.nodes[2] = std::make_shared<Node>();
        init_graph.nodes[2]->set_id(2);
        init_graph.nodes[2]->set_name("Q1");
        init_graph.nodes[2]->set_sequence("GGGGGGGG");

        init_graph.nodes[3] = std::make_shared<Node>();
        init_graph.nodes[3]->set_id(3);
        init_graph.nodes[3]->set_name("RF");
        init_graph.nodes[3]->set_sequence("AAAAAAAAAAA");

        /**
         *
         * LF           RF
         *  |           ^
         *  |           |
         *  *-> P1 -----*
         *  |           |
         *  *-> Q1 -----*
         */

        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[0]->set_from(0);
        init_graph.edges[0]->set_to(1);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[1]->set_from(0);
        init_graph.edges[1]->set_to(2);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[2]->set_from(1);
        init_graph.edges[2]->set_to(3);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[3]->set_from(2);
        init_graph.edges[3]->set_to(3);

        LOG()->set_level(spdlog::level::err);

        graph = init_graph;
    }
};

TEST_F(GraphCoordinatesTest, CanonicalPositionLookup)
{
    GraphCoordinates coordinates(graph);

    // LF has offset 0
    ASSERT_EQ(static_cast<uint64_t>(6), coordinates.canonicalPos("LF", 6));
    ASSERT_EQ(static_cast<uint64_t>(11 + 4), coordinates.canonicalPos("P1", 4));
    ASSERT_EQ(static_cast<uint64_t>(11 + 6 + 3), coordinates.canonicalPos("Q1", 3));
    ASSERT_EQ(static_cast<uint64_t>(11 + 6 + 8 + 2), coordinates.canonicalPos("RF", 2));
}

TEST_F(GraphCoordinatesTest, ReverseLookup)
{
    GraphCoordinates coordinates(graph);
    for (size_t j = 0; j < graph.node("LF")->sequence().size(); ++j)
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(j, n, offset);
        ASSERT_EQ("LF", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }

    for (size_t j = 0; j < graph.node("P1")->sequence().size(); ++j)
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + j, n, offset);
        ASSERT_EQ("P1", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
    for (size_t j = 0; j < graph.node("Q1")->sequence().size(); ++j)
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + 6 + j, n, offset);
        ASSERT_EQ("Q1", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
    for (size_t j = 0; j < graph.node("RF")->sequence().size(); ++j)
    {
        std::string n;
        uint64_t offset = static_cast<uint64_t>(-1);
        coordinates.nodeAndOffset(11 + 6 + 8 + j, n, offset);
        ASSERT_EQ("RF", n);
        ASSERT_EQ(static_cast<uint64_t>(j), offset);
    }
}

TEST_F(GraphCoordinatesTest, DistanceComputation)
{
    GraphCoordinates coordinates(graph);

    // both on LF
    ASSERT_EQ(static_cast<uint64_t>(5), coordinates.distance(10, 5));
    ASSERT_EQ(static_cast<uint64_t>(5), coordinates.distance(5, 10));

    // one on LF, one on neighbour (P1 or Q1)
    ASSERT_EQ(static_cast<uint64_t>(8), coordinates.distance(14, 6));
    ASSERT_EQ(static_cast<uint64_t>(8), coordinates.distance(20, 6));

    // LF -> RF should go via P1 because this is shorter
    ASSERT_EQ(static_cast<uint64_t>(9 + 6 + 4), coordinates.distance(2, 11 + 6 + 8 + 4));
}
