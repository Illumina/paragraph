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

#include "graphs/GraphMappingOperations.hh"

#include <string>

#include "gtest/gtest.h"

#include "graphs/Graph.hh"
#include "graphs/GraphBuilders.hh"

using std::string;
using namespace graphs;

TEST(SplitNodeCigar, ExtractsCigarAndNodeId)
{
    const string node_cigar = "1[4M5S]";
    string cigar;
    int32_t node_id;
    splitNodeCigar(node_cigar, cigar, node_id);
    EXPECT_EQ(1, node_id);
    EXPECT_EQ("4M5S", cigar);
}

TEST(DecodeGraphMapping, DecodesTypicalGraphMappings)
{
    Graph graph;
    makeDeletionGraph("AAAA", "TTGG", "TTTT", graph);
    const string read = "AAAATTCCC";
    GraphMapping graph_mapping = decodeFromString(0, "0[4M]1[2M3S]", read, graph);

    GraphMapping expected_graph_mapping(
        { 0, 1 }, { Mapping(0, "4M", "AAAA", "AAAA"), Mapping(0, "2M3S", "TTCCC", "TTGG") });
    EXPECT_EQ(expected_graph_mapping, graph_mapping);
}