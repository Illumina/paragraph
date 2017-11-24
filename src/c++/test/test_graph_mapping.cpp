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

#include "graphs/GraphMapping.hh"

#include "gtest/gtest.h"

#include "graphs/Graph.hh"
#include "graphs/GraphBuilders.hh"

using std::string;
using std::vector;

using namespace graphs;

TEST(Operation, InitializesFromString)
{
    Operation operation("3M", "ATC", "ATC");
    Operation expected_operation('M', 3, "ATC", "ATC");
    ASSERT_EQ(expected_operation, operation);
}

TEST(Operation, OutputsQueryAndReferenceSpans)
{
    {
        Operation operation("3M", "AAA", "AAA");
        EXPECT_EQ((int32_t)3, operation.querySpan());
        EXPECT_EQ((int32_t)3, operation.referenceSpan());
    }
    {
        Operation operation("4X", "AAAA", "TTTT");
        EXPECT_EQ((int32_t)4, operation.querySpan());
        EXPECT_EQ((int32_t)4, operation.referenceSpan());
    }
    {
        Operation operation("5D", "", "AAAAA");
        EXPECT_EQ((int32_t)0, operation.querySpan());
        EXPECT_EQ((int32_t)5, operation.referenceSpan());
    }
    {
        Operation operation("7I", "AAAAAAA", "");
        EXPECT_EQ((int32_t)7, operation.querySpan());
        EXPECT_EQ((int32_t)0, operation.referenceSpan());
    }
    {
        Operation operation("10S", "AAAAAAAAAA", "");
        EXPECT_EQ((int32_t)10, operation.querySpan());
        EXPECT_EQ((int32_t)0, operation.referenceSpan());
    }
    {
        Operation operation("7N", "NNNNNNN", "NNNNNNN");
        EXPECT_EQ((int32_t)7, operation.querySpan());
        EXPECT_EQ((int32_t)7, operation.referenceSpan());
    }
    {
        Operation operation("3N", "NCN", "CNN");
        EXPECT_EQ((int32_t)3, operation.querySpan());
        EXPECT_EQ((int32_t)3, operation.referenceSpan());
    }
}

TEST(Operation, ErrorsOutOnUnexpectedSequences)
{
    EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATCG"));
    EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATC"));
    EXPECT_ANY_THROW(Operation("4M", "AAA", "AAA"));

    EXPECT_ANY_THROW(Operation("4N", "NNN", "NNN"));
    EXPECT_ANY_THROW(Operation("3N", "NN", "NNN"));
    EXPECT_ANY_THROW(Operation("2N", "NT", "NT"));

    EXPECT_ANY_THROW(Operation("2X", "AT", "TT"));
    EXPECT_ANY_THROW(Operation("2X", "AT", "A"));

    EXPECT_ANY_THROW(Operation("4D", "AAA", ""));
    EXPECT_ANY_THROW(Operation("4D", "", ""));

    EXPECT_ANY_THROW(Operation("2I", "AA", "T"));

    EXPECT_ANY_THROW(Operation("2S", "TTT", ""));
    EXPECT_ANY_THROW(Operation("2S", "TT", "T"));
}

TEST(Mapping, InitializesFromCigar)
{
    // query: ---TTCGTT--TTGGGTCCCCCCCCCC
    //           ||| ||  ||   |
    //   ref: CCCTTCCNNAATT---T----------
    string query = "TTCGTTTTGGGTCCCCCCCCCC";
    string reference = "CCCTTCCNNAATTT";

    Mapping mapping(3, "3M1X2N2D2M3I1M10S", query, reference);
    const vector<Operation> operations
        = { Operation('M', 3, "TTC", "TTC"), Operation('X', 1, "G", "C"),         Operation('N', 2, "TT", "NN"),
            Operation('D', 2, "", "AA"),     Operation('M', 2, "TT", "TT"),       Operation('I', 3, "GGG", ""),
            Operation('M', 1, "T", "T"),     Operation('S', 10, "CCCCCCCCCC", "") };

    Mapping expected_mapping(3, operations, query, reference);
    ASSERT_EQ(expected_mapping, mapping);
}

TEST(Mapping, CalculatesQueryAndReferenceSpans)
{
    Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC", "CCCTTCCTTAATTT");
    EXPECT_EQ((int32_t)22, mapping.querySpan());
    EXPECT_EQ((int32_t)11, mapping.referenceSpan());
}

TEST(Mapping, StitchesQueryAndReferenceSequences)
{
    Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC", "CCCTTCCTTAATTT");
    EXPECT_EQ("TTCGTTTTGGGT", mapping.query());
    EXPECT_EQ("TTCCTTAATTT", mapping.reference());
}

TEST(NodeMapping, InitializesFromString)
{
    Graph graph;
    makeDeletionGraph("CCCC", "TTTAAAAGGGG", "TTTT", graph);
    const string query = "AAAATTTTT";
    NodeMapping node_mapping(3, "1[4M5S]", query, graph);
    EXPECT_EQ(1ull, node_mapping.node_id());
    EXPECT_EQ(Mapping(0, "4M5S", query, "AAAAGGGG"), node_mapping);
}

TEST(GraphMapping, InitializesFromString)
{
    Graph graph;
    makeDeletionGraph("AAAA", "TTGG", "TTTT", graph);
    const string query = "AAAATTCCC";
    GraphMapping graph_mapping(0, "0[4M]1[2M3S]", query, graph);
    ASSERT_EQ((size_t)2, graph_mapping.size());
    EXPECT_EQ(NodeMapping(0, "0[4M]", "AAAA", graph), graph_mapping[0]);
    EXPECT_EQ(NodeMapping(0, "1[2M3S]", "TTCCC", graph), graph_mapping[1]);
}

TEST(GraphMapping, StitchesQueryAndReferenceSequences)
{
    Graph graph;
    makeDeletionGraph("AAAA", "TTGG", "TTTT", graph);
    const string query = "AAAATTCCC";
    GraphMapping graph_mapping(0, "0[4M]1[2M3S]", query, graph);
    EXPECT_EQ("AAAATT", graph_mapping.query());
    EXPECT_EQ("AAAATT", graph_mapping.reference());
}

TEST(GraphMapping, CalculatesQueryAndReferenceSpans)
{
    Graph graph;
    makeDeletionGraph("AAAA", "TTGG", "TTTT", graph);
    const string query = "AAAATTCCC";
    GraphMapping graph_mapping(0, "0[4M]1[2M3S]", query, graph);
    EXPECT_EQ((int32_t)9, graph_mapping.querySpan());
    EXPECT_EQ((int32_t)6, graph_mapping.referenceSpan());
}