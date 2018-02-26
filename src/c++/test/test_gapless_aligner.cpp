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

#include "graphs/GaplessAligner.hh"

#include <list>
#include <string>

#include "gtest/gtest.h"

#include "graphs/GraphBuilders.hh"
#include "graphs/GraphMapping.hh"
#include "graphs/GraphMappingOperations.hh"
#include "graphs/GraphPath.hh"
#include "graphs/WalkableGraph.hh"

using std::list;
using std::string;

using namespace graphs;

TEST(SequenceAlignment, ReferenceShorterThanQuery_CausesError) { EXPECT_ANY_THROW(alignWithoutGaps("AAAA", 0, "AAA")); }

TEST(SequenceAlignment, EmptySequences_CauseError) { EXPECT_ANY_THROW(alignWithoutGaps("", 0, "")); }

TEST(SequenceAlignment, TypicalSequences_Aligned)
{
    const string query = "AGGTTTTG";
    const string reference = "NNNNATCGTTTG";
    const Mapping expected_mapping(4, "1M3X4M", query, reference);

    ASSERT_EQ(expected_mapping, alignWithoutGaps(query, 4, reference));
}

TEST(AlignmentOfSequenceToPath, SingleNodePath_Aligned)
{
    Graph graph;
    makeDeletionGraph("AAAACC", "TTTGG", "ATTT", graph);
    std::shared_ptr<WalkableGraph> wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    GraphPath path(wgraph_ptr, 1, { 1 }, 4);
    const string sequence = "ATGC";

    GraphMapping expected_graph_mapping = decodeFromString(1, "1[1X2M1X]", sequence, graph);
    GraphMapping graph_mapping = alignWithoutGaps(path, sequence);
    EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST(AlignmentOfSequenceToPath, MultiNodePath_Aligned)
{
    Graph graph;
    makeDeletionGraph("AAAACC", "TTTGG", "ATTT", graph);
    std::shared_ptr<WalkableGraph> wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    GraphPath path(wgraph_ptr, 2, { 0, 1, 2 }, 1);
    const string sequence = "TTCCTTAGGAT";

    GraphMapping expected_graph_mapping = decodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", sequence, graph);
    EXPECT_EQ(expected_graph_mapping, alignWithoutGaps(path, sequence));
}

TEST(KmerExtraction, TypicalSequence_KmersExtracted)
{
    const string sequence = "AAATTT";
    const list<string> expected_4mers = { "AAAT", "AATT", "ATTT" };
    ASSERT_EQ(expected_4mers, extractKmersFromAllPositions(sequence, 4));

    const list<string> expected_7mers = {};
    ASSERT_EQ(expected_7mers, extractKmersFromAllPositions(sequence, 7));
}

TEST(AlignmentOfSequenceToShortPath, TypicalSequence_BestAlignmentObtained)
{
    Graph graph;
    makeDeletionGraph("AAACC", "TTGGG", "TTAAA", graph);
    std::shared_ptr<WalkableGraph> wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    const GraphPath path(wgraph_ptr, 4, { 0 }, 4);
    const string sequence = "CCTTA";

    GraphMapping mapping = getBestAlignmentToShortPath(path, 1, sequence);

    GraphMapping expected_mapping = decodeFromString(3, "0[2M]2[3M]", sequence, graph);
    ASSERT_EQ(expected_mapping, mapping);
}

TEST(AlignmentOfSequenceToGraph, TypicalSequence_BestAlignmentObtained)
{
    Graph graph;
    makeDeletionGraph("AAAACC", "TTTGG", "ATTT", graph);
    std::shared_ptr<WalkableGraph> wgraph_ptr = std::make_shared<WalkableGraph>(graph);

    const int32_t kmer_len = 3;
    GaplessAligner aligner(wgraph_ptr, kmer_len);
    const string sequence = "TTCCTTAGGAT";

    GraphMapping mapping = aligner.getBestAlignment(sequence);

    GraphMapping expected_mapping = decodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", sequence, graph);
    ASSERT_EQ(expected_mapping, mapping);
}