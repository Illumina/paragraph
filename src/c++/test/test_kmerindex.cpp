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

#include "graphs/KmerIndex.hh"

#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

#include "graphs/Graph.hh"
#include "graphs/GraphBuilders.hh"
#include "graphs/GraphPath.hh"

using std::cerr;
using std::endl;
using std::list;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using namespace testing;
using namespace graphs;

class TinyDeletionGraph : public ::testing::Test
{
public:
    void SetUp()
    {
        makeDeletionGraph(left_flank, deletion, right_flank, graph);
        wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    }

    const string left_flank = "AC";
    const string deletion = "GG";
    const string right_flank = "CAG";
    Graph graph;
    std::shared_ptr<WalkableGraph> wgraph_ptr;
};

class RepetitiveDoubleSwapGraph : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        left_flank = "AAA";
        deletion1 = "TTT";
        insertion1 = "CCC";
        mid = "AAA";
        deletion2 = "TTT";
        insertion2 = "AAA";
        right_flank = "TTT";
        makeDoubleSwapGraph(left_flank, deletion1, insertion1, mid, deletion2, insertion2, right_flank, graph);
        wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    }
    Graph graph;
    string left_flank;
    string deletion1;
    string insertion1;
    string mid;
    string deletion2;
    string insertion2;
    string right_flank;

    std::shared_ptr<WalkableGraph> wgraph_ptr;
};

TEST_F(TinyDeletionGraph, InitializeKmerIndexWith1mers)
{
    const int32_t kmer_size = 1;
    KmerIndex kmer_index(wgraph_ptr, kmer_size);

    const list<GraphPath> a_paths = { GraphPath(wgraph_ptr, 0, { 0 }, 0), GraphPath(wgraph_ptr, 1, { 2 }, 1) };
    const list<GraphPath> c_paths = { GraphPath(wgraph_ptr, 1, { 0 }, 1), GraphPath(wgraph_ptr, 0, { 2 }, 0) };
    const list<GraphPath> g_paths = { GraphPath(wgraph_ptr, 0, { 1 }, 0), GraphPath(wgraph_ptr, 1, { 1 }, 1),
                                      GraphPath(wgraph_ptr, 2, { 2 }, 2) };

    const StringToPathsMap kmer_to_paths_maps = { { "A", a_paths }, { "C", c_paths }, { "G", g_paths } };

    KmerIndex expected_kmer_index(kmer_to_paths_maps);
    ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST_F(TinyDeletionGraph, InitializeKmerIndexWith2mers)
{
    const int32_t kmer_size = 2;
    KmerIndex kmer_index(wgraph_ptr, kmer_size);

    const list<GraphPath> ac_paths = { GraphPath(wgraph_ptr, 0, { 0 }, 1) };
    const list<GraphPath> cg_paths = { GraphPath(wgraph_ptr, 1, { 0, 1 }, 0) };
    const list<GraphPath> cc_paths = { GraphPath(wgraph_ptr, 1, { 0, 2 }, 0) };
    const list<GraphPath> gg_paths = { GraphPath(wgraph_ptr, 0, { 1 }, 1) };
    const list<GraphPath> gc_paths = { GraphPath(wgraph_ptr, 1, { 1, 2 }, 0) };
    const list<GraphPath> ca_paths = { GraphPath(wgraph_ptr, 0, { 2 }, 1) };
    const list<GraphPath> ag_paths = { GraphPath(wgraph_ptr, 1, { 2 }, 2) };

    const StringToPathsMap kmer_to_paths_maps
        = { { "AC", ac_paths }, { "CG", cg_paths }, { "CC", cc_paths }, { "GG", gg_paths },
            { "GC", gc_paths }, { "CA", ca_paths }, { "AG", ag_paths } };

    KmerIndex expected_kmer_index(kmer_to_paths_maps);
    ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST_F(TinyDeletionGraph, KmerIndexReportsKmersWithNonzeroCount)
{
    const int32_t kmer_size = 2;
    KmerIndex kmer_index(wgraph_ptr, kmer_size);
    const unordered_set<string> expected_kmers = { "AC", "CG", "CC", "GG", "GC", "CA", "AG" };
    ASSERT_EQ(expected_kmers, kmer_index.getKmersWithNonzeroCount());
}

TEST_F(RepetitiveDoubleSwapGraph, ExtractPathsContainingKmer)
{
    const int32_t kmer_size = 4;
    KmerIndex kmer_index(wgraph_ptr, kmer_size);
    const list<GraphPath> paths = kmer_index.getPaths("AATT");
    const list<GraphPath> expected_paths
        = { GraphPath(wgraph_ptr, 1, { 0, 1 }, 1), GraphPath(wgraph_ptr, 1, { 3, 4 }, 1),
            GraphPath(wgraph_ptr, 1, { 5, 6 }, 1) };
    ASSERT_EQ(expected_paths, paths);
}

TEST_F(RepetitiveDoubleSwapGraph, CheckIfIndexContainsKmer)
{
    const int32_t kmer_size = 6;
    KmerIndex kmer_index(wgraph_ptr, kmer_size);
    EXPECT_TRUE(kmer_index.contains("AAATTT"));
    EXPECT_FALSE(kmer_index.contains("AAATTG"));
    EXPECT_FALSE(kmer_index.contains("AAA"));
}

TEST_F(RepetitiveDoubleSwapGraph, GetNumberOfPathsContainingKmer)
{
    {
        const int32_t kmer_size = 6;
        KmerIndex kmer_index(wgraph_ptr, kmer_size);
        EXPECT_EQ(3u, kmer_index.numPaths("AAATTT"));
        EXPECT_EQ(0u, kmer_index.numPaths("AAATTG"));
        EXPECT_EQ(1u, kmer_index.numPaths("TTTTTT"));
    }
    {
        const int32_t kmer_size = 1;
        KmerIndex kmer_index(wgraph_ptr, kmer_size);
        EXPECT_EQ(9u, kmer_index.numPaths("A"));
        EXPECT_EQ(3u, kmer_index.numPaths("C"));
        EXPECT_EQ(9u, kmer_index.numPaths("T"));
        EXPECT_EQ(0u, kmer_index.numPaths("G"));
    }
}
