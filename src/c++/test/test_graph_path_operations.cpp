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

#include "graphs/GraphPathOperations.hh"

#include "gtest/gtest.h"

#include <list>
#include <string>
#include <vector>

#include "graphs/GraphBuilders.hh"
#include "graphs/GraphPath.hh"
#include "graphs/WalkableGraph.hh"

using std::list;
using std::string;
using std::vector;

using namespace graphs;

class DeletionGraphForOperations : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        left_flank = "AAAACC";
        deletion = "TTTGG";
        right_flank = "ATTT";
        makeDeletionGraph(left_flank, deletion, right_flank, graph);
        wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    }
    Graph graph;
    string left_flank;
    string deletion;
    string right_flank;
    std::shared_ptr<WalkableGraph> wgraph_ptr;
};

TEST_F(DeletionGraphForOperations, SplittingSequenceByPathOfDifferentLengthCausesError)
{
    GraphPath path(wgraph_ptr, 3, { 0, 1 }, 2);
    const string sequence = "AA";
    EXPECT_ANY_THROW(splitByPath(path, sequence));
}

TEST_F(DeletionGraphForOperations, SplittingSequenceBySingleNodePath)
{
    GraphPath path(wgraph_ptr, 1, { 1 }, 3);
    const string sequence = "AAT";
    const vector<string> expected_pieces = { sequence };
    EXPECT_EQ(expected_pieces, splitByPath(path, sequence));
}

TEST_F(DeletionGraphForOperations, SplittingSequenceByMultiNodePath)
{
    {
        GraphPath path(wgraph_ptr, 1, { 0, 1 }, 3);
        const string sequence = "AAAAAGGGG";
        const vector<string> expected_pieces = { "AAAAA", "GGGG" };
        EXPECT_EQ(expected_pieces, splitByPath(path, sequence));
    }
    {
        GraphPath path(wgraph_ptr, 3, { 0, 2 }, 1);
        const string sequence = "AAACC";
        const vector<string> expected_pieces = { "AAA", "CC" };
        EXPECT_EQ(expected_pieces, splitByPath(path, sequence));
    }
    {
        GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
        const string sequence = "AAAGGGGGCC";
        const vector<string> expected_pieces = { "AAA", "GGGGG", "CC" };
        EXPECT_EQ(expected_pieces, splitByPath(path, sequence));
    }
}

TEST(GraphPathOperations, GraphPathsOverlapDetected)
{
    Graph swap;
    makeSimpleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", swap);
    std::shared_ptr<WalkableGraph> wgraph_ptr(new WalkableGraph(swap));

    {
        const GraphPath p1(wgraph_ptr, 0, { 0, 1 }, 3);
        const GraphPath p2(wgraph_ptr, 0, { 1, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const GraphPath expected_merge(wgraph_ptr, 0, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }

    {
        const GraphPath p1(wgraph_ptr, 2, { 0, 1, 3 }, 2);
        const GraphPath p2(wgraph_ptr, 0, { 1, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const GraphPath expected_merge(wgraph_ptr, 2, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }
    {
        const GraphPath p1(wgraph_ptr, 2, { 0, 2 }, 1);
        const GraphPath p2(wgraph_ptr, 1, { 2 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const GraphPath expected_merge(wgraph_ptr, 2, { 0, 2 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }
}

TEST(GraphPathOperations, GraphPathsNoOverlapDetected)
{
    Graph swap;
    makeSimpleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", swap);
    std::shared_ptr<WalkableGraph> wgraph_ptr(new WalkableGraph(swap));

    {
        // p1 ends before p2 begins
        const GraphPath p1(wgraph_ptr, 0, { 0, 1 }, 1);
        const GraphPath p2(wgraph_ptr, 2, { 1, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }

    {
        // no shared nodes
        const GraphPath p1(wgraph_ptr, 0, { 0 }, 3);
        const GraphPath p2(wgraph_ptr, 2, { 1, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }

    {
        // incompatible
        const GraphPath p1(wgraph_ptr, 0, { 0, 1, 3 }, 3);
        const GraphPath p2(wgraph_ptr, 2, { 0, 2, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }
    {
        // incompatible 2
        const GraphPath p1(wgraph_ptr, 0, { 0, 1 }, 3);
        const GraphPath p2(wgraph_ptr, 2, { 2, 3 }, 3);

        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }
}

TEST(GraphPathOperations, GraphsMergedExhaustively)
{
    Graph swap;
    makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA", swap);
    std::shared_ptr<WalkableGraph> wgraph_ptr(new WalkableGraph(swap));

    {
        const GraphPath p0(wgraph_ptr, 0, { 1, 3 }, 3);
        const GraphPath p1(wgraph_ptr, 0, { 2, 3 }, 3);
        const GraphPath p2(wgraph_ptr, 0, { 3, 4 }, 3);
        const GraphPath p3(wgraph_ptr, 0, { 3, 5 }, 3);

        ASSERT_TRUE(p0.isValid());
        ASSERT_TRUE(p1.isValid());
        ASSERT_TRUE(p2.isValid());
        ASSERT_TRUE(p3.isValid());

        const std::list<GraphPath> expected_merge{
            GraphPath(wgraph_ptr, 0, { 1, 3, 4 }, 3),
            GraphPath(wgraph_ptr, 0, { 2, 3, 5 }, 3),
            GraphPath(wgraph_ptr, 0, { 2, 3, 4 }, 3),
            GraphPath(wgraph_ptr, 0, { 1, 3, 5 }, 3),
        };
        std::list<GraphPath> merged{ p0, p1, p2, p3 };
        graphs::exhaustiveMerge(merged);

        /*
        for (const auto & p : merged)
        {
            std::cerr << p << std::endl;
        }
        */
        ASSERT_EQ(merged.size(), expected_merge.size());

        for (auto it = merged.cbegin(), ex_it = expected_merge.cbegin(); it != merged.cend(); ++it, ++ex_it)
        {
            ASSERT_EQ(*it, *ex_it);
        }
    }
}
