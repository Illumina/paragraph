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

#include "graphs/GraphPath.hh"

#include <list>

#include "gtest/gtest.h"

#include "graphs/Graph.hh"
#include "graphs/GraphBuilders.hh"
#include "graphs/WalkableGraph.hh"

using std::list;
using std::string;

using namespace graphs;

class DeletionGraph : public ::testing::Test
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

class SwapGraph : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        left_flank = "AAAAA";
        deletion = "CCCC";
        insertion = "TTT";
        right_flank = "GGGG";
        makeSimpleSwapGraph(left_flank, deletion, insertion, right_flank, graph);
        wgraph_ptr = std::make_shared<WalkableGraph>(graph);
    }
    Graph graph;
    string left_flank;
    string deletion;
    string insertion;
    string right_flank;

    std::shared_ptr<WalkableGraph> wgraph_ptr;
};

class DoubleSwapGraph : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        left_flank = "AAAAA";
        deletion1 = "CCCC";
        insertion1 = "TTT";
        mid = "CCCC";
        deletion2 = "GGGGG";
        insertion2 = "AAA";
        right_flank = "TTTT";
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

TEST_F(DeletionGraph, PathReturnsItsSequence)
{
    {
        GraphPath path(wgraph_ptr, 3, { 0 }, 3);
        EXPECT_EQ("A", path.seq());
    }
    {
        GraphPath path(wgraph_ptr, 3, { 1 }, 4);
        EXPECT_EQ("GG", path.seq());
    }

    {
        GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
        EXPECT_EQ("ACCTTTGGAT", path.seq());
    }
}

TEST_F(DeletionGraph, PathReturnsItsLength)
{
    {
        GraphPath path(wgraph_ptr, 3, { 0 }, 3);
        EXPECT_EQ((size_t)1, path.length());
    }
    {
        GraphPath path(wgraph_ptr, 3, { 1 }, 4);
        EXPECT_EQ((size_t)2, path.length());
    }

    {
        GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
        EXPECT_EQ((size_t)10, path.length());
    }
}

TEST_F(DeletionGraph, PathReturnsLengthOfNodeOverlaps)
{
    {
        GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
        EXPECT_EQ((size_t)3, path.lengthOnNode(0));
        EXPECT_EQ((size_t)5, path.lengthOnNode(1));
        EXPECT_EQ((size_t)2, path.lengthOnNode(2));
    }
    {
        GraphPath path(wgraph_ptr, 1, { 2 }, 2);
        EXPECT_EQ((size_t)2, path.lengthOnNode(2));
    }
}

TEST_F(DeletionGraph, PathReturnsSequenceOfNodeItOverlaps)
{
    GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
    EXPECT_EQ("ACC", path.seqOnNode(0));
    EXPECT_EQ("TTTGG", path.seqOnNode(1));
    EXPECT_EQ("AT", path.seqOnNode(2));
}

TEST_F(DeletionGraph, WellFormedPathIsValid)
{
    GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
    ASSERT_TRUE(path.isValid());
}

TEST_F(DeletionGraph, PathStartingOutsideOfNodeSequenceIsInvalid)
{

    GraphPath path(wgraph_ptr, 6, { 0, 1, 2 }, 1);
    ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathEndingOutsideOfNodeSequenceIsInvalid)
{
    GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 10);
    ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathWithUnorderedNodesIsInvalid)
{
    GraphPath path(wgraph_ptr, 3, { 2, 1 }, 1);
    ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, SingleNodePathWithEndBeforeStartIsInvalid)
{
    GraphPath path(wgraph_ptr, 3, { 0 }, 1);
    ASSERT_FALSE(path.isValid());
}

TEST_F(SwapGraph, DisconnectedPathIsInvalid)
{
    GraphPath path(wgraph_ptr, 0, { 0, 3 }, 0);
    ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathIsEncodedAsString)
{
    {
        GraphPath path(wgraph_ptr, 0, { 0 }, 1);
        ASSERT_EQ("(LF@0)-(LF@1)", path.encode());
    }
    {
        GraphPath path(wgraph_ptr, 3, { 0, 1, 2 }, 1);
        ASSERT_EQ("(LF@3)-(DEL)-(RF@1)", path.encode());
    }
}

TEST_F(DeletionGraph, PathStartPositionIsExtended)
{
    {
        GraphPath path(wgraph_ptr, 3, { 0, 1 }, 1);
        GraphPath extended_path = path.extendStartPosition(3);
        GraphPath expected_path(wgraph_ptr, 0, { 0, 1 }, 1);
        ASSERT_EQ(expected_path, extended_path);
    }
    {
        GraphPath path(wgraph_ptr, 3, { 0, 1 }, 1);
        GraphPath extended_path = path.extendEndPosition(2);
        GraphPath expected_path(wgraph_ptr, 3, { 0, 1 }, 3);
        ASSERT_EQ(expected_path, extended_path);
    }
}

TEST_F(DeletionGraph, InvalidPathExtensionCausesError)
{
    GraphPath path(wgraph_ptr, 3, { 0, 1 }, 1);
    ASSERT_ANY_THROW(path.extendStartPosition(4));
    ASSERT_ANY_THROW(path.extendEndPosition(4));
}

TEST_F(SwapGraph, PathIsExtendedToAdjacentNode)
{
    {
        GraphPath path(wgraph_ptr, 1, { 1, 3 }, 2);
        GraphPath extended_path = path.extendStartNodeTo(0);
        const int32_t new_left_pos = left_flank.length() - 1;
        GraphPath expected_path(wgraph_ptr, new_left_pos, { 0, 1, 3 }, 2);
        ASSERT_EQ(expected_path, extended_path);
    }
    {
        GraphPath path(wgraph_ptr, 1, { 0 }, 2);
        GraphPath extended_path = path.extendEndNodeTo(1);
        GraphPath expected_path(wgraph_ptr, 1, { 0, 1 }, 0);
        ASSERT_EQ(expected_path, extended_path);
    }
}

TEST_F(SwapGraph, ExtendingPathToNonadjacentNodeCausesError)
{
    {
        GraphPath path(wgraph_ptr, 1, { 2, 3 }, 3);
        ASSERT_ANY_THROW(path.extendStartNodeTo(1));
    }
    {
        GraphPath path(wgraph_ptr, 1, { 0 }, 3);
        ASSERT_ANY_THROW(path.extendEndNodeTo(3));
    }
}

TEST_F(DeletionGraph, EndOfPathIsExtendedBySpecifiedLength)
{
    GraphPath path(wgraph_ptr, 2, { 0 }, 4);
    {
        const int32_t start_extension = 0;
        const int32_t end_extension = 1;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        const list<GraphPath> expected_path_extensions = { GraphPath(wgraph_ptr, 2, { 0 }, 5) };
        ASSERT_EQ(expected_path_extensions, path_extensions);
    }

    {
        const int32_t start_extension = 0;
        const int32_t end_extension = 2;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        const list<GraphPath> expected_path_extensions
            = { GraphPath(wgraph_ptr, 2, { 0, 1 }, 0), GraphPath(wgraph_ptr, 2, { 0, 2 }, 0) };
        ASSERT_EQ(expected_path_extensions, path_extensions);
    }
}

TEST_F(SwapGraph, EndOfPathIsExtendedBySpecifiedLength)
{
    const GraphPath path(wgraph_ptr, 0, { 0 }, 4);
    const int32_t start_extension = 0;
    const int32_t end_extension = 5;
    const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions
        = { GraphPath(wgraph_ptr, 0, { 0, 1, 3 }, 0), GraphPath(wgraph_ptr, 0, { 0, 2, 3 }, 1) };
    ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST_F(DeletionGraph, PathIsExtendedBySpecifiedLength)
{
    GraphPath path(wgraph_ptr, 2, { 0 }, 4);
    {
        const int32_t start_extension = 2;
        const int32_t end_extension = 1;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        const list<GraphPath> expected_path_extensions = { GraphPath(wgraph_ptr, 0, { 0 }, 5) };
        ASSERT_EQ(expected_path_extensions, path_extensions);
    }
}

TEST_F(DoubleSwapGraph, PathIsExtendedBySpecifiedLength)
{
    GraphPath path(wgraph_ptr, 1, { 3 }, 3);
    {
        const int32_t start_extension = 5;
        const int32_t end_extension = 1;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        const list<GraphPath> expected_path_extensions
            = { GraphPath(wgraph_ptr, 0, { 1, 3, 4 }, 0), GraphPath(wgraph_ptr, 0, { 1, 3, 5 }, 0),
                GraphPath(wgraph_ptr, 4, { 0, 2, 3, 4 }, 0), GraphPath(wgraph_ptr, 4, { 0, 2, 3, 5 }, 0) };
        ASSERT_EQ(expected_path_extensions, path_extensions);
    }
}

TEST_F(DoubleSwapGraph, PathExtensionsThatAreTooLongAreNotOutput)
{

    GraphPath path(wgraph_ptr, 1, { 3 }, 3);
    {
        const int32_t start_extension = 50;
        const int32_t end_extension = 1;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        ASSERT_TRUE(path_extensions.empty());
    }
    {
        const int32_t start_extension = 10;
        const int32_t end_extension = 9;
        const list<GraphPath> path_extensions = path.extendBy(start_extension, end_extension);
        const list<GraphPath> expected_path_extensions = { GraphPath(wgraph_ptr, 0, { 0, 1, 3, 4, 6 }, 3) };
        ASSERT_EQ(expected_path_extensions, path_extensions);
    }
}