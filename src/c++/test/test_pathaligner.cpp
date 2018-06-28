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

#include "graphcore/GraphBuilders.hh"

#include <list>

#include "gtest/gtest.h"

using common::Read;
using graphtools::Graph;
using graphtools::Path;
using grm::PathAligner;
using std::list;

using namespace testing;
using namespace common;

TEST(PathAligner, Aligns_ExactMatch)
{
    Graph g = graphtools::makeDeletionGraph("AAAAAAAAA", "CCCC", "GGGGGGGGG");

    list<Path> paths{};

    PathAligner aligner(16);
    aligner.setGraph(&g, paths);

    {
        Read read;
        read.setCoreInfo("f1", "AAAAAAAAGGGGGGGG", "################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(1, read.graph_pos());
        ASSERT_EQ("0[8M]2[8M]", read.graph_cigar());
        ASSERT_EQ(16, read.graph_alignment_score());
        ASSERT_FALSE(read.is_graph_reverse_strand());
    }

    {
        Read read;
        read.setCoreInfo("f1", "CCCCCCCCTTTTTTTT", "################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(1, read.graph_pos());
        ASSERT_EQ("0[8M]2[8M]", read.graph_cigar());
        ASSERT_EQ(16, read.graph_alignment_score());
        ASSERT_TRUE(read.is_graph_reverse_strand());
    }

    {
        Read read;
        read.setCoreInfo("f1", "AAAAAAAACCCCGGGG", "################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(1, read.graph_pos());
        ASSERT_EQ("0[8M]1[4M]2[4M]", read.graph_cigar());
        ASSERT_EQ(16, read.graph_alignment_score());
        ASSERT_FALSE(read.is_graph_reverse_strand());
    }

    {
        Read read;
        read.setCoreInfo("f1", "CCCCGGGGTTTTTTTT", "################");
        read.set_is_reverse_strand(true);
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(1, read.graph_pos());
        ASSERT_EQ("0[8M]1[4M]2[4M]", read.graph_cigar());
        ASSERT_EQ(16, read.graph_alignment_score());
        ASSERT_TRUE(read.is_graph_reverse_strand());
    }
}

TEST(PathAligner, Aligns_ExactMatchLongMEM)
{
    Graph g = graphtools::makeDeletionGraph("AAAAAAAAA", "CCCC", "GGGGGGGGG");

    list<Path> paths{};

    PathAligner aligner(16);
    aligner.setGraph(&g, paths);

    {
        Read read;
        read.setCoreInfo("f1", "AAAAAAAAGGGGGGGGG", "################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(1, read.graph_pos());
        ASSERT_EQ("0[8M]2[9M]", read.graph_cigar());
        ASSERT_EQ(17, read.graph_alignment_score());
        ASSERT_FALSE(read.is_graph_reverse_strand());
    }
    {
        Read read;
        read.setCoreInfo("f1", "CCCCCCCCCTTTTTTTTT", "################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(0, read.graph_pos());
        ASSERT_EQ("0[9M]2[9M]", read.graph_cigar());
        ASSERT_EQ(18, read.graph_alignment_score());
        ASSERT_TRUE(read.is_graph_reverse_strand());
    }
}

TEST(PathAligner, Aligns_MultipleMatches)
{
    Graph g = graphtools::makeDeletionGraph("GGGGGGGGGGGG", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGTGGG");

    list<Path> paths{};

    PathAligner aligner(16);
    aligner.setGraph(&g, paths);

    {
        Read read;
        read.setCoreInfo("f1", "CCCCCCCCCCCCGGGGGGGGGGGG", "#####################################");
        aligner.alignRead(read);
        ASSERT_EQ(common::Read::MAPPED, read.graph_mapping_status());
        ASSERT_EQ(4, read.graph_pos());
        ASSERT_EQ("1[12M]2[12M]", read.graph_cigar());
        ASSERT_EQ(24, read.graph_alignment_score());
        ASSERT_FALSE(read.is_graph_reverse_strand());
        ASSERT_FALSE(read.is_graph_alignment_unique());
        ASSERT_EQ(0, read.graph_mapq());
    }
}
