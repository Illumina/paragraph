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
