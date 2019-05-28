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

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"

#include "paragraph/ReadFilter.hh"

using std::string;
using std::vector;

using namespace graphtools;

TEST(ReadFilter, FilterNonUniq)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");

    const string query1 = "AAAATTCCC";
    common::Read read1("read1", query1, string(query1.size(), '#'));
    read1.set_graph_cigar("0[4M]1[2M3S]");
    read1.set_graph_alignment_score(6);
    read1.set_graph_mapping_status(common::Read::MAPPED);
    read1.set_is_graph_alignment_unique(false);

    const string query2 = "AAAATTGG";
    common::Read read2("read2", query2, string(query2.size(), '#'));
    read2.set_graph_cigar("0[4M]1[4M]");
    read2.set_graph_alignment_score(8);
    read2.set_graph_mapping_status(common::Read::MAPPED);
    read2.set_is_graph_alignment_unique(true);

    auto read_filter = paragraph::createReadFilter(&graph, true, 0.0, 0);

    const auto read1_result = read_filter->filterRead(read1);
    const auto read2_result = read_filter->filterRead(read2);

    ASSERT_TRUE(read1_result.first);
    ASSERT_EQ("nonuniq", read1_result.second);
    ASSERT_FALSE(read2_result.first);
    ASSERT_EQ("", read2_result.second);
}

TEST(ReadFilter, FilterBadAlign)
{
    Graph graph = makeDeletionGraph("AAAA", "GGGG", "TTTT");

    const string query1 = "AAAACCCCCCCC";
    common::Read read1("read1", query1, string(query1.size(), '#'));
    read1.set_graph_cigar("0[4M8S]");
    read1.set_graph_alignment_score(4);
    read1.set_graph_mapping_status(common::Read::MAPPED);
    read1.set_is_graph_alignment_unique(true);

    const string query2 = "AAAAGCCCCCCC";
    common::Read read2("read2", query2, string(query2.size(), '#'));
    read2.set_graph_cigar("0[4M]1[1M7S]");
    read2.set_graph_alignment_score(5);
    read2.set_graph_mapping_status(common::Read::MAPPED);
    read2.set_is_graph_alignment_unique(true);

    auto read_filter = paragraph::createReadFilter(&graph, true, 0.4, 0);

    const auto read1_result = read_filter->filterRead(read1);
    const auto read2_result = read_filter->filterRead(read2);

    ASSERT_TRUE(read1_result.first);
    ASSERT_EQ("bad_align", read1_result.second);
    ASSERT_FALSE(read2_result.first);
    ASSERT_EQ("", read2_result.second);
}

TEST(ReadFilter, FilterKmers)
{
    Graph graph = makeDeletionGraph("AGAG", "TTGG", "TTT");
    auto read_filter = paragraph::createReadFilter(&graph, false, 0.0, 3);
    {
        const string query = "AGAGTT";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]1[2M]");
        read.set_graph_alignment_score(6);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_TRUE(read_result.first);
        ASSERT_EQ("kmer_uncov_1", read_result.second);
    }
    {
        const string query = "AGAGTTT";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]2[3M]");
        read.set_graph_alignment_score(7);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_FALSE(read_result.first);
        ASSERT_EQ("", read_result.second);
    }
}

TEST(ReadFilter, FilterKmersSnpMismatch)
{
    Graph graph = makeSwapGraph("AGAG", "T", "C", "ACAC");
    auto read_filter = paragraph::createReadFilter(&graph, false, 0.0, 4);
    {
        const string query = "AGAGGACAC";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]1[1X]3[4M]");
        read.set_graph_alignment_score(8);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_EQ(true, read_result.first);
        ASSERT_EQ("kmer_uncov_1", read_result.second);
    }
    {
        const string query = "AGAGTACAC";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]1[1M]3[4M]");
        read.set_graph_alignment_score(8);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_FALSE(read_result.first);
        ASSERT_EQ("", read_result.second);
    }
    {
        const string query = "AGAGTACAC";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]2[1X]3[4M]");
        read.set_graph_alignment_score(8);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_TRUE(read_result.first);
        ASSERT_EQ("kmer_uncov_2", read_result.second);
    }
    {
        const string query = "AGAGCACAC";
        common::Read read("read", query, string(query.size(), '#'));
        read.set_graph_cigar("0[4M]2[1M]3[4M]");
        read.set_graph_alignment_score(8);
        read.set_graph_mapping_status(common::Read::MAPPED);

        const auto read_result = read_filter->filterRead(read);
        ASSERT_FALSE(read_result.first);
        ASSERT_EQ("", read_result.second);
    }
}
