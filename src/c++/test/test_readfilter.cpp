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
