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

#include <vector>

#include "gmock/gmock.h"

#include "common/ReadPairs.hh"

using common::Read;
using common::ReadPair;
using common::ReadPairs;
using std::vector;

class AReadPairContainer : public ::testing::Test
{
public:
    ReadPairs read_pairs;
    Read read1_from_fragment_1;
    Read read2_from_fragment_1;
    Read read2_from_fragment_2;
    virtual void SetUp()
    {
        read1_from_fragment_1.setCoreInfo("frag_1", "ATCG", "####");
        read1_from_fragment_1.set_is_first_mate(true);

        read2_from_fragment_1.setCoreInfo("frag_1", "ATCG", "####");
        read2_from_fragment_1.set_is_first_mate(false);

        read2_from_fragment_2.setCoreInfo("frag_2", "ATCG", "####");
        read2_from_fragment_2.set_is_first_mate(false);
    }
};

TEST_F(AReadPairContainer, ThrowsExceptionWhenAccessingPairThatDoesntExist)
{
    read_pairs.add(read1_from_fragment_1);
    ASSERT_ANY_THROW(read_pairs["missing_fragment"]);
}

TEST_F(AReadPairContainer, AddsReadsFromDistinctFragments)
{
    read_pairs.add(read1_from_fragment_1);
    read_pairs.add(read2_from_fragment_2);
    EXPECT_EQ(read1_from_fragment_1, read_pairs["frag_1"].first_mate());
    EXPECT_EQ(read2_from_fragment_2, read_pairs["frag_2"].second_mate());
}

TEST_F(AReadPairContainer, IsEmptyInitially) { EXPECT_EQ(0, read_pairs.num_reads()); }

TEST_F(AReadPairContainer, UpdatesCountWhenPlacingReadsToNewSpots)
{
    read_pairs.add(read1_from_fragment_1);
    read_pairs.add(read2_from_fragment_2);
    EXPECT_EQ(2, read_pairs.num_reads());
}

TEST_F(AReadPairContainer, KeepsCountWhenPlacingReadsToOccupiedSpots)
{
    read_pairs.add(read1_from_fragment_1);
    read_pairs.add(read1_from_fragment_1);
    EXPECT_EQ(1, read_pairs.num_reads());
}

TEST_F(AReadPairContainer, ReturnsVectorOfReadsItStores)
{
    vector<Read> expected_reads = { read1_from_fragment_1, read2_from_fragment_1, read2_from_fragment_2 };
    for (const auto& read : expected_reads)
    {
        read_pairs.add(read);
    }

    vector<common::Read> observed_reads;
    read_pairs.getReads(observed_reads);

    EXPECT_EQ(expected_reads, observed_reads);
}

TEST_F(AReadPairContainer, ClearesItsContent)
{
    read_pairs.add(read1_from_fragment_1);
    read_pairs.clear();
    EXPECT_EQ(0, read_pairs.num_reads());
}