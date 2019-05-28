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