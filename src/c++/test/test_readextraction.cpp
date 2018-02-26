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

#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "common/Read.hh"
#include "common/ReadExtraction.hh"
#include "common/ReadPairs.hh"
#include "common/ReadReader.hh"
#include "common/Region.hh"

using std::vector;
using namespace testing;
using namespace common;

class MockReader : public ReadReader
{
public:
    MOCK_METHOD1(getAlign, bool(Read&));
    MOCK_METHOD2(getAlignedMate, bool(const Read&, Read&));
    MOCK_METHOD1(setRegion, void(const std::string&));
};

class ExtractReads : public Test
{
public:
    MockReader reader;
    ReadPairs read_pairs;
    Read read1;
    Read read2;
    virtual void SetUp()
    {
        read1.set_chrom_id(1);
        read1.set_pos(100);
        read1.set_fragment_id("Fragment_1");
        read1.set_bases("AAAA");
        read1.set_quals("####");
        read2.set_chrom_id(1);
        read2.set_pos(100);
        read2.set_fragment_id("Fragment_2");
        read2.set_bases("AAAA");
        read2.set_quals("####");
    }
};

class RecoverMissingMates : public Test
{
public:
    ReadPairs read_pairs;
    Read read_with_anomalous_mate_a;
    Read mate_of_read_a;
    Read read_with_normal_mate_b;
    Read read_with_anomalous_mate_c;
    Read mate_of_read_c;

    virtual void SetUp()
    {
        read_with_anomalous_mate_a.setCoreInfo("Fragment_1", "AAAA", "####");
        read_with_anomalous_mate_a.set_is_first_mate(true);
        read_with_anomalous_mate_a.set_chrom_id(1);
        read_with_anomalous_mate_a.set_pos(100);
        read_with_anomalous_mate_a.set_mate_chrom_id(1);
        read_with_anomalous_mate_a.set_mate_pos(1600);

        mate_of_read_a.setCoreInfo("Fragment_1", "TTTT", "####");
        mate_of_read_a.set_is_first_mate(false);

        read_with_normal_mate_b.setCoreInfo("Fragment_2", "CCCC", "####");
        read_with_normal_mate_b.set_is_first_mate(true);
        read_with_normal_mate_b.set_chrom_id(3);
        read_with_normal_mate_b.set_pos(500);
        read_with_normal_mate_b.set_mate_chrom_id(3);
        read_with_normal_mate_b.set_mate_pos(800);

        read_with_anomalous_mate_c.setCoreInfo("Fragment_3", "AAAA", "####");
        read_with_anomalous_mate_c.set_is_first_mate(false);
        read_with_anomalous_mate_c.set_chrom_id(5);
        read_with_anomalous_mate_c.set_pos(500);
        read_with_anomalous_mate_c.set_mate_chrom_id(3);
        read_with_anomalous_mate_c.set_mate_pos(500);

        mate_of_read_c.setCoreInfo("Fragment_3", "GGGG", "####");
        mate_of_read_c.set_is_first_mate(true);
    }
};

TEST_F(ExtractReads, ExtractsAllReadsFromReader)
{
    EXPECT_CALL(reader, getAlign(_))
        .WillOnce(DoAll(SetArgReferee<0>(read1), Return(true)))
        .WillOnce(DoAll(SetArgReferee<0>(read2), Return(true)))
        .WillOnce(Return(false));

    int max_reads = 10;
    const std::string chrom = std::string("1");
    const Region region(chrom, 0, 1800);
    vector<Read> observed_reads;
    extractMappedReadsFromRegion(read_pairs, max_reads, reader, region);
    read_pairs.getReads(observed_reads);

    vector<Read> expected_reads = { read1, read2 };
    ASSERT_EQ(expected_reads, observed_reads);
}

TEST_F(ExtractReads, ExtractsMaxAllowedReadsFromReader)
{
    const int max_reads = 1;
    EXPECT_CALL(reader, getAlign(_)).WillOnce(DoAll(SetArgReferee<0>(read1), Return(true)));

    const std::string chrom = std::string("1");
    const Region region(chrom, 0, 1800);
    vector<Read> observed_reads;
    extractMappedReadsFromRegion(read_pairs, max_reads, reader, region);
    read_pairs.getReads(observed_reads);

    vector<Read> expected_reads = { read1 };
    ASSERT_EQ(expected_reads, observed_reads);
}

TEST_F(RecoverMissingMates, RecoversAnomalousMates)
{
    read_pairs.add(read_with_anomalous_mate_a);
    read_pairs.add(read_with_normal_mate_b);
    read_pairs.add(read_with_anomalous_mate_c);

    MockReader reader;
    EXPECT_CALL(reader, getAlignedMate(read_with_anomalous_mate_a, _))
        .WillOnce(DoAll(SetArgReferee<1>(mate_of_read_a), Return(true)));
    EXPECT_CALL(reader, getAlignedMate(read_with_anomalous_mate_c, _))
        .WillOnce(DoAll(SetArgReferee<1>(mate_of_read_c), Return(true)));

    vector<Read> expected_reads = { read_with_anomalous_mate_a, mate_of_read_a, read_with_normal_mate_b, mate_of_read_c,
                                    read_with_anomalous_mate_c };

    recoverMissingMates(reader, read_pairs);
    vector<Read> observed_reads;
    read_pairs.getReads(observed_reads);

    ASSERT_EQ(expected_reads, observed_reads);
}

TEST_F(ExtractReads, isReadOrItsMateInRegion)
{
    const std::string chrom = std::string("1");
    const Region region_left(chrom, 0, 50);
    const Region region_overlap(chrom, 101, 103);
    const Region region_overlap_mate(chrom, 1550, 1650);
    const Region region_right(chrom, 110, 200);
    ASSERT_FALSE(isReadOrItsMateInRegion(read1, region_left));
    ASSERT_TRUE(isReadOrItsMateInRegion(read1, region_overlap));
    ASSERT_FALSE(isReadOrItsMateInRegion(read1, region_right));

    read1.set_mate_chrom_id(1);
    read1.set_mate_pos(1600);
    ASSERT_TRUE(isReadOrItsMateInRegion(read1, region_overlap_mate));
}