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

/**
 *  \brief Phred conversion
 *
 * \file Phred.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "gtest/gtest.h"

#include "common/ReadPileup.hh"

#include "common/Error.hh"

using common::ReadPileup;
using common::ReadPileupInfo;
using common::SimpleRead;

TEST(ReadPileup, CountsReads)
{
    ReadPileup pileup;

    pileup.addRead(SimpleRead::make(10, 20));
    pileup.addRead(SimpleRead::make(20, 20));
    pileup.addRead(SimpleRead::make(30, 20));
    pileup.addRead(SimpleRead::make(40, 21));
    pileup.addRead(SimpleRead::make(50, 20));
    pileup.addRead(SimpleRead::make(60, 20));

    ASSERT_THROW(pileup.addRead(SimpleRead::make(59)), std::runtime_error);

    const std::set<int64_t> expected_0{ 10, 20 };
    size_t count = 0;
    pileup.pileup(20, [&expected_0, &count](ReadPileupInfo* info) {
        ASSERT_EQ(1ull, expected_0.count(info->pos()));
        ++count;
    });
    ASSERT_EQ(expected_0.size(), count);

    const std::set<int64_t> expected_1{ 30, 40 };
    count = 0;
    pileup.pileup(45, [&expected_1, &count](ReadPileupInfo* info) {
        ASSERT_EQ(1ull, expected_1.count(info->pos()));
        ++count;
    });
    ASSERT_EQ(expected_1.size(), count);

    const std::set<int64_t> expected_2{ 40, 50, 60 };
    count = 0;
    pileup.pileup(60, [&expected_2, &count](ReadPileupInfo* info) {
        ASSERT_EQ(1ull, expected_2.count(info->pos()));
        ++count;
    });
    ASSERT_EQ(expected_2.size(), count);
}

TEST(ReadPileup, FlushesReads)
{
    ReadPileup pileup;

    pileup.addRead(SimpleRead::make(10, 20));
    pileup.addRead(SimpleRead::make(20, 20));
    pileup.addRead(SimpleRead::make(30, 20));
    pileup.addRead(SimpleRead::make(40, 21));
    pileup.addRead(SimpleRead::make(50, 20));
    pileup.addRead(SimpleRead::make(60, 20));

    // this should remove the first two reads
    pileup.flush(41);

    size_t count = 0;
    pileup.pileup(20, [&count](ReadPileupInfo* info) { ++count; });
    ASSERT_EQ((size_t)0, count);

    const std::set<int64_t> expected_1{ 30, 40 };
    count = 0;
    pileup.pileup(45, [&expected_1, &count](ReadPileupInfo* info) {
        ASSERT_EQ(1ull, expected_1.count(info->pos()));
        ++count;
    });
    ASSERT_EQ(expected_1.size(), count);
}
