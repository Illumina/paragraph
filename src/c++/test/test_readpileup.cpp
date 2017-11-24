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
