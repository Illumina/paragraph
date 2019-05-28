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

#include "gmock/gmock.h"

#include "common/Read.hh"
#include "common/ReadPair.hh"

using ::testing::Test;

using common::Read;
using common::ReadPair;

class AReadPair : public Test
{
public:
    ReadPair mate_pair;
    Read read;

    virtual void SetUp()
    {
        read.set_fragment_id("Fragment1");
        read.set_bases("ATCG");
        read.set_quals("####");
    }
};

TEST_F(AReadPair, AddsFirstMate)
{
    read.set_is_first_mate(true);
    mate_pair.add(read);
    ASSERT_EQ(read, mate_pair.first_mate());
}

TEST_F(AReadPair, AddsSecondMate)
{
    read.set_is_first_mate(false);
    mate_pair.add(read);
    ASSERT_EQ(read, mate_pair.second_mate());
}
