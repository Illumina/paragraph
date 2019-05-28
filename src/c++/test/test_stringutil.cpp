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
 *
 * \file test_stringutil.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "gtest/gtest.h"

#include "common/StringUtil.hh"

using namespace common::stringutil;

TEST(String, StringSplit)
{
    std::vector<std::string> v;

    split("a, b, c", v);

    ASSERT_TRUE(v.size() == 3);
    ASSERT_TRUE(v[0] == "a");
    ASSERT_TRUE(v[1] == "b");
    ASSERT_TRUE(v[2] == "c");

    v.clear();

    split("a, b, c", v, " ,", true);

    ASSERT_TRUE(v.size() == 5);
    ASSERT_TRUE(v[0] == "a");
    ASSERT_TRUE(v[1] == "");
    ASSERT_TRUE(v[2] == "b");
    ASSERT_TRUE(v[3] == "");
    ASSERT_TRUE(v[4] == "c");

    v.clear();

    split("a: b- c", v, " :-");

    ASSERT_TRUE(v.size() == 3);
    ASSERT_TRUE(v[0] == "a");
    ASSERT_TRUE(v[1] == "b");
    ASSERT_TRUE(v[2] == "c");
}

TEST(String, StringEndsWith)
{
    ASSERT_TRUE(endsWith("abc", "bc"));
    ASSERT_TRUE(!endsWith("abc", "cb"));
}

TEST(String, StringReplaceAll)
{
    ASSERT_EQ(replaceAll("1,000,000", ",", ""), "1000000");
    ASSERT_EQ(replaceAll("1,000,000", ",0", "2"), "1200200");
}

TEST(String, StringFormatPos)
{
    ASSERT_EQ(formatPos("chr1"), "chr1");
    ASSERT_EQ(formatPos("chr1", 999), "chr1:1000");
    ASSERT_EQ(formatPos("chr1", 999, 1999), "chr1:1000-2000");
}

TEST(String, StringParsePos)
{
    std::string chr = "-";
    int64_t start = -1, end = -1;

    parsePos("chr1", chr, start, end);

    ASSERT_EQ(chr, "chr1");
    ASSERT_EQ(start, -1);
    ASSERT_EQ(end, -1);

    chr = "-";
    start = -1;
    end = -1;

    parsePos("chr1:1,000", chr, start, end);

    ASSERT_EQ(chr, "chr1");
    ASSERT_EQ(start, 999);
    ASSERT_EQ(end, -1);

    chr = "-";
    start = -1;
    end = -1;

    parsePos("chr1:1,000-2000", chr, start, end);

    ASSERT_EQ(chr, "chr1");
    ASSERT_EQ(start, 999);
    ASSERT_EQ(end, 1999);
}
