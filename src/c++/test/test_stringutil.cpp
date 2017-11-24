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
