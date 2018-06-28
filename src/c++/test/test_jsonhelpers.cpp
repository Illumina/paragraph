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
 * \brief Test JSON I/O
 *
 * \file test_jsonhelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/JsonHelpers.hh"
#include "gtest/gtest.h"

#include <limits>
#include <string>

using std::string;

TEST(JsonHelpers, ReadJson)
{
    Json::Value expected;
    expected["a"] = 0;
    expected["b"] = 1.2;
    expected["c"] = "abc";
    expected["d"] = Json::arrayValue;
    expected["d"].append(1);
    expected["d"].append(2);

    const string json_text = "{ \"a\": 0,"
                             "  \"b\": 1.2,"
                             "  \"c\": \"abc\","
                             "  \"d\": [1, 2] "
                             "}";
    const Json::Value observed = common::getJSON(json_text);
    ASSERT_EQ(expected, observed);
}

TEST(JsonHelpers, WritesInfAndNaN)
{
    Json::Value input;
    input["a"] = std::numeric_limits<double>::infinity();
    input["b"] = -std::numeric_limits<double>::infinity();
    input["c"] = -std::numeric_limits<double>::quiet_NaN();

    const string expected = "{\n\t\"a\" : 1e+9999,\n\t\"b\" : -1e+9999,\n\t\"c\" : null\n}";
    const string observed = common::writeJson(input);
    ASSERT_EQ(expected, observed);
}
