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
