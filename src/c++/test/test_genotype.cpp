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

#include "genotyping/Genotype.hh"
#include "gtest/gtest.h"
#include <vector>

using std::vector;
using namespace genotyping;

TEST(Genotype, Relabel)
{
    std::vector<uint64_t> new_labels{ 1, 3 };

    Genotype variant;
    variant.gt = { 0, 1 };
    variant.gl_name.emplace_back(std::initializer_list<uint64_t>{ 0, 0 });
    variant.gl_name.emplace_back(std::initializer_list<uint64_t>{ 0, 1 });
    variant.gl_name.emplace_back(std::initializer_list<uint64_t>{ 1, 1 });

    variant.relabel(new_labels);

    EXPECT_EQ("1/3", variant.toString());
    EXPECT_EQ(variant.gl_name[0][0], 1ull);
    EXPECT_EQ(variant.gl_name[0][1], 1ull);
    EXPECT_EQ(variant.gl_name[1][0], 1ull);
    EXPECT_EQ(variant.gl_name[1][1], 3ull);
    EXPECT_EQ(variant.gl_name[2][0], 3ull);
    EXPECT_EQ(variant.gl_name[2][1], 3ull);
}
