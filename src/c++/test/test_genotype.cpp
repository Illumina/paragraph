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
