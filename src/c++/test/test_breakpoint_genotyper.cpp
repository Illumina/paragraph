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

#include "genotyping/BreakpointGenotyper.hh"
#include "gmock/gmock.h"
#include <string>
#include <vector>

using std::string;
using std::vector;
using namespace genotyping;

TEST(BeakpointGenotyper, ErrorsOutIfFewerThanTwoCountsAreGiven)
{
    const double read_depth = 40.0;
    const int32_t read_length = 100;
    BreakpointGenotyper genotyper(read_depth, read_length);
    ASSERT_ANY_THROW(genotyper.genotype({}));
    ASSERT_ANY_THROW(genotyper.genotype({ 10 }));
}

TEST(BeakpointGenotyper, GenotypesWellCoveredBreakpoints)
{
    const double read_depth = 40.0;
    const int32_t read_length = 100;
    BreakpointGenotyper genotyper(read_depth, read_length);

    EXPECT_EQ("0/0", (string)genotyper.genotype({ 20, 0 }));
    EXPECT_EQ("0/1", (string)genotyper.genotype({ 20, 20 }));
    EXPECT_EQ("1/1", (string)genotyper.genotype({ 0, 20 }));
    EXPECT_EQ("1/1", (string)genotyper.genotype({ 0, 20 }, { 0.4, 0.6 }, { 0.1, 0.001, 0.1 }));
    EXPECT_EQ("1/3", (string)genotyper.genotype({ 1, 20, 2, 20, 2 }));
}