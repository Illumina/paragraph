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
#include "genotyping/GenotypingParameters.hh"
#include "gmock/gmock.h"
#include <math.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
using namespace genotyping;

TEST(BreakpointGenotyper, ThrowsWhenWrongNumberOfReadcounts)
{
    double read_depth = 40.0;
    int32_t read_length = 100;
    double depth_sd = sqrt(read_depth * 5);
    const vector<string> alleles = { "REF", "ALT" };

    auto param = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles, 2));
    BreakpointGenotyper genotyper(param);
    const BreakpointGenotyperParameter b_param(read_depth, read_length, depth_sd, false);
    ASSERT_ANY_THROW(genotyper.genotype(b_param, {}));
    ASSERT_ANY_THROW(genotyper.genotype(b_param, { 10 }));
}

TEST(BreakpointGenotyper, GenotypesWellCoveredBreakpoints)
{
    const vector<string> alleles = { "REF", "ALT" };
    auto param = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles, 2));
    BreakpointGenotyper genotyper(param);

    double read_depth = 40.0;
    int32_t read_length = 100;
    double depth_sd = 20;
    const BreakpointGenotyperParameter b_param(read_depth, read_length, depth_sd, false);

    EXPECT_EQ("0/0", (string)genotyper.genotype(b_param, { 20, 0 }));
    EXPECT_EQ("0/1", (string)genotyper.genotype(b_param, { 20, 20 }));
    EXPECT_EQ("1/1", (string)genotyper.genotype(b_param, { 0, 20 }));

    // test chrX
    auto haploid_param = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles, 1));
    BreakpointGenotyper haploid_genotyper(haploid_param);
    EXPECT_EQ("1", (string)haploid_genotyper.genotype(b_param, { 0, 20 }));

    EXPECT_FLOAT_EQ(0.24825223, (float)genotyper.genotype(b_param, { 0, 20 }).coverage_test_pvalue);

    const BreakpointGenotyperParameter b_poisson_param(read_depth, read_length, depth_sd, true);
    EXPECT_FLOAT_EQ(0.0080560343, (float)genotyper.genotype(b_poisson_param, { 0, 20 }).coverage_test_pvalue);

    const vector<string> alleles2 = { "REF", "ALT1", "ALT2", "ALT3", "ALT4" };
    auto param2 = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles2, 2));
    BreakpointGenotyper genotyper_q(param2);
    EXPECT_EQ("1/3", (string)genotyper_q.genotype(b_param, { 1, 20, 2, 20, 2 }));
}