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