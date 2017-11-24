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

#include "genotyping/PopStats.hh"
#include "gtest/gtest.h"

using namespace genotyping;
using std::vector;

TEST(PopStats, SimpleHWE)
{
    vector<VariantGenotype> pop_genotypes;
    VariantGenotype ref({ 0, 0 });
    VariantGenotype het({ 0, 1 });
    VariantGenotype alt({ 1, 1 });

    for (size_t i = 0; i < 10; i++)
    {
        pop_genotypes.push_back(ref);
    }

    // with only one genotype
    PopStats ps0(pop_genotypes);
    EXPECT_EQ(1, ps0.getHWE());

    // with one missing genotype
    for (size_t i = 0; i < 10; i++)
    {
        pop_genotypes.push_back(het);
    }
    PopStats ps1(pop_genotypes);
    EXPECT_LT(0.20, ps1.getHWE());
    EXPECT_GT(0.21, ps1.getHWE());

    // regular case
    pop_genotypes.push_back(alt);
    PopStats ps2(pop_genotypes);
    EXPECT_LT(0.40, ps2.getHWE());
    EXPECT_GT(0.41, ps2.getHWE());
}
