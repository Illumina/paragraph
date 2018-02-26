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

#include "genotyping/PopulationStatistics.hh"
#include "gtest/gtest.h"

using namespace genotyping;
using std::string;
using std::vector;

TEST(PopStats, SimpleHWE)
{
    GenotypeSet pop_genotypes;
    Genotype ref({ 0, 0 });
    Genotype het({ 0, 1 });
    Genotype alt({ 1, 1 });
    vector<string> alleles = { "REF", "ALT" };

    for (size_t i = 0; i < 83; i++)
    {
        pop_genotypes.add(alleles, ref);
    }

    PopulationStatistics ps0(pop_genotypes);
    EXPECT_EQ(1.0, ps0.getChisqPvalue());

    for (size_t i = 0; i < 13; i++)
    {
        pop_genotypes.add(alleles, het);
    }
    for (size_t i = 0; i < 4; i++)
    {
        pop_genotypes.add(alleles, alt);
    }
    PopulationStatistics ps1(pop_genotypes);
    EXPECT_DOUBLE_EQ(0.0020474148859159769, ps1.getChisqPvalue());
    EXPECT_DOUBLE_EQ(0.010293433548874801, ps1.getFisherExactPvalue());

    // test for multi alleles
    Genotype extra02({ 0, 2 });
    Genotype extra12({ 1, 2 });
    Genotype extra22({ 2, 2 });
    vector<string> multi_alleles = { "REF", "ALT1", "ALT2" };
    GenotypeSet pop_genotypes2;
    for (size_t i = 0; i < 24; i++)
    {
        pop_genotypes2.add(multi_alleles, ref);
    }
    for (size_t i = 0; i < 31; i++)
    {
        pop_genotypes2.add(multi_alleles, het);
    }
    for (size_t i = 0; i < 10; i++)
    {
        pop_genotypes2.add(multi_alleles, alt);
    }
    for (size_t i = 0; i < 19; i++)
    {
        pop_genotypes2.add(multi_alleles, extra02);
    }
    for (size_t i = 0; i < 11; i++)
    {
        pop_genotypes2.add(multi_alleles, extra12);
    }
    for (size_t i = 0; i < 5; i++)
    {
        pop_genotypes2.add(multi_alleles, extra22);
    }
    PopulationStatistics ps2(pop_genotypes2);
    EXPECT_DOUBLE_EQ(0.50000945615245529, ps2.getChisqPvalue());
}
