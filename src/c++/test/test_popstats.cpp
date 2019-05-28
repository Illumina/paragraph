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
