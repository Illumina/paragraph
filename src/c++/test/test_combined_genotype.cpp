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
#include "genotyping/CombinedGenotype.hh"
#include "genotyping/Genotype.hh"
#include "genotyping/GenotypingParameters.hh"
#include "gtest/gtest.h"

#include <string>
#include <vector>

using namespace genotyping;

using std::string;
using std::vector;

TEST(CombinedGenotype, SimplePass)
{
    Genotype gt1;

    const vector<string> alleles{ "REF", "ALT" };

    gt1.gt = { 1, 1 };
    gt1.gl_name = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    gt1.gl = { -10, -10, -0.1 };

    GenotypeSet gs;
    for (size_t i = 0; i < 2; i++)
    {
        gs.add(alleles, gt1);
    }
    const Genotype combined_genotype = combinedGenotype(gs);

    EXPECT_EQ("1/1", combined_genotype.toString());
    EXPECT_EQ("ALT/ALT", combined_genotype.toString(alleles));
}

TEST(CombinedGenotype, GenotypeUnphasedMatch)
{
    const vector<string> alleles{ "REF", "ALT" };

    Genotype gt1;
    gt1.gt = { 0, 1 };
    gt1.gl_name = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    gt1.gl = { -10, -0.1, -10 };

    Genotype gt2;
    gt2.gt = { 1, 0 };
    gt2.gl_name = { { 1, 0 }, { 1, 1 }, { 0, 0 } };
    gt2.gl = { -0.1, -10, -10 };

    GenotypeSet gs;
    gs.add(alleles, gt1);
    gs.add(alleles, gt2);

    const Genotype combined_genotype = combinedGenotype(gs);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("PASS", combined_genotype.filter);
}

TEST(CombinedGenotype, GenotypeConflictNoConsensus)
{
    const vector<string> alleles{ "REF", "ALT" };

    Genotype gt1;
    gt1.gt = { 0, 1 };
    gt1.num_reads = 10;
    gt1.allele_fractions = { 0.5, 0.5 };

    Genotype gt2;
    gt2.gt = { 1, 1 };
    gt2.num_reads = 10;
    gt2.allele_fractions = { 0, 1 };

    GenotypeSet gs;
    gs.add(alleles, gt1);
    gs.add(alleles, gt2);

    GenotypingParameters param(alleles);
    BreakpointGenotyper genotyper(&param);

    const double read_depth = 10.0;
    const int32_t read_length = 100;
    const Genotype combined_genotype = combinedGenotype(gs, &genotyper, read_depth, read_length);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("CONFLICT", combined_genotype.filter);
}

TEST(CombinedGenotype, GenotypeMissing)
{
    const vector<string> alleles{ "REF", "ALT" };

    Genotype gt1;

    Genotype gt2;
    gt2.gt = { 0, 1 };
    gt2.gl_name = { { 0, 1 }, { 1, 1 }, { 0, 0 } };
    gt2.gl = { 0.1, -10, -10 };

    GenotypeSet gs;
    gs.add(alleles, gt1);
    gs.add(alleles, gt2);

    const Genotype combined_genotype = combinedGenotype(gs);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("MISSING", combined_genotype.filter);
}
