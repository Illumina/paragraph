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
    EXPECT_EQ("ALT/ALT", combined_genotype.toString(&alleles));
}

TEST(CombinedGenotype, GenotypeUnphasedMatch)
{
    const vector<string> alleles{ "REF", "ALT" };

    Genotype gt1;
    gt1.gt = { 0, 1 };
    gt1.gl_name = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    gt1.gl = { -10, -0.1, -10 };
    gt1.gq = 20;

    Genotype gt2;
    gt2.gt = { 1, 0 };
    gt2.gl_name = { { 1, 0 }, { 1, 1 }, { 0, 0 } };
    gt2.gl = { -0.1, -10, -10 };
    gt2.gq = 30;

    GenotypeSet gs;
    gs.add(alleles, gt1);
    gs.add(alleles, gt2);

    const Genotype combined_genotype = combinedGenotype(gs);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("PASS", combined_genotype.filterString());
    EXPECT_EQ(20, combined_genotype.gq);
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

    auto param = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles, 2));
    BreakpointGenotyper genotyper(param);

    double read_depth = 10.0;
    int32_t read_length = 100;
    double depth_sd = 50;
    const BreakpointGenotyperParameter b_param(read_depth, read_length, depth_sd, false);
    const Genotype combined_genotype = combinedGenotype(gs, &b_param, &genotyper);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("CONFLICT", combined_genotype.filterString());
    EXPECT_EQ(8, combined_genotype.gq);

    // additional test for chrX conflict genotypes
    auto haploid_param = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleles, 1));
    BreakpointGenotyper haploid_genotyper(haploid_param);
    const BreakpointGenotyperParameter b_param2(read_depth, read_length, depth_sd, false);
    Genotype gtX1;
    gtX1.gt = { 0 };
    gtX1.num_reads = 10;
    gtX1.allele_fractions = { 1, 0 };
    Genotype gtX2;
    gtX2.gt = { 1 };
    gtX2.num_reads = 2;
    gtX2.allele_fractions = { 0, 1 };
    GenotypeSet gsX;
    gsX.add(alleles, gtX1);
    gsX.add(alleles, gtX2);
    const Genotype combined_haploid_genotype = combinedGenotype(gsX, &b_param2, &haploid_genotyper);
    EXPECT_EQ("0", combined_haploid_genotype.toString());
}

TEST(CombinedGenotype, GenotypeMissing)
{
    const vector<string> alleles{ "REF", "ALT" };

    Genotype gt1;

    Genotype gt2;
    gt2.gt = { 0, 1 };
    gt2.gl_name = { { 0, 1 }, { 1, 1 }, { 0, 0 } };
    gt2.gl = { -1, -10, -10 };
    gt2.gq = 36;

    GenotypeSet gs;
    gs.add(alleles, gt1);
    gs.add(alleles, gt2);

    const Genotype combined_genotype = combinedGenotype(gs);

    EXPECT_EQ("0/1", combined_genotype.toString());
    EXPECT_EQ("BP_NO_GT", combined_genotype.filterString());
    EXPECT_EQ(36, combined_genotype.gq);
}
