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

/**
 * \Storage for pooled genotyping stats
 *
 * \author Sai Chen
 * \email schen6@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/VariantGenotype.hh"
#include <map>
#include <vector>

namespace genotyping
{
class PopStats
{
public:
    /**
     * Population stats should be constructed from a vector of genotypes
     * The genotype can be variant genotype (votted result from breakpoints)
     *      Or simple breakpoint genotypes
     */
    PopStats(std::vector<Genotype>& breakpoint_genotypes);
    PopStats(std::vector<VariantGenotype>& variant_genotypes);

    /**
     *  simple getters
     */
    std::map<uint64_t, int> alleleCounts() { return allele_counts; };

    /**
     *  advanced population-scale stats
     */
    double callRate() { return (double)num_valid_samples / num_total_samples; };

    /**
     *  pearson chisq test for HWE
     */
    double getHWE();

private:
    std::map<GenotypeVector, int> genotype_counts; // genotype -> #samples with this genotype
    std::map<uint64_t, int> allele_counts; // allele -> count of this allele in all samples
    int num_total_samples;
    int num_valid_samples;
};
};