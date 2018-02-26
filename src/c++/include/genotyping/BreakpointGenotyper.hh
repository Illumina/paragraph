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
 * Simple Poisson model for genotyping a single diploid breakpoint from a single sample
 *
 *  A full description of the model can be found in doc/graph-models.md
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/GenotypingParameters.hh"

#include <map>
#include <vector>

namespace genotyping
{
class BreakpointGenotyper
{
public:
    /**
     * Initialize a breakpoint genotyper through genotype parameters
     */
    explicit BreakpointGenotyper(GenotypingParameters* param);

    /**
     * public function to do breakpoint genotyping from read counts on edges
     * @param read_depth expected / mean read depth
     * @param read_length read length
     * @param read counts for each allele
     */
    Genotype genotype(double read_depth, int32_t read_length, const std::vector<int32_t>& read_counts_per_allele) const;

private:
    /**
     * return genotype likelihood given genotype (using Poisson model with internal parameters)
     * @param lambda Poisson distribution parameter
     * @param gv Genotype vector
     * @param read_counts Read count vector for each allele
     */
    double genotypeLikelihood(double lambda, const GenotypeVector& gv, const std::vector<int32_t>& read_counts) const;

    unsigned int n_alleles_;
    unsigned int ploidy_;

    /**
     *  genotyping specific stats
     */
    int32_t min_overlap_bases_; // Minimum number of bases that a high-confidence alignment must overlap a node.

    std::vector<GenotypeVector> possible_genotypes; // all possible genotypes

    std::vector<double> allele_error_rate_; // error rate for each allele

    std::vector<double> haplotype_read_fraction_; // mean expected haplotype fraction for each allele

    std::map<GenotypeVector, double> genotype_prior_; // genotype prior log probabilities
};
};