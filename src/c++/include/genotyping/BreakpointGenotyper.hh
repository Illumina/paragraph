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

/**
 * structure for passing parameter into breakpoint genotyper
 */
struct BreakpointGenotyperParameter
{
    BreakpointGenotyperParameter(double _read_depth, int32_t _read_length, double _depth_sd, bool _use_poisson_depth)
        : read_depth(_read_depth)
        , read_length(_read_length)
        , depth_sd(_depth_sd)
        , use_poisson_depth(_use_poisson_depth){};

    double read_depth;
    int32_t read_length;
    double depth_sd;
    bool use_poisson_depth;
};

class BreakpointGenotyper
{
public:
    /**
     * Initialize a breakpoint genotyper through genotype parameters
     * @param param Pointer to genotyping parameters
     */
    explicit BreakpointGenotyper(std::unique_ptr<GenotypingParameters> const& param);

    /**
     * public function to do breakpoint genotyping from read counts on edges
     * @param read_depth Expected / mean read depth
     * @param read_length Read length
     * @param read counts Number of reads for each allele
     * @param depth variance Variance of depth across genome. If non-positive, Poission test will be used
     */
    Genotype
    genotype(const BreakpointGenotyperParameter& param, const std::vector<int32_t>& read_counts_per_allele) const;

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
     * cutoff for coverage test p value
     * 1. lower end
     * 2. upper end
     */
    std::pair<double, double> coverage_test_cutoff_;

    /**
     * minimum GQ for a PASS event
     */
    int min_pass_gq_;

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