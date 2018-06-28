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
#include "common/Error.hh"
#include <algorithm>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <limits>
#include <numeric>

using boost::math::poisson_distribution;
using std::map;
using std::vector;

namespace genotyping
{

/**
 * @param genotype_parameters GenotypeParameter class
 */
BreakpointGenotyper::BreakpointGenotyper(std::unique_ptr<GenotypingParameters> const& param)
    : n_alleles_(param->numAlleles())
    , ploidy_(param->ploidy())
    , coverage_test_cutoff_(param->coverageTestCutoff())
    , min_overlap_bases_(param->minOverlapBases())
    , possible_genotypes(param->possibleGenotypes())
{
    if (param->alleleErrorRates().empty())
    {
        allele_error_rate_.push_back(param->otherAlleleErrorRate());
    }
    else
    {
        allele_error_rate_ = param->alleleErrorRates();
    }

    if (param->hetHaplotypeFractions().empty())
    {
        haplotype_read_fraction_.push_back(param->otherHetHaplotypeFraction());
    }
    else
    {
        haplotype_read_fraction_ = param->hetHaplotypeFractions();
    }

    if (!param->genotypeFractions().empty())
    {
        genotype_prior_ = param->genotypeFractions();
        for (auto& phi : genotype_prior_)
        {
            if (phi.first.size() < ploidy_)
            {
                error("Error: genotype and ploidy does not match.");
            }
            if (phi.second < 0 || phi.second > 1)
            {
                error("Error: genotype prior should be between 0~1.");
            }
            phi.second = log(phi.second);
        }
    }
};

Genotype BreakpointGenotyper::genotype(double read_depth, int32_t read_length, const vector<int32_t>& read_counts) const
{
    if (read_counts.size() != n_alleles_)
    {
        error("Error: number of read counts and alleles mismatches. %i != %i.", (int)read_counts.size(), n_alleles_);
    }
    Genotype result;

    // compute adjusted depth
    const double multiplier = (read_length - min_overlap_bases_) / (double)read_length;
    assert(multiplier > 0);
    const double lambda = read_depth * multiplier;
    const int32_t total_num_reads = std::accumulate(read_counts.begin(), read_counts.end(), 0);
    if (total_num_reads == 0)
    {
        result.filter = "NO_READS";
        return result;
    }
    result.num_reads = total_num_reads;

    // compute GL and GT
    double best_gl = -std::numeric_limits<double>::max();

    for (const auto& igt : possible_genotypes)
    {
        const double gl = genotypeLikelihood(lambda, igt, read_counts);
        result.gl_name.push_back(igt);
        result.gl.push_back(gl);

        // update GT if GL is better
        if (gl > best_gl)
        {
            best_gl = gl;
            result.gt = igt;
        }
    }

    // compute allele fractions
    result.allele_fractions.resize(n_alleles_, 0.0);
    for (unsigned int al = 0; al < n_alleles_; ++al)
    {
        result.allele_fractions[al] = ((double)read_counts[al]) / total_num_reads;
    }

    // compute coverage test p value
    const poisson_distribution<> coverage_distribution(lambda);
    double coverage_test_pvalue = cdf(coverage_distribution, total_num_reads);
    if (coverage_test_pvalue > 0.5)
    {
        coverage_test_pvalue = 1 - coverage_test_pvalue;
    }
    result.coverage_test_pvalue = coverage_test_pvalue;

    if (coverage_test_pvalue < coverage_test_cutoff_)
    {
        result.filter = "DEPTH";
    }

    return result;
}

/**
 * return genotype likelihood given genotype (using Poisson model with internal parameters)
 * @param lambda Poisson distribution parameter
 * @param gv Genotype vector
 * @param read_counts Read count vector for each allele
 */
double BreakpointGenotyper::genotypeLikelihood(
    double lambda, const GenotypeVector& gv, const vector<int32_t>& read_counts) const
{
    double log_phi;
    auto it = genotype_prior_.find(gv);
    if (it == genotype_prior_.end())
    {
        log_phi = 0;
    }
    else
    {
        log_phi = it->second;
    }

    // compute how many copies of each allele we have in this GT
    vector<int> allele_ploidy;
    allele_ploidy.resize(n_alleles_, 0);
    for (unsigned int al = 0; al < n_alleles_; ++al)
    {
        for (const auto g : gv)
        {
            if (al == g)
            {
                ++allele_ploidy[al];
            }
        }
    }

    // compute GL by summing all allele contributions
    double gl = log_phi;
    for (unsigned int al = 0; al < n_alleles_; ++al)
    {
        if (allele_ploidy[al] == 0)
        {
            // no copies -> all reads supporting this allele will be errors
            const double eps = (allele_error_rate_.size() == 1 ? allele_error_rate_[0] : allele_error_rate_[al]);
            const poisson_distribution<> error_distribution(lambda * eps);
            gl += log(pdf(error_distribution, read_counts[al]));
        }
        else
        {
            const double mu
                = (haplotype_read_fraction_.size() == 1 ? haplotype_read_fraction_[0] : haplotype_read_fraction_[al]);
            const double mu_with_ploidy = allele_ploidy[al] * mu;

            const poisson_distribution<> allele_count_distribution(lambda * mu_with_ploidy);
            gl += log(pdf(allele_count_distribution, read_counts[al]));
        }
        if (std::isinf(gl))
        {
            return -std::numeric_limits<double>::max();
        }
    }

    return gl;
}
}
