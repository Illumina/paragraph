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
#include "common/Error.hh"
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <limits>
#include <math.h>
#include <numeric>

using boost::math::normal_distribution;
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
    , min_pass_gq_(param->minPassGQ())
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

Genotype BreakpointGenotyper::genotype(
    const BreakpointGenotyperParameter& param, const std::vector<int32_t>& read_counts_per_allele) const
{
    if (read_counts_per_allele.size() != n_alleles_)
    {
        error(
            "Error: number of read counts and alleles mismatches. %i != %i.", (int)read_counts_per_allele.size(),
            n_alleles_);
    }
    Genotype result;

    // compute adjusted depth
    const double multiplier = (param.read_length - min_overlap_bases_) / (double)param.read_length;
    assert(multiplier > 0);
    const double lambda = param.read_depth * multiplier;
    const int32_t total_num_reads = std::accumulate(read_counts_per_allele.begin(), read_counts_per_allele.end(), 0);
    if (total_num_reads == 0)
    {
        result.filters.insert("NO_READS");
        return result;
    }
    result.num_reads = total_num_reads;

    // compute GL and GT
    double best_gl = -std::numeric_limits<double>::max();

    for (const auto& igt : possible_genotypes)
    {
        const double gl = genotypeLikelihood(lambda, igt, read_counts_per_allele);
        result.gl_name.push_back(igt);
        result.gl.push_back(gl);

        // update GT if GL is better
        if (gl > best_gl)
        {
            best_gl = gl;
            result.gt = igt;
        }
    }

    // compute GQ and set filter
    double sum_gl = 0;
    for (auto l : result.gl)
    {
        sum_gl += exp(l);
    }
    double pr_gt_error = (double)1.0 - exp(best_gl) / sum_gl;
    if (pr_gt_error == 0)
    {
        result.gq = 100;
    }
    else
    {
        double gq_log10 = log10(pr_gt_error);
        if (gq_log10 < -10)
        {
            result.gq = 100;
        }
        else
        {
            result.gq = -10 * gq_log10;
        }
    }
    if (result.gq < min_pass_gq_)
    {
        result.filters.insert("GQ");
    }

    // compute allele fractions
    result.allele_fractions.resize(n_alleles_, 0.0);
    for (unsigned int al = 0; al < n_alleles_; ++al)
    {
        result.allele_fractions[al] = ((double)read_counts_per_allele[al]) / total_num_reads;
    }

    // compute coverage test p value
    double coverage_test_pvalue;
    if (param.use_poisson_depth) // use poisson test for depth (more stringent)
    {
        const poisson_distribution<> poisson_coverage_distribution(lambda);
        coverage_test_pvalue = cdf(poisson_coverage_distribution, total_num_reads);
    }
    else // use normal test for depth (default)
    {
        const normal_distribution<> normal_coverage_distribution(lambda, param.depth_sd);
        coverage_test_pvalue = cdf(normal_coverage_distribution, total_num_reads);
    }

    if (coverage_test_pvalue > 0.5)
    {
        coverage_test_pvalue = 1 - coverage_test_pvalue;
        if (coverage_test_pvalue < coverage_test_cutoff_.first)
        {
            result.filters.insert("BP_DEPTH");
        }
    }
    else
    {
        if (coverage_test_pvalue < coverage_test_cutoff_.second)
        {
            result.filters.insert("BP_DEPTH");
        }
    }
    result.coverage_test_pvalue = coverage_test_pvalue;

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
