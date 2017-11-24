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
#include <limits>
#include <math.h>

using boost::math::poisson_distribution;
using std::vector;

namespace genotyping
{

Genotype BreakpointGenotyper::genotype(
    const vector<int32_t>& read_counts, std::vector<double> haplotypeReadFraction,
    std::vector<double> genotypeErrorRate) const
{
    if (read_counts.size() < 2)
    {
        error("Error: haploid genotypes are not supported yet.");
    }
    if (!haplotypeReadFraction.empty())
    {
        if (read_counts.size() != haplotypeReadFraction.size())
        {
            error(
                "Error: use-specified haplotypeReadFraction should have size equal to the number of haplotypes. Error "
                "at %d, %d.",
                read_counts.size(), haplotypeReadFraction.size());
        }
        for (double mu : haplotypeReadFraction)
        {
            if (mu <= 0 || mu >= 1)
            {
                error("Error: use-specified haplotypeReadFraction should be between 0~1.");
            }
        }
    }
    if (!genotypeErrorRate.empty())
    {
        uint64_t num_haplotypes = (uint64_t)read_counts.size();
        uint64_t num_genotypes = num_haplotypes * (num_haplotypes + 1) / 2;
        if ((uint64_t)genotypeErrorRate.size() != num_genotypes)
        {
            error("Error: use-specified genotypeErrorRate should have size equal to the number of possible "
                  "genotypes.");
        }
    }
    return genotypeDiploidBreakpoint(read_counts, haplotypeReadFraction, genotypeErrorRate);
}

Genotype BreakpointGenotyper::genotypeDiploidBreakpoint(
    const vector<int32_t>& read_counts, std::vector<double> haplotypeReadFraction,
    std::vector<double> genotypeErrorRate) const
{
    const int32_t num_haplotypes = (int32_t)read_counts.size();
    const int32_t num_genotypes = num_haplotypes * (num_haplotypes + 1) / 2;
    if (haplotypeReadFraction.empty())
    {
        haplotypeReadFraction.resize(num_haplotypes, 0.5);
    }
    if (genotypeErrorRate.empty())
    {
        genotypeErrorRate.resize(num_genotypes, genotype_error_rate_);
    }
    const double phi = (double)1 / num_genotypes; // Fixed prior for each genotype.
    const double multiplier = (read_length_ - min_overlap_bases_) / (double)read_length_;
    const double lambda = read_depth_ * multiplier;

    int32_t total_num_reads = 0;
    for (int32_t count : read_counts)
    {
        total_num_reads += count;
    }

    Genotype genotype;
    genotype.num_reads = total_num_reads;
    double max_log_genotype_prob = std::numeric_limits<double>::lowest();
    int32_t most_likely_haplotype1_index = -1;
    int32_t most_likely_haplotype2_index = -1;

    int genotype_index = 0;
    for (int32_t haplotype1_index = 0; haplotype1_index != num_haplotypes; haplotype1_index++)
    {
        int32_t num_haplotype1_reads = read_counts[haplotype1_index];
        for (int32_t haplotype2_index = haplotype1_index; haplotype2_index < num_haplotypes; haplotype2_index++)
        {
            double log_genotype_prob = 0;
            if (haplotype2_index == haplotype1_index)
            {
                double lambda_error = lambda * genotypeErrorRate[genotype_index];
                double lambda_mapped = lambda * (1 - genotypeErrorRate[genotype_index]);
                if (lambda_mapped == 0)
                {
                    lambda_mapped = std::numeric_limits<double>::lowest() * (-1);
                }
                if (lambda_error == 0)
                {
                    lambda_error = std::numeric_limits<double>::lowest() * (-1);
                }
                const poisson_distribution<> hom_hap1_distribution(lambda_mapped);
                const poisson_distribution<> hom_hap1_error_distribution(lambda_error);
                int32_t num_spurious_reads = total_num_reads - num_haplotype1_reads;
                log_genotype_prob = log(phi) + log(pdf(hom_hap1_distribution, num_haplotype1_reads))
                    + log(pdf(hom_hap1_error_distribution, num_spurious_reads));
            }
            else
            {
                double lambda_hap1 = lambda * haplotypeReadFraction[haplotype1_index];
                double lambda_hap2 = lambda * haplotypeReadFraction[haplotype2_index];
                double lambda_het_error = lambda * genotypeErrorRate[genotype_index];
                if (lambda_hap1 == 0)
                {
                    lambda_hap1 = std::numeric_limits<double>::lowest() * (-1);
                }
                if (lambda_hap2 == 0)
                {
                    lambda_hap2 = std::numeric_limits<double>::lowest() * (-1);
                }
                if (lambda_het_error == 0)
                {
                    lambda_het_error = std::numeric_limits<double>::lowest() * (-1);
                }
                const poisson_distribution<> single_hap1_distribution(lambda_hap1);
                const poisson_distribution<> single_hap2_distribution(lambda_hap2);
                const poisson_distribution<> het_error_distribution(lambda_het_error);
                int32_t num_haplotype2_reads = read_counts[haplotype2_index];
                int32_t num_spurious_reads = total_num_reads - num_haplotype1_reads - num_haplotype2_reads;
                log_genotype_prob = log(phi) + log(pdf(single_hap1_distribution, num_haplotype1_reads))
                    + log(pdf(single_hap2_distribution, num_haplotype2_reads))
                    + log(pdf(het_error_distribution, num_spurious_reads));
            }
            genotype.gl_name.push_back({ (uint64_t)haplotype1_index, (uint64_t)haplotype2_index });
            if (isinf(log_genotype_prob))
            {
                log_genotype_prob = -std::numeric_limits<double>::max();
            }
            genotype.gl.push_back(log_genotype_prob);

            if (max_log_genotype_prob < log_genotype_prob)
            {
                max_log_genotype_prob = log_genotype_prob;
                most_likely_haplotype1_index = haplotype1_index;
                most_likely_haplotype2_index = haplotype2_index;
            }
            genotype_index++;
        }
    }

    if (genotype_index > 0)
    {
        if (most_likely_haplotype1_index == -1) // all -Inf
        {
            return genotype;
        }
        genotype.gt = { (uint64_t)most_likely_haplotype1_index, (uint64_t)most_likely_haplotype2_index };
        int32_t sum_count = 0;
        for (auto r_count : read_counts)
        {
            sum_count += r_count;
        }
        genotype.allele_fractions.resize(num_haplotypes, 0);
        for (auto i = 0; i < num_haplotypes; i++)
        {
            genotype.allele_fractions[i] = (double)read_counts[i] / sum_count;
        }
    }
    return genotype;
}
}