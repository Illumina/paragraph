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
 * Population statistics for genotype sets
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/PopulationStatistics.hh"

#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <numeric>

using std::vector;

namespace genotyping
{
PopulationStatistics::PopulationStatistics(GenotypeSet const& genotypes)
{
    num_valid_samples = 0;
    num_total_samples = static_cast<int>(genotypes.size());

    for (auto const& genotype : genotypes)
    {
        if (genotype.gt.empty())
        {
            continue;
        }
        num_valid_samples++;
        if (genotype_counts.find(genotype.gt) == genotype_counts.end())
        {
            genotype_counts[genotype.gt] = 1;
        }
        else
        {
            genotype_counts[genotype.gt]++;
        }

        for (auto const& gt : genotype.gt)
        {
            if (allele_counts.size() <= gt)
            {
                allele_counts.resize(gt + 1, 0);
            }
            allele_counts[gt]++;
        }
    }
}

Json::Value PopulationStatistics::toJson() const
{
    double hwe_p_chisq = getChisqPvalue();
    double hwe_p_fisher = -1;
    if (needFisherExactHWE())
    {
        hwe_p_fisher = getFisherExactPvalue();
    }
    double call_rate = getCallrate();

    Json::Value json_result;
    json_result["hwe"] = hwe_p_chisq;
    if (hwe_p_fisher == -1)
    {
        json_result["hwe_fisher"] = "";
    }
    else
    {
        json_result["hwe_fisher"] = hwe_p_fisher;
    }
    json_result["call_rate"] = call_rate;

    auto allele_frequencies = getAlleleFrequencies();
    json_result["allele_frequencies"] = Json::arrayValue;
    for (auto& a : allele_frequencies)
    {
        json_result["allele_frequencies"].append(a);
    }

    return json_result;
}

double PopulationStatistics::getChisqPvalue() const
{
    double chisq_val = 0;
    for (auto& gv : genotype_counts)
    {
        if (gv.first.size() != 2)
        {
            continue;
        }
        uint64_t h1 = gv.first[0];
        uint64_t h2 = gv.first[1];
        if (allele_counts[h1] == 0 || allele_counts[h2] == 0) // skip unobserved alleles
        {
            continue;
        }
        double e_count;
        if (h1 == h2)
        {
            e_count = ((double)allele_counts[h1] / num_valid_samples / 2)
                * ((double)allele_counts[h1] / num_valid_samples / 2) * num_valid_samples;
        }
        else
        {
            e_count = 2 * ((double)allele_counts[h1] / num_valid_samples / 2)
                * ((double)allele_counts[h2] / num_valid_samples / 2) * num_valid_samples;
        }

        double diff = e_count - gv.second;
        double norm_diff_square = (double)(diff * diff) / e_count;
        chisq_val += norm_diff_square;
    }
    boost::math::chi_squared chisq_distribution(1);
    double hwe_pval = 1 - boost::math::cdf(chisq_distribution, chisq_val);
    return hwe_pval;
}

/**
 * @return True if need to use Fisher's exact test for HWE P
 * For multi-alleles, always use chisq
 * For bi-allelic:
 *      if
 *          N<=30, or count of the rarest genotype <= 20, or rarest expected count <= 20, use fisher's exact
 *      else
 *          use chisq
 */
bool PopulationStatistics::needFisherExactHWE() const
{
    int num_observed_alleles = 0;
    for (auto a : allele_counts)
    {
        if (a > 0)
        {
            num_observed_alleles++;
        }
    }
    if (num_observed_alleles <= 1)
    {
        return false;
    }
    if (num_observed_alleles > 2)
    {
        return false;
    }

    if (num_valid_samples <= 30)
    {
        return true;
    }

    for (auto& g : genotype_counts)
    {
        if (g.second > 0 && g.second <= 20)
        {
            return true;
        }
    }

    auto minor_allele_index = minNonZeroAlleleIndex();
    double minor_allele_freq = (double)allele_counts[minor_allele_index] / 2 / num_valid_samples;
    if (minor_allele_freq * minor_allele_freq * num_valid_samples <= 20)
    {
        return true;
    }
    return false;
}

/**
 * from JE Wigginton 2005 AJHG
 */
double PopulationStatistics::getFisherExactPvalue() const
{
    auto minor_allele_index = minNonZeroAlleleIndex();
    auto p_major = std::max_element(allele_counts.begin(), allele_counts.end());
    int minor_allele_count = allele_counts[minor_allele_index];
    int major_allele_count = *p_major;

    GenotypeVector het_gv;
    het_gv.push_back(static_cast<uint64_t>(p_major - allele_counts.begin()));
    het_gv.push_back(static_cast<uint64_t>(minor_allele_index));
    std::sort(het_gv.begin(), het_gv.end());
    int observed_num_het = 0;
    for (auto& g : genotype_counts)
    {
        if (g.first.size() != 2)
        {
            continue;
        }
        if (g.first[0] == het_gv[0] && g.first[1] == het_gv[1])
        {
            observed_num_het = g.second;
            break;
        }
    }

    int num_expect_het = std::round(
        2 * ((double)minor_allele_count / num_valid_samples / 2) * ((double)major_allele_count / num_valid_samples / 2)
        * num_valid_samples);

    std::vector<double> scaled_pvals;
    double observe_scaled_pval = -1;

    // for num_het > expected het
    int prev_num_ref_hom = (minor_allele_count - num_expect_het) / 2;
    int prev_num_alt_hom = num_valid_samples - prev_num_ref_hom - num_expect_het;
    double prev_scaled_pval = 1;
    for (int num_het = num_expect_het; num_het <= minor_allele_count; num_het += 2)
    {
        if (num_het == num_expect_het)
        {
            scaled_pvals.push_back(1);
            continue;
        }
        int prev_num_het = num_het - 2;
        double iscale
            = prev_scaled_pval * (4 * prev_num_ref_hom * prev_num_alt_hom) / ((prev_num_het + 2) * (prev_num_het + 1));
        scaled_pvals.push_back(iscale);
        prev_scaled_pval = iscale;
        prev_num_ref_hom--;
        prev_num_alt_hom--;
        if (observe_scaled_pval == -1 && num_het == observed_num_het)
        {
            observe_scaled_pval = iscale;
        }
    }

    // for num_het < expected het
    prev_num_ref_hom = (minor_allele_count - num_expect_het) / 2;
    prev_num_alt_hom = num_valid_samples - prev_num_ref_hom - num_expect_het;
    prev_scaled_pval = 1;
    for (int num_het = num_expect_het; num_het >= 0; num_het -= 2)
    {
        if (num_het == num_expect_het)
        {
            continue;
        }
        int prev_num_het = num_het + 2;
        double iscale = prev_scaled_pval / 4 * prev_num_het / (prev_num_ref_hom + 1) * (prev_num_het - 1)
            / (prev_num_alt_hom + 1);
        scaled_pvals.push_back(iscale);
        prev_scaled_pval = iscale;
        prev_num_ref_hom++;
        prev_num_alt_hom++;
        if (observe_scaled_pval == -1 && num_het == observed_num_het)
        {
            observe_scaled_pval = iscale;
        }
    }

    // calculate HWE
    double hwe_scale_sum = 0;
    for (auto s : scaled_pvals)
    {
        if (s <= observe_scaled_pval)
        {
            hwe_scale_sum += s;
        }
    }
    double pval = hwe_scale_sum / std::accumulate(scaled_pvals.begin(), scaled_pvals.end(), (double)0);
    return pval;
}

/**
 * @return allele frequencies
 */
std::vector<double> PopulationStatistics::getAlleleFrequencies() const
{
    std::vector<double> result;
    const uint32_t sum = std::accumulate(allele_counts.begin(), allele_counts.end(), (uint32_t)0);
    result.reserve(allele_counts.size());
    for (auto ac : allele_counts)
    {
        if (sum > 0)
        {
            result.push_back(((double)ac) / sum);
        }
        else
        {
            result.push_back(0);
        }
    }
    return result;
}

/**
 * @return iterator to the lowest allele non-zero minor allele count
 */
size_t PopulationStatistics::minNonZeroAlleleIndex() const
{
    auto p_minor = std::min_element(allele_counts.begin(), allele_counts.end());
    if (*p_minor > 0)
    {
        return static_cast<size_t>(p_minor - allele_counts.begin());
    }
    p_minor = std::max_element(allele_counts.begin(), allele_counts.end());
    if (*p_minor == 0)
    {
        return static_cast<size_t>(0);
    }
    for (auto it = allele_counts.begin(); it != allele_counts.end(); it++)
    {
        if (*it < *p_minor)
        {
            p_minor = it;
        }
    }
    return static_cast<size_t>(p_minor - allele_counts.begin());
}
}