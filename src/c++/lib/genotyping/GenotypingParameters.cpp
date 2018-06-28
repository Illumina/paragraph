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
 * Interface for genotyping parameters
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/GenotypingParameters.hh"

#include <algorithm>

using std::string;
using std::vector;

namespace genotyping
{
GenotypingParameters::GenotypingParameters(const vector<string>& _allele_names, unsigned int ploidy)
    : ploidy_(ploidy)
    , num_alleles(static_cast<const unsigned int>(_allele_names.size()))
    , coverage_test_cutoff(-1.0)
    , allele_names(_allele_names)
    , min_overlap_bases(16)
    , reference_allele("REF")
    , reference_allele_error_rate(0.05)
    , other_allele_error_rate(0.05)
    , other_het_haplotype_fraction(0.5)
    , other_genotype_fraction(1)
{
    setPossibleGenotypes();
}

void GenotypingParameters::setPossibleGenotypes()
{
    vector<GenotypeVector> gts;
    if (num_alleles)
    {
        // generate all possible GTs (as described in the VCF SPEC)
        const std::function<void(unsigned int, unsigned int, std::vector<uint64_t>)> makeGenotypes
            = [&makeGenotypes, &gts](unsigned int p, unsigned int n, vector<uint64_t> suffix) {
                  for (unsigned int a = 0; a <= n; ++a)
                  {
                      if (p == 1)
                      {
                          auto new_suffix = suffix;
                          new_suffix.insert(new_suffix.begin(), a);
                          gts.push_back(new_suffix);
                      }
                      else if (p > 1)
                      {
                          auto new_suffix = suffix;
                          new_suffix.insert(new_suffix.begin(), a);
                          makeGenotypes(p - 1, a, new_suffix);
                      }
                  }
              };
        makeGenotypes(ploidy_, num_alleles - 1, {});
    }
    possible_genotypes = std::move(gts);
}

void GenotypingParameters::setFromJson(Json::Value& param_json)
{
    auto logger = LOG();
    // simple statistic fields
    bool uniform_het_haplotype_fraction = false;
    for (auto& key : param_json.getMemberNames())
    {
        auto& field = param_json[key];
        if (key == "min_overlap_bases")
        {
            min_overlap_bases = field.asUInt();
        }
        if (key == "coverage_test_cutoff")
        {
            coverage_test_cutoff = field.asDouble();
        }
        else if (key == "reference_allele")
        {
            reference_allele = field.asString();
        }
        else if (key == "reference_allele_error_rate")
        {
            reference_allele_error_rate = field.asDouble();
        }
        else if (key == "other_allele_error_rate")
        {
            other_allele_error_rate = field.asDouble();
        }
        else if (key == "het_haplotype_fraction")
        {
            if (field.asString()[0] == '[')
            {
                other_het_haplotype_fraction = field.asDouble();
                uniform_het_haplotype_fraction = true;
            }
        }
        else if (key == "other_genotype_fraction")
        {
            other_genotype_fraction = field.asDouble();
        }
        else if (key == "ploidy")
        {
            ploidy_ = (unsigned int)field.asInt();
        }
    }

    if (param_json.isMember("allele_error_rates")
        || (param_json.isMember("het_haplotype_fractions") && !uniform_het_haplotype_fraction)
        || param_json.isMember("genotype_fractions"))
    {
        if (!param_json.isMember("allele_names"))
        {
            error("Error: with allele_error_rates/het_haplotype_fractions/genotype_fractions specified in JSON, "
                  "allele_names must be specified as well.");
        }
        auto conversion_index = alleleNameConversionIndex(param_json["allele_names"]);

        if (std::all_of(conversion_index.begin(), conversion_index.end(), [](int x) { return x == -1; }))
        {
            logger->warn("None of the allele names in JSON match those in graph. Parameters in "
                         "allele_error_rates/het_haplotype_fractions/genotype_fractions are not used.");
        }
        else
        {
            if (param_json.isMember("allele_error_rates"))
            {
                setAlleleErrorRate(param_json["allele_error_rates"], conversion_index);
            }

            if (param_json.isMember("het_haplotype_fractions") && !uniform_het_haplotype_fraction)
            {
                setHetHaplotypeFractions(param_json["het_haplotype_fractions"], conversion_index);
            }

            if (param_json.isMember("genotype_fractions"))
            {
                setGenotypeFractions(param_json["genotype_fractions"], conversion_index);
            }
        }
    }
}

vector<int> GenotypingParameters::alleleNameConversionIndex(Json::Value& param_json)
{
    vector<int> conversion_index;
    for (auto& key : param_json)
    {
        auto it = std::find(allele_names.begin(), allele_names.end(), key.asString());
        if (it != allele_names.end())
        {
            const int dist = static_cast<int>(it - allele_names.begin());
            conversion_index.push_back(dist);
        }
        else
        {
            conversion_index.push_back(-1);
        }
    }
    return conversion_index;
}

void GenotypingParameters::setAlleleErrorRate(Json::Value& param_json, vector<int>& conversion_index)
{
    allele_error_rates.resize(num_alleles, other_allele_error_rate);
    auto it_ref = std::find(allele_names.begin(), allele_names.end(), reference_allele);
    if (it_ref != allele_names.end())
    {
        const size_t ref_index = (size_t)(it_ref - allele_names.begin());
        allele_error_rates[ref_index] = reference_allele_error_rate;
    }
    size_t index = 0;
    for (auto& v : param_json)
    {
        auto new_index = conversion_index[index];
        if (new_index != -1)
        {
            allele_error_rates[new_index] = v.asDouble();
        }
        index++;
    }
}

void GenotypingParameters::setHetHaplotypeFractions(Json::Value& param_json, vector<int>& conversion_index)
{
    het_haplotype_fractions.resize(num_alleles, other_het_haplotype_fraction);
    size_t index = 0;
    for (auto& v : param_json)
    {
        auto new_index = conversion_index[index];
        if (new_index != -1)
        {
            het_haplotype_fractions[new_index] = v.asDouble();
        }
        index++;
    }
}

void GenotypingParameters::setGenotypeFractions(Json::Value& param_json, vector<int>& conversion_index)
{
    for (auto& gt_str : param_json.getMemberNames())
    {
        GenotypeVector gv = genotypeVectorFromString(gt_str);
        if (gv.empty())
        {
            error("Error: Empty or illegal genotype in parameter JSON: %s", gt_str.c_str());
        }

        GenotypeVector new_gt;
        for (auto& g : gv)
        {
            const auto new_index = (size_t)conversion_index[g];
            if (new_index != (size_t)-1)
            {
                new_gt.push_back(new_index);
            }
            else
            {
                break;
            }
        }
        if (new_gt.size() == ploidy_)
        {
            genotype_fractions[new_gt] = param_json["genotype_fractions"][gt_str].asDouble();
        }
    }

    // complete the map
    for (auto& gt : possible_genotypes)
    {
        if (genotype_fractions.find(gt) == genotype_fractions.end())
        {
            genotype_fractions[gt] = other_genotype_fraction;
        }
    }
}
}
