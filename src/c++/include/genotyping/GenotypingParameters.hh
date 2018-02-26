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
 * Class for genotyping parameters
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */
#pragma once

#include "common/Error.hh"
#include "genotyping/Genotype.hh"

#include <json/json.h>
#include <map>
#include <string>
#include <vector>

namespace genotyping
{

class GenotypingParameters
{
public:
    /**
     * Given allele names from graph, set default parameters
     * @param _allele_names Vector of allele names
     */
    explicit GenotypingParameters(const std::vector<std::string>& _allele_names);

    /**
     * Use parameters from external JSON to override default parameters
     * @param param_json Input JSON value for parameters
     */
    void setFromJson(Json::Value& param_json);

    /**
     * getters for genotypers
     */
    unsigned int ploidy() { return ploidy_; }

    unsigned int numAlleles() { return num_alleles; }

    unsigned int minOverlapBases() { return min_overlap_bases; }

    std::vector<double>& alleleErrorRates() { return allele_error_rates; }

    std::vector<double>& hetHaplotypeFractions() { return het_haplotype_fractions; }

    std::map<GenotypeVector, double>& genotypeFractions() { return genotype_fractions; }

    double otherAlleleErrorRate() { return other_allele_error_rate; }

    double otherHetHaplotypeFraction() { return other_het_haplotype_fraction; }

    std::vector<GenotypeVector>& possibleGenotypes() { return possible_genotypes; };

private:
    /**
     * Set the vector of all possible unphased genotypes from alleles and ploidy.
     */
    void setPossibleGenotypes();

    /**
     * Build allele index: allele index in JSON -> allele index in original allele names
     * @param param_json JSON array of allele names
     */
    std::vector<int> alleleNameConversionIndex(Json::Value& param_json);

    /**
     * Assign parameters from JSON to EMPTY parameter vectors
     * @param param_json JSON value containing the parameters
     * @param conversion_index Conversion to allele index for alleles in JSON
     */

    void setAlleleErrorRate(Json::Value& param_json, std::vector<int>& conversion_index);

    void setHetHaplotypeFractions(Json::Value& param_json, std::vector<int>& conversion_index);

    void setGenotypeFractions(Json::Value& param_json, std::vector<int>& conversion_index);

    /**
     * ploidy of variants
     */
    unsigned int ploidy_;

    /**
     * number of alleles. Same as length of allele_names
     */
    const unsigned int num_alleles;

    /**
     * allele names from graph input
     */
    const std::vector<std::string> allele_names;

    /**
     * number of offset bases for poission lambda calculation
     */
    unsigned int min_overlap_bases;

    /**
     * all possible genotypes
     */
    std::vector<GenotypeVector> possible_genotypes;

    /**
     * error rate for alleles used in breakpint genotyping
     */
    std::vector<double> allele_error_rates;

    /**
     * fraction of non-reference allele in a REF/ALT genotype
     * a representation of read extraction efficiency
     */
    std::vector<double> het_haplotype_fractions;

    /**
     * fraction of genenotypes (the prior)
     */
    std::map<GenotypeVector, double> genotype_fractions;

    /**
     * parameters of other no-show alleles
     */
    std::string reference_allele;

    double reference_allele_error_rate;

    double other_allele_error_rate;

    double other_het_haplotype_fraction;

    double other_genotype_fraction;
};
};