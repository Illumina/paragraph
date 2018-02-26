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

#pragma once

#include "genotyping/GenotypeSet.hh"
#include "json/json.h"

#include <map>
#include <vector>

namespace genotyping
{
class PopulationStatistics
{
public:
    explicit PopulationStatistics(const GenotypeSet& genotypes);

    /**
     * perform all calculations and output everything to json.
     */
    Json::Value toJson() const;

    /**
     * @return pval Chisq HWE P value for this event
     * default method of HWE P
     */
    double getChisqPvalue() const;

    /**
     * @return True if need to use Fisher's exact test to calculate HWE P
     */
    bool needFisherExactHWE() const;

    /**
     * @return Fisher's exact HWE P value for this event
     */
    double getFisherExactPvalue() const;

    /**
     * @return fraction of samples with unfiltered calls
     */
    double getCallrate() const { return (double)num_valid_samples / num_total_samples; };

    /**
     * @return allele frequencies
     */
    std::vector<double> getAlleleFrequencies() const;

    /**
     *  @return allele counts
     */
    std::vector<uint32_t> const& alleleCounts() const { return allele_counts; };

private:
    /**
     * @return iterator to the lowest allele non-zero minor allele count
     */
    size_t minNonZeroAlleleIndex() const;

    std::map<GenotypeVector, int> genotype_counts; // genotype -> #samples with this genotype
    std::vector<uint32_t> allele_counts; // allele -> count of this allele in all samples
    int num_total_samples;
    int num_valid_samples;
};
};