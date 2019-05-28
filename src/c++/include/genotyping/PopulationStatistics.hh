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