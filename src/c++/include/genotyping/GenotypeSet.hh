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
 * Container for multiple genotypes
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "BreakpointStatistics.hh"
#include "Genotype.hh"

#include <string>
#include <vector>

namespace genotyping
{
/**
 * Container for genotypes in multiple samples
 */
class GenotypeSet
{
public:
    /**
     * Add a genotype for a particular sample
     * @param allele_names names of all alleles in GT
     * @param gt genotype of the sample
     * @return index
     */
    size_t add(std::vector<std::string> const& allele_names, Genotype const& gt);

    /**
     * Get the genotype for a sample
     * @param sample_index index of the sample
     * @return  genotype in the sample
     */
    Genotype const& operator[](size_t sample_index) const { return genotypes[sample_index]; }

    /**
     * @return The number of samples we have genotypes for
     */
    size_t size() const { return genotypes.size(); }

    /**
     * make range-based for loop compatible
     */
    typedef std::vector<Genotype>::const_iterator const_iterator;
    const_iterator begin() const { return genotypes.cbegin(); }
    const_iterator end() const { return genotypes.cend(); }

    /**
     * Get the names of all alleles
     * @return vector of allele names
     */
    std::vector<std::string> getAlleleNames() const { return merged_allele_names; }

private:
    std::vector<Genotype> genotypes;
    std::vector<std::string> merged_allele_names;
};
}
