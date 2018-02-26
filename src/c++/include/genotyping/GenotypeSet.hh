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
