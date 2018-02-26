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
 * Container for multiple genotypes for one locus and for population-level statistics
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/GenotypeSet.hh"

#include <algorithm>

namespace genotyping
{
/**
 * Add a genotype for a particular sample
 * @param allele_names names of all alleles in GT
 * @param gt genotype of the sample
 * @return sample index
 */
size_t GenotypeSet::add(std::vector<std::string> const& allele_names, Genotype const& gt)
{
    genotypes.push_back(gt);
    Genotype& remapped_gt = genotypes.back();

    std::vector<uint64_t> gt_remapping(allele_names.size());
    size_t j = 0;
    for (auto const& a : allele_names)
    {
        auto a_it = std::find(merged_allele_names.begin(), merged_allele_names.end(), a);
        if (a_it == merged_allele_names.end())
        {
            gt_remapping[j] = merged_allele_names.size();
            merged_allele_names.push_back(a);
        }
        else
        {
            gt_remapping[j] = static_cast<unsigned long>(a_it - merged_allele_names.begin());
        }
        ++j;
    }
    remapped_gt.relabel(gt_remapping);
    return genotypes.size() - 1;
}
}
