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
