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
 * Class to combine breakpoint genotypes to variant genotypes
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <list>

#include "genotyping/BreakpointGenotyper.hh"
#include "genotyping/Genotype.hh"
#include "genotyping/GenotypeSet.hh"

namespace genotyping
{

/**
 * Function to combine breakpoint genotypes to variant genotypes with voting scheme:
 *
 * If all PASS GTs agree, will use consensus
 * If all FAIL GTs agree, will use consensus and tag filter as ALL_BAD_BP
 * If PASS GTs disagree, use total counts from PASS GTs to re-genotype and tag as CONFLICT
 * If all FAIL GTs and disagree, use all total counts to re-genotype and tag as ALL_BAD_BP + CONFLICT
 *
 * @param genotypes Set of genotypes
 * @param genotyper Breakpoint genotyper for additional genotyping
 * @param depth Depth of this sample
 * @param read_length Read length of this sample
 * @param depth_sd Depth standard deviation
 * @return final GT
 */
Genotype combinedGenotype(
    GenotypeSet const& genotypes, const BreakpointGenotyperParameter* b_param = NULL,
    const BreakpointGenotyper* p_genotyper = NULL);

/**
 * Count number of unique genotypes in the genotype set
 * @param genotypes Set of genotypes
 * @param pass_only If true, only use pass genotypes
 * @return Count
 */
size_t countUniqGenotypes(GenotypeSet const& genotypes, bool pass_only);

/**
 * @param genotypes Set of genotypes
 * @param pass_only If ture, only use pass genotypes
 * @return Final GT
 */
Genotype reportConsensusGenotypes(GenotypeSet const& genotypes, bool pass_only);

/**
 * @param genotypes Set of genotypes
 * @param pass_only If ture, only use pass genotypes
 * @param genotyper Breakpoint genotyper class
 * @param depth Depth of this sample
 * @param read_length Read length of this sample
 * @return Final GT from total counts
 */
Genotype genotypeByTotalCounts(
    GenotypeSet const& genotypes, bool use_pass_only, const BreakpointGenotyper* p_genotyper,
    const BreakpointGenotyperParameter* b_param);
};
