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
 * @return final GT
 */
Genotype combinedGenotype(
    GenotypeSet const& genotypes, const BreakpointGenotyper* p_genotyper = NULL, double depth = 0, int read_length = 0);

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
    GenotypeSet const& genotypes, bool use_pass_only, const BreakpointGenotyper* p_genotyper, double depth,
    int read_length);
};
