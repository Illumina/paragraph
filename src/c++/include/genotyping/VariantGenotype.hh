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
 * \Storage for breakpoint-voted genotypes
 *
 * \author Sai Chen
 * \email schen6@illumina.com
 *
 */

#pragma once

#include "common/Phred.hh"
#include "genotyping/Genotype.hh"
#include "json/json.h"

namespace genotyping
{
// store vote information for breakpoints
struct genotypeVoteInfo
{
    int num_support;
    int num_unsupport;
    genotypeVoteInfo()
        : num_support(0)
        , num_unsupport(0)
    {
    }
};

class VariantGenotype
{
public:
    explicit VariantGenotype(GenotypeVector external_gt = {})
        : gt_(external_gt){};

    // simple getters
    std::string filter() const { return filter_; }
    bool empty() const { return gt_.empty(); }
    bool zeroCount() const { return filter_ == "ZERO_COUNT"; }
    bool pass() const { return filter_ == "PASS"; }
    GenotypeVector getGenotype() { return gt_; }

    // add information from one exteranl breakpoint
    void addInfoFromSingleBreakpoint(const Genotype& breakpoint_genotype);

    // perform genotyping. Should be called when all breakpoints are added
    void genotype(GenotypeVector reference_homozygote = { 0, 0 });

    /**
     * output genotype information only
     */
    std::string genotypeString() const;
    // convert index to allele names in genotype output
    std::string genotypeString(std::vector<std::string>& allele_names) const;

    /**
     * export everything to JSON
     * Pluse convert index to name string
     */
    Json::Value toJson(std::vector<std::string> allele_names) const;

    /**
     * export key result into a single string
     */
    explicit operator std::string() const;

private:
    /**
     * if 2 best ones and one is reference homozygous, report the non-reference one
     * otherwise return empty genotypes
     */
    GenotypeVector genotypeFromDecisionTree(
        GenotypeVector& reference_homozygote, std::vector<std::vector<uint64_t>>& best_genotype_vec);

    // Most likely genotype.
    GenotypeVector gt_;

    // other most likely genotypes (with same GL)
    std::vector<GenotypeVector> equivalent_gt_;

    // # breakpoints with missing read counts
    int num_missing_ = 0;

    // other possible genotypes and their vote info
    std::map<GenotypeVector, genotypeVoteInfo> gt_vote_info_;

    // possible filter flag (include PASS)
    std::string filter_;
};
};
