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
 * \brief Genotyper interface
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <set>
#include <utility>

#include "json/json.h"

namespace genotyping
{
typedef std::vector<uint64_t> GenotypeVector;

/**
 * genotype vector from given string such as 0/1
 */
GenotypeVector genotypeVectorFromString(std::string& gt_str);

struct Genotype
{
    Genotype() = default;
    Genotype(Genotype const&) = default;
    Genotype(Genotype&&) = default;
    explicit Genotype(GenotypeVector rhs)
        : gt(std::move(rhs))
    {
    }
    virtual ~Genotype() = default;

    Genotype& operator=(Genotype const&) = default;
    Genotype& operator=(Genotype&&) = default;

    /**
     *  Most likely genotype.
     */
    GenotypeVector gt;

    /**
     *  Log-likelihood for each possible genotype.
     */
    std::vector<GenotypeVector> gl_name;
    std::vector<double> gl;

    /**
     * Genotype quality
     * -10log10(1 - gl/sum(gl))
     */
    int gq = -1;

    /**
     *  fractions of edge count
     */
    std::vector<double> allele_fractions;

    /**
     * tail probability of coverage (all alleles)
     */
    double coverage_test_pvalue = -1;

    /**
     *  count of reads (from all edges)
     */
    int num_reads = 0;

    /**
     * filter description (include PASS)
     */
    std::set<std::string> filters;

    /**
     * relabel index-based genotypes with given new labels
     */
    void relabel(std::vector<uint64_t> const& new_labels);

    /**
     * everything to json.
     * plus convert node index into node name
     */
    Json::Value toJson(std::vector<std::string> const& allele_names) const;

    /**
     *  output string of most likely genotype
     */
    std::string toString(const std::vector<std::string>* p_allele_names = NULL) const;

    /**
     * output string of filter
     */
    std::string filterString() const;

    /**
     * get major information as a single string (simplified)
     */
    explicit operator std::string() const;
};
};
