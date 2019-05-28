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
