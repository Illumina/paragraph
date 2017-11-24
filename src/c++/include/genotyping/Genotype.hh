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

#include "common/Phred.hh"
#include "json/json.h"

namespace genotyping
{
typedef std::vector<uint64_t> GenotypeVector;

struct Genotype
{
    /**
     *  Most likely genotype.
     */
    GenotypeVector gt;

    /**
     *  output string of most likely genotype
     */
    std::string genotypeString() const;
    std::string genotypeString(std::vector<std::string>& allele_names) const;

    /**
     * other most likely genotypes (with same GL)
     */
    std::vector<GenotypeVector> equivalent_gt;

    /**
     *  Log-probability for each possible genotype.
     */
    std::vector<GenotypeVector> gl_name;
    std::vector<double> gl;

    /**
     *  fractions of edge count
     */
    std::vector<double> allele_fractions;

    /**
     *  count of reads (from all edges)
     */
    int num_reads = 0;

    /**
     * filter flag (include PASS)
     */
    std::string filter;

    /**
     * get major information as a single string (simplified)
     */
    explicit operator std::string() const;

    /**
     * recode genotypes using external index
     */
    void recode(std::vector<std::vector<uint64_t>>& alternative_codes);

    /**
     * everything to json.
     * plus convert node index into node name
     */
    Json::Value toJson(std::vector<std::string>& allele_names) const;
};
};
