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
 * Simple Poisson model for genotyping a single diploid breakpoint from a single sample
 *
 *  A full description of the model can be found in doc/graph-models.md
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"

namespace genotyping
{
class BreakpointGenotyper
{
public:
    BreakpointGenotyper(
        const double read_depth, const int32_t read_length, const double genotype_error_rate = 0.01,
        const int32_t min_overlap_bases = 16)
        : read_depth_(read_depth)
        , read_length_(read_length)
        , genotype_error_rate_(genotype_error_rate)
        , min_overlap_bases_(min_overlap_bases){};

    /**
     * public function to do breakpoint genotyping from read counts on edges
     * after sanity check, invode genotypeDiploidBrekpoint
     */
    Genotype genotype(
        const std::vector<int32_t>& read_counts, std::vector<double> haplotypeReadFraction = {},
        std::vector<double> genotypeErrorRate = {}) const;

private:
    /**
     * core funciton for doing breakpoint genotyping
     */
    Genotype genotypeDiploidBreakpoint(
        const std::vector<int32_t>& read_counts, std::vector<double> haplotypeReadFraction,
        std::vector<double> genotypeErrorRate) const;

    /**
     * necessary bam/cram stats
     */
    const double read_depth_;
    const int32_t read_length_;

    // genotyping specific stats
    const double genotype_error_rate_;
    const int32_t min_overlap_bases_; // Minimum number of bases that a high-confidence alignment must overlap a node.
};
};