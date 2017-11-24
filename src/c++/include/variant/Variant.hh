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

#pragma once

#include "variant.pb.h"

#include <graphs/GraphCoordinates.hh>
#include <list>
#include <memory>
#include <string>

#include "RefVar.hh"
#include "common/Phred.hh"
#include "json/json.h"

namespace variant
{

/**
 * Pileup information for a single base
 */
struct PileupData
{
    // DP on forward and reverse strand
    int stranded_DP[2] = { 0, 0 };
    // qual-weighted DP on forward and reverse strand
    float qual_weighted_DP[2] = { 0, 0 };

    /**
     * add an observation
     * @param is_rev true if obs is on reverse strand
     */
    void addObs(bool is_rev = false, int pqual = 60)
    {
        stranded_DP[is_rev ? 1 : 0]++;
        qual_weighted_DP[is_rev ? 1 : 0] += 1.0f - (float)common::phred::phredToErrorProb(pqual);
    }

    PileupData& operator+=(PileupData const& rhs)
    {
        stranded_DP[0] += rhs.stranded_DP[0];
        stranded_DP[1] += rhs.stranded_DP[1];
        qual_weighted_DP[0] += rhs.qual_weighted_DP[0];
        qual_weighted_DP[1] += rhs.qual_weighted_DP[1];
        return *this;
    }

    PileupData& operator-=(PileupData const& rhs)
    {
        stranded_DP[0] = std::max(0, stranded_DP[0] - rhs.stranded_DP[0]);
        stranded_DP[1] = std::max(0, stranded_DP[1] - rhs.stranded_DP[1]);
        qual_weighted_DP[0] = std::max(0.0f, qual_weighted_DP[0] - rhs.qual_weighted_DP[0]);
        qual_weighted_DP[1] = std::max(0.0f, qual_weighted_DP[1] - rhs.qual_weighted_DP[1]);
        return *this;
    }

    PileupData& operator/=(float val)
    {
        stranded_DP[0] /= val;
        stranded_DP[1] /= val;
        qual_weighted_DP[0] /= val;
        qual_weighted_DP[1] /= val;
        return *this;
    }
};

/**
 * A variant candidate list collects variant candidates relative
 * to a reference sequence.
 *
 * It also resolves these candidates into a canonical set and outputs
 * reference coverage for each reference position.
 */
class VariantCandidateList
{
public:
    explicit VariantCandidateList(std::string const& reference);

    VariantCandidateList();
    VariantCandidateList(VariantCandidateList const& rhs);
    VariantCandidateList(VariantCandidateList&& rhs) noexcept;

    ~VariantCandidateList();

    VariantCandidateList& operator=(VariantCandidateList const& rhs);
    VariantCandidateList& operator=(VariantCandidateList&& rhs) noexcept;

    /**
     * @return the reference sequence passed at construction
     */
    std::string const& getReference() const;

    /**
     * Add observation for RefVar allele from a read.
     *
     * @param rv the allele
     * @param is_rev true if observation is on reverse strand
     * @param left_boundary position of preceding variant to the left (if any)
     * @param pqual phred-scaled base quality for rv
     * @return the rightmost position this refvar could be placed
     */
    int64_t addRefVarObservation(RefVar rv, bool is_rev = false, int64_t left_boundary = -1, int pqual = 60);

    /**
     * Get piled up depth summary by position.
     */
    PileupData const& getRefPileup(int pos) const;

    PileupData const& getNonrefPileup(int pos) const;

    /**
     * Get variants
     *
     * @return consolidated list of variants + depths
     */
    std::list<Variant*> getVariants() const;

    /**
     * Append coverage values to JSON value
     *
     * This will append to arrays like this:
     *
     * {
     *    "ref": <ref-matching depth>
     *    "ref:FWD": <...>
     *    "other": non-ref-matching
     *    ...
     * }
     *
     *
     * @param coords coordinate system on the graph to add position information
     * @param node_name which node name to write + use for the Graph coordinates
     * @param coverage coverage JSON. If entries are present already, this function will append
     */
    void
    appendCoverage(graphs::GraphCoordinates const& coords, std::string const& node_name, Json::Value& coverage) const;

private:
    struct VariantCandidateListImpl;
    std::unique_ptr<VariantCandidateListImpl> _impl;
};
};
