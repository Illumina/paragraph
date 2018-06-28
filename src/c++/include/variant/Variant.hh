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

#include <list>
#include <memory>
#include <string>

#include "RefVar.hh"
#include "common/Phred.hh"
#include "graphcore/GraphCoordinates.hh"
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

class Variant
{
public:
    Variant() = default;
    Variant(Variant const&) = default;
    Variant(Variant&&) = default;
    Variant& operator=(Variant const&) = default;
    Variant& operator=(Variant&&) = default;
    ~Variant() = default;

    std::string const& contig() { return contig_; }
    void set_contig(std::string value) { contig_ = value; }

    int32_t start() { return start_; }
    void set_start(int32_t value) { start_ = value; }

    int32_t end() const { return end_; }
    void set_end(int32_t value) { end_ = value; }

    // STR shifting: left / rightmost positions
    int32_t leftmost() const { return leftmost_; }
    void set_leftmost(int32_t value) { leftmost_ = value; }
    int32_t rightmost() const { return rightmost_; }
    void set_rightmost(int32_t value) { rightmost_ = value; }

    std::string const& alt() const { return alt_; }
    void set_alt(std::string value) { alt_ = std::move(value); }

    // allele depths for reference
    int32_t adr_forward() const { return adr_forward_; }
    void set_adr_forward(int32_t value) { adr_forward_ = value; }
    int32_t adr_backward() const { return adr_backward_; }
    void set_adr_backward(int32_t value) { adr_backward_ = value; }

    // allele depths for alt
    int32_t ada_forward() const { return ada_forward_; }
    void set_ada_forward(int32_t value) { ada_forward_ = value; }
    int32_t ada_backward() const { return ada_backward_; }
    void set_ada_backward(int32_t value) { ada_backward_ = value; }

    // allele depths for others
    int32_t ado_forward() const { return ado_forward_; }
    void set_ado_forward(int32_t value) { ado_forward_ = value; }
    int32_t ado_backward() const { return ado_backward_; }
    void set_ado_backward(int32_t value) { ado_backward_ = value; }

    // qual-weighted allele depths for reference
    float wadr_forward() const { return wadr_forward_; }
    void set_wadr_forward(float value) { wadr_forward_ = value; }
    float wadr_backward() const { return wadr_backward_; }
    void set_wadr_backward(float value) { wadr_backward_ = value; }

    // qual-weighted allele depths for alt
    float wada_forward() const { return wada_forward_; }
    void set_wada_forward(float value) { wada_forward_ = value; }
    float wada_backward() const { return wada_backward_; }
    void set_wada_backward(float value) { wada_backward_ = value; }

    // qual-weighted allele depths for others
    float wado_forward() const { return wado_forward_; }
    void set_wado_forward(float value) { wado_forward_ = value; }
    float wado_backward() const { return wado_backward_; }
    void set_wado_backward(float value) { wado_backward_ = value; }

    Json::Value toJson() const
    {
        Json::Value val;
        if (!contig_.empty())
            val["contig"] = contig_;
        if (start_)
            val["start"] = start_;
        if (end_)
            val["end"] = end_;

        // STR shifting: left / rightmost positions
        if (leftmost_)
            val["leftmost"] = leftmost_;
        if (rightmost_)
            val["rightmost"] = rightmost_;

        if (!alt_.empty())
            val["alt"] = alt_;

        // allele depths for reference
        if (adr_forward_)
            val["adrForward"] = adr_forward_;
        if (adr_backward_)
            val["adrBackward"] = adr_backward_;

        // allele depths for alt
        if (ada_forward_)
            val["adaForward"] = ada_forward_;
        if (ada_backward_)
            val["adaBackward"] = ada_backward_;

        // allele depths for others
        if (ado_forward_)
            val["adoForward"] = ado_forward_;
        if (ado_backward_)
            val["adoBackward"] = ado_backward_;

        // qual-weighted allele depths for reference
        if (wadr_forward_)
            val["wadrForward"] = wadr_forward_;
        if (wadr_backward_)
            val["wadrBackward"] = wadr_backward_;

        // qual-weighted allele depths for alt
        if (wada_forward_)
            val["wadaForward"] = wada_forward_;
        if (wada_backward_)
            val["wadaBackward"] = wada_backward_;

        // qual-weighted allele depths for others
        if (wado_forward_)
            val["wadoForward"] = wado_forward_;
        if (wado_backward_)
            val["wadoBackward"] = wado_backward_;
        return val;
    }

private:
    // Each node is tagged with a set of sequences
    std::string contig_;
    int32_t start_ = 0;
    int32_t end_ = 0;

    // STR shifting: left / rightmost positions
    int32_t leftmost_ = 0;
    int32_t rightmost_ = 0;

    std::string alt_;

    // allele depths for reference
    int32_t adr_forward_ = 0;
    int32_t adr_backward_ = 0;

    // allele depths for alt
    int32_t ada_forward_ = 0;
    int32_t ada_backward_ = 0;

    // allele depths for others
    int32_t ado_forward_ = 0;
    int32_t ado_backward_ = 0;

    // qual-weighted allele depths for reference
    float wadr_forward_ = 0.f;
    float wadr_backward_ = 0.f;

    // qual-weighted allele depths for alt
    float wada_forward_ = 0.f;
    float wada_backward_ = 0.f;

    // qual-weighted allele depths for others
    float wado_forward_ = 0.f;
    float wado_backward_ = 0.f;
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
    void appendCoverage(
        graphtools::GraphCoordinates const& coords, std::string const& node_name, Json::Value& coverage) const;

private:
    struct VariantCandidateListImpl;
    std::unique_ptr<VariantCandidateListImpl> _impl;
};
};
