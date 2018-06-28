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
 * \summary graph alignment statistics helper class
 *
 * \file AlignmentStatistics.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#pragma once

#include <cstddef>

#include "graphalign/GraphAlignment.hh"
#include "json/json.h"

namespace paragraph
{

/**
 * Alignment statistics for ONE graph element: can be node, edge or allele
 */
class AlignmentStatistics
{
public:
    AlignmentStatistics() = default;

    explicit AlignmentStatistics(const size_t _length)
        : length(_length)
        , num_match_bases(0)
        , num_mismatch_bases(0)
        , num_gap_bases(0)
        , num_clip_bases(0)
        , num_fwd_strand_reads(0)
        , num_rev_strand_reads(0)
    {
    }

    /**
     * simple getters
     */
    int num_reads() const { return num_fwd_strand_reads + num_rev_strand_reads; }

    /* Add mapping stats for a node. will count the read as forward / reverse and
     * add all types of bases as mapped to the node */
    void addNodeMapping(const graphtools::Alignment& alignment, bool is_graph_reverse_strand, bool count_clipped_bases);

    /**
     * Add mapping for an edge: adds read + counts for an edge
     */
    void addEdgeMapping(
        const graphtools::Alignment& from_alignment, const graphtools::Alignment& to_alignment,
        bool is_graph_reverse_strand, bool count_clipped_bases_from, bool count_clipped_bases_to);

    /**
     * add mapping stats of a read to an allele, don't count clipping for source / sink where appropriate
     */
    void addAlleleMapping(
        const graphtools::GraphAlignment& graph_alignment, bool is_graph_reverse_strand, bool has_source_and_sink);

    /**
     * output alignment statistics to a JSON value
     */
    Json::Value toJson();

private:
    /**
     * description for this element
     */
    size_t length;

    /**
     * internal data structure for base-level statistics
     */
    int num_match_bases;
    int num_mismatch_bases;
    int num_gap_bases;
    int num_clip_bases;

    /**
     * internal data structure for read-level statistics
     */
    int num_fwd_strand_reads;
    int num_rev_strand_reads;
};
}
