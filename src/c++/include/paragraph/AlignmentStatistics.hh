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
