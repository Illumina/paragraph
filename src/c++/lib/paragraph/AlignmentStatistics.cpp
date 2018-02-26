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

#include "paragraph/AlignmentStatistics.hh"
#include "graphs/GraphMappingOperations.hh"

namespace paragraph
{

using graphs::GraphMapping;
using graphs::NodeMapping;
using std::map;
using std::string;
using std::vector;

void AlignmentStatistics::addNodeMappingBases(const NodeMapping& node_mapping)
{
    num_match_bases += node_mapping.mapping.matched();
    num_mismatch_bases += node_mapping.mapping.mismatched();
    num_gap_bases += node_mapping.mapping.inserted() + node_mapping.mapping.deleted();
}

void AlignmentStatistics::addNodeMapping(const NodeMapping& node_mapping, bool is_graph_reverse_strand)
{
    addNodeMappingBases(node_mapping);
    is_graph_reverse_strand ? num_rev_strand_reads++ : num_fwd_strand_reads++;
}

void AlignmentStatistics::addEdgeMapping(
    const NodeMapping& node_mapping1, const NodeMapping& node_mapping2, bool is_graph_reverse_strand)
{
    addNodeMappingBases(node_mapping1);
    addNodeMapping(node_mapping2, is_graph_reverse_strand);
}

void AlignmentStatistics::addAlleleMapping(
    const GraphMapping& graph_mapping, bool is_graph_reverse_strand, uint64_t source, uint64_t sink)
{
    bool is_first_node = true;
    for (auto& node_mapping : graph_mapping)
    {
        if (is_first_node)
        {
            addNodeMapping(node_mapping, is_graph_reverse_strand);
            is_first_node = false;
        }
        else
        {
            addNodeMappingBases(node_mapping);
        }
        if ((uint64_t)node_mapping.node_id != source && (uint64_t)node_mapping.node_id != sink)
        {
            addClipBases(node_mapping);
        }
    }
}

Json::Value AlignmentStatistics::toJson()
{
    Json::Value output;
    output["num_fwd_reads"] = num_fwd_strand_reads;
    output["num_rev_reads"] = num_rev_strand_reads;
    output["mismatch_rate"] = (double)num_mismatch_bases / (num_match_bases + num_mismatch_bases + num_gap_bases);
    output["gap_rate"] = (double)num_gap_bases / (num_match_bases + num_mismatch_bases + num_gap_bases);
    output["clip_rate"] = (double)num_clip_bases / (num_match_bases + num_mismatch_bases + num_gap_bases);
    output["match_base_depth"] = (double)num_match_bases / length;
    output["contig_length"] = (int)length;
    return output;
}
}
