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

using graphtools::Alignment;
using graphtools::GraphAlignment;
using graphtools::NodeId;

namespace paragraph
{

using std::map;
using std::string;
using std::vector;

/* Add mapping stats for a node. will count the read as forward / reverse and
 * add all types of bases as mapped to the node */
void AlignmentStatistics::addNodeMapping(
    const Alignment& alignment, bool is_graph_reverse_strand, bool count_clipped_bases)
{
    num_match_bases += alignment.numMatched();
    num_mismatch_bases += alignment.numMismatched();
    num_gap_bases += alignment.numInserted() + alignment.numDeleted();

    if (count_clipped_bases)
    {
        num_clip_bases += alignment.numClipped();
    }

    if (is_graph_reverse_strand)
    {
        num_rev_strand_reads++;
    }
    else
    {
        num_fwd_strand_reads++;
    }
}

void AlignmentStatistics::addEdgeMapping(
    const Alignment& from_alignment, const Alignment& to_alignment, bool is_graph_reverse_strand,
    bool count_clipped_bases_from, bool count_clipped_bases_to)
{
    num_match_bases += from_alignment.numMatched();
    num_mismatch_bases += from_alignment.numMismatched();
    num_gap_bases += from_alignment.numInserted() + from_alignment.numDeleted();
    num_match_bases += to_alignment.numMatched();
    num_mismatch_bases += to_alignment.numMismatched();
    num_gap_bases += to_alignment.numInserted() + to_alignment.numDeleted();

    if (count_clipped_bases_from)
    {
        num_clip_bases += from_alignment.numClipped();
    }

    if (count_clipped_bases_to)
    {
        num_clip_bases += to_alignment.numClipped();
    }

    if (is_graph_reverse_strand)
    {
        num_rev_strand_reads++;
    }
    else
    {
        num_fwd_strand_reads++;
    }
}

/**
 * add mapping stats of a read to an allele: adds reads + counts if mapping contains a path from source to sink;
 * bases before source and sink are added as clipped.
 */
void AlignmentStatistics::addAlleleMapping(
    const GraphAlignment& graph_alignment, bool is_graph_reverse_strand, bool has_source_and_sink)
{
    const NodeId source = 0;
    const auto sink = static_cast<const NodeId>(graph_alignment.path().graphRawPtr()->numNodes() - 1);

    for (NodeId node_index = 0; node_index != graph_alignment.size(); ++node_index)
    {
        const auto node_id = static_cast<const NodeId>(graph_alignment.getNodeIdByIndex(node_index));
        const auto& alignment = graph_alignment[node_index];

        num_match_bases += alignment.numMatched();
        num_mismatch_bases += alignment.numMismatched();
        num_gap_bases += alignment.numInserted() + alignment.numDeleted();

        if (has_source_and_sink && ((node_id == source) || (node_id == sink)))
        {
            continue;
        }
        num_clip_bases += alignment.numClipped();
    }
    if (is_graph_reverse_strand)
    {
        num_rev_strand_reads++;
    }
    else
    {
        num_fwd_strand_reads++;
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
    if (length > 0)
    {
        output["match_base_depth"] = (double)num_match_bases / length;
    }
    output["contig_length"] = (int)length;
    return output;
}
}
