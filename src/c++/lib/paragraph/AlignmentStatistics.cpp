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
