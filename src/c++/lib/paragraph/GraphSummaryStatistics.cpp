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
 * \Summary of graph alignment statistics for validation & filtering purpose
 *
 * \file GraphSummaryStatistics.cpp
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#include "paragraph/GraphSummaryStatistics.hh"
#include "graphs/GraphMappingOperations.hh"
#include "paragraph/ReadFilter.hh"

namespace paragraph
{

using graphs::GraphMapping;
using graphs::NodeMapping;
using std::map;
using std::string;
using std::vector;

/**
 * add summary statistics on graph alignments into JSON output
 */
void summarizeAlignments(graphs::WalkableGraph& wgraph, vector<common::p_Read>& reads, Json::Value& output)
{
    std::vector<std::string> gkeys = { "nodes", "edges", "alleles" }; // keys for output items
    std::map<std::string, std::map<std::string, AlignmentStatistics>> gstats; // type -> node/edge/allele -> stats

    std::map<std::string, int> allele_score_sum; // allele -> sum of graph mapping scores
    std::map<std::string, int> broken_path; // allele -> #reads support broken path
    map<uint64_t, size_t> allele_lengths; // sequence id -> length

    for (auto& n_id : wgraph.allNodes())
    {
        for (auto s : wgraph.node(n_id)->sequence_ids())
        {
            if (allele_lengths.find(s) == allele_lengths.end())
            {
                allele_lengths[s] = 0;
            }
            allele_lengths[s] += wgraph.node(n_id)->sequence().size();
        }
    }

    for (auto& read : reads)
    {
        if (read->graph_mapping_status() != reads::MappingStatus::MAPPED)
        {
            continue;
        }
        GraphMapping graph_mapping
            = graphs::decodeFromString(read->graph_pos(), read->graph_cigar(), read->bases(), wgraph);

        const NodeMapping* ptr_prev_node = nullptr;
        string previous_node_name;

        for (const auto& node_mapping : graph_mapping)
        {
            // node stat
            string node_name = wgraph.nodeName(static_cast<uint64_t>(node_mapping.node_id));
            if (gstats["nodes"].find(node_name) == gstats["nodes"].end())
            {
                gstats["nodes"][node_name]
                    = AlignmentStatistics(wgraph.node((uint64_t)node_mapping.node_id)->sequence().size());
            }
            gstats["nodes"][node_name].addNodeMapping(node_mapping, read->is_graph_reverse_strand());
            if ((uint64_t)node_mapping.node_id != wgraph.source() && (uint64_t)node_mapping.node_id != wgraph.sink())
            {
                gstats["nodes"][node_name].addClipBases(node_mapping);
            }

            // edge stats
            if (ptr_prev_node != NULL)
            {
                string edge_name = previous_node_name + "_" + node_name;
                if (gstats["edges"].find(edge_name) == gstats["edges"].end())
                {
                    size_t edge_length = wgraph.node(static_cast<uint64_t>(ptr_prev_node->node_id))->sequence().size()
                        + wgraph.node(static_cast<uint64_t>(node_mapping.node_id))->sequence().size();
                    gstats["edges"][edge_name] = AlignmentStatistics(edge_length);
                }
                gstats["edges"][edge_name].addEdgeMapping(
                    *ptr_prev_node, node_mapping, read->is_graph_reverse_strand());
                if ((uint64_t)node_mapping.node_id != wgraph.source()
                    && (uint64_t)node_mapping.node_id != wgraph.sink())
                {
                    gstats["edges"][edge_name].addClipBases(node_mapping);
                }
                if ((uint64_t)ptr_prev_node->node_id != wgraph.source()
                    && (uint64_t)ptr_prev_node->node_id != wgraph.sink())
                {
                    gstats["edges"][edge_name].addClipBases(*ptr_prev_node);
                }
            }
            else
            {
                ptr_prev_node = &node_mapping;
            }
            previous_node_name = node_name;
        }

        for (auto& allele : read->graph_sequences_supported())
        {
            if (gstats["alleles"].find(allele) == gstats["alleles"].end())
            {
                gstats["alleles"][allele] = AlignmentStatistics(allele_lengths[wgraph.sequenceId(allele)]);
            }
            if (allele_score_sum.find(allele) == allele_score_sum.end())
            {
                allele_score_sum[allele] = 0;
            }
            gstats["alleles"][allele].addAlleleMapping(
                graph_mapping, read->is_graph_reverse_strand(), wgraph.source(), wgraph.sink());
            allele_score_sum[allele] += read->graph_alignment_score();
        }
        for (auto& allele : read->graph_sequences_broken())
        {
            if (broken_path.find(allele) == broken_path.end())
            {
                broken_path[allele] = 0;
            }
            broken_path[allele]++;
        }
    }

    output["alignment_statistics"] = Json::Value(Json::ValueType::objectValue);
    for (auto& gkey : gkeys)
    {
        output["alignment_statistics"][gkey] = Json::Value(Json::ValueType::objectValue);
        auto& p_json = output["alignment_statistics"][gkey];
        for (auto& k2 : gstats[gkey])
        {
            p_json[k2.first] = k2.second.toJson();
            if (gkey == "alleles")
            {
                p_json[k2.first]["avr_score"]
                    = k2.second.num_reads() == 0 ? 0 : (double)allele_score_sum[k2.first] / k2.second.num_reads();
                if (broken_path.find(k2.first) != broken_path.end())
                {
                    p_json[k2.first]["num_reads_for_broken_path"] = broken_path[k2.first];
                }
            }
        }
    }
}
}
