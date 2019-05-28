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
 * \Summary of graph alignment statistics for validation & filtering purpose
 *
 * \file GraphSummaryStatistics.cpp
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#include "paragraph/GraphSummaryStatistics.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "paragraph/ReadFilter.hh"

namespace paragraph
{

using graphtools::Alignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::map;
using std::string;
using std::vector;

/**
 * add summary statistics on graph alignments into JSON output
 */
void summarizeAlignments(Graph const& graph, common::ReadBuffer const& reads, Json::Value& output)
{
    std::vector<std::string> gkeys = { "nodes", "edges", "alleles" }; // keys for output items
    std::map<std::string, std::map<std::string, AlignmentStatistics>> gstats; // type -> node/edge/allele -> stats

    std::map<std::string, int> allele_score_sum; // allele -> sum of graph mapping scores
    std::map<std::string, int> broken_path; // allele -> #reads support broken path
    map<std::string, size_t> allele_lengths; // sequence id -> length

    for (NodeId n_id = 0; n_id < (NodeId)graph.numNodes(); ++n_id)
    {
        std::set<std::string> node_predecessor_sequence_ids;
        for (const auto& pred : graph.predecessors(n_id))
        {
            if (graph.hasEdge(pred, n_id))
            {
                for (const auto& label : graph.edgeLabels(pred, n_id))
                {
                    node_predecessor_sequence_ids.insert(label);
                }
            }
        }
        std::set<std::string> node_successor_sequence_ids;
        for (const auto& succ : graph.successors(n_id))
        {
            if (graph.hasEdge(n_id, succ))
            {
                for (const auto& label : graph.edgeLabels(n_id, succ))
                {
                    node_successor_sequence_ids.insert(label);
                }
            }
        }
        std::list<std::string> node_sequence_ids;
        std::set_intersection(
            node_predecessor_sequence_ids.begin(), node_predecessor_sequence_ids.end(),
            node_successor_sequence_ids.begin(), node_successor_sequence_ids.end(),
            std::back_inserter(node_sequence_ids));

        for (auto const& s : node_sequence_ids)
        {
            if (allele_lengths.find(s) == allele_lengths.end())
            {
                allele_lengths[s] = 0;
            }
            allele_lengths[s] += graph.nodeSeq(n_id).size();
        }
    }

    const bool has_source_or_sink
        = (graph.nodeName(0) == "source") || (graph.nodeName(static_cast<NodeId>(graph.numNodes() - 1)) == "sink");
    for (auto& read : reads)
    {
        if (read->graph_mapping_status() != common::Read::MAPPED)
        {
            continue;
        }
        GraphAlignment graph_alignment
            = graphtools::decodeGraphAlignment(read->graph_pos(), read->graph_cigar(), &graph);

        NodeId pred_node_id;
        for (size_t alignment_index = 0; alignment_index != graph_alignment.size(); ++alignment_index)
        {
            NodeId current_node_id = graph_alignment.getNodeIdByIndex(alignment_index);
            const bool is_source_or_sink
                = has_source_or_sink && (current_node_id == 0 || current_node_id == graph.numNodes() - 1);

            // node statistics
            const auto& node_name = graph.nodeName(current_node_id);
            if (gstats["nodes"].find(node_name) == gstats["nodes"].end())
            {
                gstats["nodes"][node_name] = AlignmentStatistics(graph.nodeSeq(current_node_id).size());
            }
            gstats["nodes"][node_name].addNodeMapping(
                graph_alignment[alignment_index], read->is_graph_reverse_strand(), !is_source_or_sink);

            // edge statistics
            if (alignment_index > 0)
            {
                const string edge_name = graph.nodeName(pred_node_id) + "_" + node_name;
                if (gstats["edges"].find(edge_name) == gstats["edges"].end())
                {
                    const size_t edge_length
                        = graph.nodeSeq(pred_node_id).size() + graph.nodeSeq(current_node_id).size();
                    gstats["edges"][edge_name] = AlignmentStatistics(edge_length);
                }
                gstats["edges"][edge_name].addEdgeMapping(
                    graph_alignment[alignment_index - 1], graph_alignment[alignment_index],
                    read->is_graph_reverse_strand(), has_source_or_sink && (current_node_id - 1 == 0),
                    is_source_or_sink);
            }
            pred_node_id = current_node_id;
        }

        for (auto& allele : read->graph_sequences_supported())
        {
            if (gstats["alleles"].find(allele) == gstats["alleles"].end())
            {
                gstats["alleles"][allele] = AlignmentStatistics(allele_lengths[allele]);
            }
            if (allele_score_sum.find(allele) == allele_score_sum.end())
            {
                allele_score_sum[allele] = 0;
            }
            gstats["alleles"][allele].addAlleleMapping(
                graph_alignment, read->is_graph_reverse_strand(), has_source_or_sink);
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
