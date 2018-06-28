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
 * \brief Graph read disambiguation code
 *
 * \file Disambiguation.cpp
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <fstream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string/join.hpp>

#include "common/Fragment.hh"
#include "common/Phred.hh"
#include "common/ReadExtraction.hh"
#include "common/ReadPairs.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphCoordinates.hh"
#include "grm/Align.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/GraphSummaryStatistics.hh"
#include "paragraph/GraphVariants.hh"
#include "paragraph/HaplotypePaths.hh"
#include "paragraph/ReadCounting.hh"
#include "paragraph/ReadFilter.hh"
#include "variant/Variant.hh"

#include "spdlog/spdlog.h"

// Error.hh always needs to go last
#include "common/Error.hh"

using common::Read;
using common::p_Read;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::GraphCoordinates;
using graphtools::NodeId;
using graphtools::decodeGraphAlignment;
using std::vector;

//#define DEBUG_DISAMBIGUATION

namespace paragraph
{

/**
 * Update sequence labels in read according to nodes the read has traversed
 * @param g graph structure
 * @param reads list of aligned reads
 * @param nodefilter filter to check if a read supports a particular node
 * @param edgefilter filter to check if a read supports a particular edge
 * @param paths
 */
void disambiguateReads(
    Graph* g, std::vector<common::p_Read>& reads, ReadSupportsNode nodefilter, ReadSupportsEdge edgefilter)
{
    for (auto& read : reads)
    {
        read->clear_graph_sequences_supported();
        read->clear_graph_nodes_supported();
        read->clear_graph_edges_supported();
        if (read->graph_mapping_status() == common::Read::MAPPED)
        {
            bool has_previous = false;
            NodeId pnode = 0;

            std::set<std::pair<std::string, std::string>> edges_supported_by_read;
            std::set<NodeId> nodes_supported_by_read;
            std::set<std::string> overlapped_pfams;

            GraphAlignment gm = decodeGraphAlignment(read->graph_pos(), read->graph_cigar(), g);
            auto const& path = gm.path();
            for (auto node = path.begin(); node != path.end(); ++node)
            {
                if (has_previous
                    && (edgefilter == nullptr || edgefilter(*read, g->nodeName(pnode), g->nodeName(*node))))
                {
                    edges_supported_by_read.emplace(g->nodeName(pnode), g->nodeName(*node));
                    for (const auto& s : g->edgeLabels(pnode, *node))
                    {
                        overlapped_pfams.insert(s);
                    }
                }
                has_previous = true;
                pnode = *node;

                // check if node is rejected
                if (nodefilter == nullptr || nodefilter(*read, g->nodeName(*node)))
                {
                    nodes_supported_by_read.emplace(*node);
                }
            }

            for (auto n : nodes_supported_by_read)
            {
                read->add_graph_nodes_supported(g->nodeName(n));
            }

            for (auto const& e : edges_supported_by_read)
            {
                read->add_graph_edges_supported(e.first + "_" + e.second);
            }

            for (auto const& label : overlapped_pfams)
            {
                graphtools::PathFamily pfam(g, label);
                if (pfam.containsPath(path))
                {
                    read->add_graph_sequences_supported(label);
                }
            }
        }
    }
}

/**
 * Align reads from single BAM file to graph and disambiguate reads
 * to produce counts.
 *
 * @param parameters alignment parameters
 * @param output_reads pass a pointer to a vector to retrieve all reads
 * @return results as JSON value
 */
Json::Value alignAndDisambiguate(const Parameters& parameters, common::ReadBuffer& all_reads)
{
    auto logger = LOG();

    // Initialize the graph aligner.
    graphtools::Graph graph = grm::graphFromJson(parameters.description(), parameters.reference_path());

    std::unordered_map<std::string, NodeId> node_id_map;
    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        node_id_map[graph.nodeName(node_id)] = node_id;
    }

    Json::Value output = parameters.description();
    output["reference"] = parameters.reference_path();

    common::ReadBuffer output_reads;

    if (parameters.output_enabled(Parameters::ALIGNMENTS) || parameters.output_enabled(Parameters::FILTERED_ALIGNMENTS))
    {
        output["alignments"] = Json::Value();
    }

    auto read_filter = createReadFilter(
        &graph, parameters.remove_nonuniq_reads(), parameters.bad_align_frac(), parameters.kmer_len());
    // remember total number of reads for later
    const size_t total_reads_input = all_reads.size();
    std::map<std::string, size_t> read_filter_counts;
    std::mutex output_mutex;
    auto read_filter_function
        = [&read_filter, &read_filter_counts, &parameters, &output, &output_reads, &output_mutex](Read& r) -> bool {
        const auto result_and_error = read_filter->filterRead(r);
        if (result_and_error.first && parameters.output_enabled(Parameters::FILTERED_ALIGNMENTS))
        {
            r.set_graph_mapping_status(common::Read::BAD_ALIGN);
            Json::Value r_json = r.toJson();
            r_json["error"] = result_and_error.second;
            {
                std::lock_guard<std::mutex> output_guard(output_mutex);
                auto count_it = read_filter_counts.find(result_and_error.second);
                if (count_it == read_filter_counts.end())
                {
                    read_filter_counts[result_and_error.second] = 1;
                }
                else
                {
                    count_it->second++;
                }
                output["alignments"].append(r_json);
                output_reads.emplace_back(new Read(r));
            }
        }
        return result_and_error.first;
    };

    grm::alignReads(
        &graph, grm::pathsFromJson(&graph, parameters.description()["paths"]), all_reads, read_filter_function,
        parameters.path_sequence_matching(), parameters.graph_sequence_matching(), parameters.klib_sequence_matching(),
        parameters.kmer_sequence_matching(), parameters.validate_alignments(), parameters.threads());

    auto nodefilter = [&graph, &node_id_map](Read& read, const std::string& node) -> bool {
        try
        {
            GraphAlignment alignment = decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), &graph);

            const auto node_id = node_id_map[node];
            const bool is_short_node = graph.nodeSeq(node_id).size() < read.bases().size() / 2;

            int32_t index = 0;
            for (const auto& node_alignment : alignment)
            {
                if (node_id == alignment.getNodeIdByIndex(index))
                {
                    const size_t nonmatch = node_alignment.numMismatched() + node_alignment.numClipped();
                    const size_t indel = node_alignment.numInserted() + node_alignment.numDeleted();

                    if (is_short_node && (nonmatch > 0 || indel > 0))
                    {
                        return false;
                    }
                    return nonmatch + indel <= read.bases().size() / 2;
                }
                ++index;
            }
        }
        catch (std::exception const&)
        {
            LOG()->warn("Invalid read mapping for {} : {}", read.fragment_id(), read.graph_cigar());
        }
        return false; // node not covered by read
    };

    auto edgefilter = [&graph, &node_id_map](Read& read, const std::string& node1, const std::string& node2) -> bool {
        try
        {
            GraphAlignment alignment = decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), &graph);

            const auto node_id1 = node_id_map[node1];
            const auto node_id2 = node_id_map[node2];

            const graphtools::Alignment* previous_alignment = nullptr;
            auto previous_node_id = static_cast<NodeId>(-1); // Large positive number
            int32_t index = 0;
            for (const auto& node_alignment : alignment)
            {
                const auto node_alignment_node_id = alignment.getNodeIdByIndex(index);
                if (previous_alignment != nullptr && previous_node_id == node_id1 && node_alignment_node_id == node_id2)
                {
                    auto const min_node_overlap = static_cast<int32_t>(read.bases().length() / 10 + 1);
                    bool status
                        = (previous_alignment->numMatched()
                               >= (unsigned)std::min(previous_alignment->referenceLength(), (uint32_t)min_node_overlap)
                           && node_alignment.numMatched()
                               >= (unsigned)std::min(node_alignment.referenceLength(), (uint32_t)min_node_overlap));

#ifdef DISABLE_ADDITIONAL_EDGE_FILTER
                    return status;
#else
                    if (status) // require softclip shorter than half of query length. Filter for cigars like 2[1M15S]
                    {
                        status = (previous_alignment->queryLength() < previous_alignment->referenceLength() * 2)
                            && (node_alignment.queryLength() < node_alignment.referenceLength() * 2);
                    }
                    if (status) // require minimum overlap for one node. Filter for cigars like 2[1M]
                    {
                        const auto node1_length = static_cast<int32_t>(graph.nodeSeq(node_id1).size());
                        const auto node2_length = static_cast<int32_t>(graph.nodeSeq(node_id2).size());
                        status = ((int32_t)previous_alignment->numMatched() >= std::min(node1_length, min_node_overlap))
                            && ((int32_t)node_alignment.numMatched() >= std::min(node2_length, min_node_overlap));
                    }
                    return status;
#endif
                }

                previous_alignment = &node_alignment;
                previous_node_id = node_alignment_node_id;
                ++index;
            }
        }
        catch (std::exception const&)
        {
            LOG()->warn("Invalid read mapping for {} : {}", read.fragment_id(), read.graph_cigar());
        }
        return false; // edge not covered by read
    };

    // add supported read haplotypes
    Json::Value paths = parameters.description()["paths"];
    if (parameters.output_enabled(Parameters::HAPLOTYPES))
    {
        addHaplotypePaths(all_reads, graph, paths, output);

        // update all edge labels -- we do this here because addHaplotypePaths
        // doesn't need to know about nodeIdMap
        for (auto& json_edge : output["edges"])
        {
            const NodeId from = node_id_map[json_edge["from"].asString()];
            const NodeId to = node_id_map[json_edge["to"].asString()];
            json_edge["sequences"] = Json::arrayValue;
            for (const auto& sequence : graph.edgeLabels(from, to))
            {
                json_edge["sequences"].append(sequence);
            }
        }
    }

    disambiguateReads(&graph, all_reads, nodefilter, edgefilter);

    graphtools::GraphCoordinates coordinates(&graph);
    countReads(
        coordinates, all_reads, output, parameters.output_enabled(Parameters::NODE_READ_COUNTS),
        parameters.output_enabled(Parameters::EDGE_READ_COUNTS),
        parameters.output_enabled(Parameters::PATH_READ_COUNTS),
        parameters.output_enabled(Parameters::DETAILED_READ_COUNTS));

    getVariants(
        coordinates, all_reads, output, parameters.min_reads_for_variant(), parameters.min_frac_for_variant(), paths,
        parameters.output_enabled(Parameters::VARIANTS), parameters.output_enabled(Parameters::NODE_COVERAGE),
        parameters.output_enabled(Parameters::PATH_COVERAGE));

    summarizeAlignments(graph, all_reads, output);
    double bad_alignment_pct = 0;
    if (total_reads_input > 0)
    {
        auto bad_align_count_it = read_filter_counts.find("bad_align");
        if (bad_align_count_it != read_filter_counts.end())
        {
            bad_alignment_pct = ((double)bad_align_count_it->second) / total_reads_input;
        }
    }
    output["alignment_statistics"]["bad_alignment_pct"] = bad_alignment_pct;
    for (auto const& read_filter_type : read_filter_counts)
    {
        output["alignment_statistics"]["read_filter_" + read_filter_type.first] = (Json::UInt64)read_filter_type.second;
    }

    if (parameters.output_enabled(Parameters::ALIGNMENTS))
    {
        output_reads.reserve(all_reads.size() + output_reads.size());
        for (auto& r : all_reads)
        {
            Json::Value r_json = r->toJson();
            output["alignments"].append(r_json);
            output_reads.emplace_back(std::move(r));
        }
    }
    all_reads = std::move(output_reads);

    return output;
}
}
