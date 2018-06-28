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
#include <vector>

#include "spdlog/spdlog.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string/join.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphCoordinates.hh"
#include "paragraph/GraphSummaryStatistics.hh"
#include "paragraph/GraphVariants.hh"
#include "paragraph/HaplotypePaths.hh"
#include "paragraph/ReadFilter.hh"
#include "variant/Variant.hh"

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
 * Take apart CIGAR string and collect variant candidates
 * @param read read after alignment
 * @param target vector of candidate lists
 */
void updateVariantCandidateLists(graphtools::Graph const* g, common::Read const& read, NodeCandidates& target)
{
    const std::string& graph_cigar = read.graph_cigar();
    int pos_in_node = read.graph_pos();
    size_t pos = 0;
    std::string remaining_read = read.bases();
    while (pos < graph_cigar.size())
    {
        std::string nodenum;
        while (graph_cigar[pos] >= '0' && graph_cigar[pos] <= '9')
        {
            nodenum += graph_cigar[pos];
            ++pos;
        }
        assert(!nodenum.empty());

        std::string nodecigar;
        assert(graph_cigar[pos] == '[');
        pos++;
        while (graph_cigar[pos] != ']')
        {
            nodecigar += graph_cigar[pos];
            ++pos;
        }
        assert(graph_cigar[pos] == ']');
        pos++;

        const auto i_nodenum = static_cast<const NodeId>(atoll(nodenum.c_str()));
        int ref_left = 0;
        int alt_left = 0;
        std::list<variant::RefVar> vars_this_node = variant::cigarToRefVar(
            g->nodeSeq(i_nodenum).substr((unsigned long)pos_in_node), remaining_read, nodecigar, ref_left, alt_left,
            true);

#ifdef DEBUG_DISAMBIGUATION
        std::cerr << "Read " << read.fragment_id() << " adds the following variants to node " << i_nodenum;
        for (const auto& var : vars_this_node)
        {
            std::cerr << " " << var;
        }
        std::cerr << std::endl;
#endif

        remaining_read = remaining_read.substr(remaining_read.size() - alt_left);

        auto vcl_it = target.find(i_nodenum);
        if (vcl_it == target.end())
        {
            vcl_it = target.emplace(i_nodenum, variant::VariantCandidateList(g->nodeSeq(i_nodenum))).first;
        }

        int64_t last_end = -1;
        for (auto& var : vars_this_node)
        {
            var.start += pos_in_node;
            var.end += pos_in_node;

            int mean_qual = 0;
            // var.flags gives pos in read
            if (var.flags >= 0 && var.flags < (signed)read.bases().size())
            {
                std::string qual_substr;
                if (!var.alt.empty()) // insertion or substitution: use mean qual across bases
                {
                    qual_substr = read.quals().substr((unsigned long)var.flags, var.alt.size());
                }
                else // deletion: use bases before and after
                {
                    const int64_t vstart = std::max((int64_t)0, var.flags - 1);
                    const int64_t vend = std::max((int64_t)0, var.flags);
                    qual_substr = read.quals().substr((unsigned long)vstart, (unsigned long)(vend - vstart + 1));
                }

                double fqual = 0.0f;
                for (auto x : qual_substr)
                {
                    fqual += (common::phred::phredToErrorProb(x - 33));
                }
                if (qual_substr.size() > 1)
                {
                    fqual /= qual_substr.size();
                }
                mean_qual = (int)common::phred::errorProbToPhred(fqual);
            }

            last_end = std::max(
                last_end,
                vcl_it->second.addRefVarObservation(var, read.is_graph_reverse_strand(), last_end, mean_qual));
        }

        pos_in_node = 0;
    }
}

/**
 * Extract on-graph variants
 * @param coordinates graph coordinates and graph information
 * @param reads list of reads
 * @param output  output JSON node
 * @param min_reads_for_variant minumum number of reads that must support a variant
 * @param min_frac_for_variant minimum fraction of reads that must support a variant
 * @param paths set of paths to compute coverage over
 * @param write_variants output variants
 * @param write_node_coverage output coverage for nodes
 * @param write_node_coverage output coverage for paths
 */
void getVariants(
    graphtools::GraphCoordinates const& coordinates, common::ReadBuffer const& reads, Json::Value& output,
    int min_reads_for_variant, float min_frac_for_variant, Json::Value const& paths, bool write_variants,
    bool write_node_coverage, bool write_path_coverage)
{
    graphtools::Graph const& graph(coordinates.getGraph());
    std::unordered_map<std::string, NodeId> node_id_map;
    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        node_id_map[graph.nodeName(node_id)] = node_id;
    }
    // collect variant candidates for every node
    NodeCandidates candidates;
    std::unordered_map<std::string, NodeCandidates> candidates_by_sequence;
    for (const auto& r : reads)
    {
        try
        {
            if (write_variants || write_node_coverage)
            {
                updateVariantCandidateLists(&graph, *r, candidates);
            }
            if (write_path_coverage)
            {
                for (const auto& seq : r->graph_sequences_supported())
                {
                    auto candidate_list = candidates_by_sequence.find(seq);
                    if (candidate_list == candidates_by_sequence.end())
                    {
                        candidate_list = candidates_by_sequence.emplace(seq, NodeCandidates()).first;
                    }
                    updateVariantCandidateLists(&graph, *r, candidate_list->second);
                }
            }
        }
        catch (std::exception const& e)
        {
            LOG()->warn(
                "Read {} cigar {} could not be used to produce candidate lists: {}", r->fragment_id(), r->graph_cigar(),
                e.what());
        }
    }

    if (write_variants)
    {
        output["variants"] = Json::Value(Json::ValueType::objectValue);
        for (auto const& node_candidates : candidates)
        {
            const std::string& node_name = graph.nodeName(node_candidates.first);
            output["variants"][node_name] = Json::Value(Json::ValueType::arrayValue);
            for (auto variant : node_candidates.second.getVariants())
            {
                const int variant_alt_count = variant->ada_backward() + variant->ada_forward();
                const int variant_total_count = variant->adr_backward() + variant->adr_forward()
                    + variant->ada_backward() + variant->ada_forward() + variant->ado_backward()
                    + variant->ado_forward();

                if (variant_alt_count < min_reads_for_variant
                    || float(variant_alt_count) / variant_total_count < min_frac_for_variant)
                {
                    continue;
                }

                Json::Value val = variant->toJson();
                output["variants"][node_name].append(val);
            }
        }
    }
    if (write_node_coverage)
    {
        output["node_coverage"] = Json::Value(Json::ValueType::objectValue);
        for (auto const& node_candidates : candidates)
        {
            const std::string node_name = graph.nodeName(node_candidates.first);
            output["node_coverage"][node_name] = Json::Value(Json::ValueType::objectValue);
            node_candidates.second.appendCoverage(coordinates, node_name, output["node_coverage"][node_name]);
        }
    }
    if (write_path_coverage)
    {
        output["path_coverage"] = Json::Value(Json::ValueType::objectValue);
        for (auto const& p : paths)
        {
            const auto& path_id = p["path_id"].asString();
            const auto& sequence_id = p["sequence"].asString();

            const auto this_sequence_candidates = candidates_by_sequence.find(sequence_id);

            output["path_coverage"][path_id] = Json::Value(Json::ValueType::objectValue);
            for (auto const& node_name : p["nodes"])
            {
                auto const& node_name_str = node_name.asString();
                assert(node_id_map.count(node_name_str));
                auto const node_id = node_id_map[node_name_str];
                if (this_sequence_candidates == candidates_by_sequence.end())
                {
                    variant::VariantCandidateList vcl(graph.nodeSeq(node_id));
                    vcl.appendCoverage(coordinates, node_name_str, output["path_coverage"][path_id]);
                }
                else
                {
                    auto nc = this_sequence_candidates->second.find(node_id);
                    if (nc == this_sequence_candidates->second.end())
                    {
                        variant::VariantCandidateList vcl(graph.nodeSeq(node_id));
                        vcl.appendCoverage(coordinates, node_name_str, output["path_coverage"][path_id]);
                    }
                    else
                    {
                        nc->second.appendCoverage(coordinates, node_name_str, output["path_coverage"][path_id]);
                    }
                }
            }
        }
    }
}
}
