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
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string/join.hpp>

#include <google/protobuf/util/json_util.h>

#include "common/Fragment.hh"
#include "common/Phred.hh"
#include "common/ReadExtraction.hh"
#include "common/ReadPairs.hh"
#include "graphs/GraphMapping.hh"
#include "grm/Align.hh"
#include "paragraph/Disambiguation.hh"
#include "variant/Variant.hh"

#include "spdlog/spdlog.h"

// Error.hh always needs to go last
#include "common/Error.hh"

using common::Read;
using common::p_Read;
using graphs::GraphMapping;
using std::vector;

//#define DEBUG_DISAMBIGUATION

namespace paragraph
{

/**
 * Add fragment to counts
 * @param fragment the read
 * @param target target count JSON
 * @param nodes counts for all nodes
 * @param edges counts for all edges (node sequences of length 2)
 * @param paths counts for every node sequence
 */
void countFragment(
    common::Fragment const& fragment, Json::Value& target, bool nodes = false, bool edges = false, bool paths = false)
{
    std::list<std::string> to_append_unstranded = { "total" };

    if (nodes || edges || paths)
    {
        if (nodes)
        {
            for (const auto& n : fragment.graph_nodes_supported())
            {
                to_append_unstranded.push_back(n);
            }
        }
        if (edges)
        {
            for (const auto& e : fragment.graph_edges_supported())
            {
                to_append_unstranded.push_back(e);
            }
        }
        if (paths)
        {
            to_append_unstranded.push_back(boost::algorithm::join(fragment.graph_nodes_supported(), "_"));
        }
    }

    std::set<std::string> to_append;
    for (auto const& s : to_append_unstranded)
    {
        to_append.insert(s);
    }

    for (auto const& n : to_append)
    {
        if (!target.isMember(n))
        {
            target[n] = 1ull;
            target[n + ":READS"] = fragment.get_n_reads();
            target[n + ":FWD"] = fragment.get_n_graph_forward_reads();
            target[n + ":REV"] = fragment.get_n_graph_reverse_reads();
        }
        else
        {
            target[n] = target[n].asUInt64() + 1ull;
            target[n + ":READS"] = target[n + ":READS"].asUInt64() + fragment.get_n_reads();
            target[n + ":FWD"] = target[n + ":FWD"].asUInt64() + fragment.get_n_graph_forward_reads();
            target[n + ":REV"] = target[n + ":REV"].asUInt64() + fragment.get_n_graph_reverse_reads();
        }
    }
}

/**
 * Count disambiguated reads by fragment
 * @param coordinates graph and coordinate information
 * @param reads list of reads
 * @param output output JSON node
 * @param by_node output per-node counts
 * @param by_edge output per-edge counts
 * @param by_path output per-path counts
 * @param by_path_detailed output node and edge counts for all paths
 */
void countReads(
    graphs::GraphCoordinates const& coordinates, vector<common::p_Read>& reads, Json::Value& output,
    bool by_node = true, bool by_edge = true, bool by_path = true, bool by_path_detailed = false)
{
    if (by_node)
    {
        output["read_counts_by_node"] = Json::ValueType::objectValue;
    }

    if (by_edge)
    {
        output["read_counts_by_edge"] = Json::ValueType::objectValue;
    }

    if (by_path)
    {
        output["read_counts_by_sequence"] = Json::ValueType::objectValue;
    }

    common::FragmentList fragments;
    common::readsToFragments(coordinates, reads, fragments);

    {
        using namespace boost;
        using namespace boost::accumulators;

        typedef accumulator_set<double, features<tag::mean, tag::median, tag::variance, tag::density>> acc;
        typedef iterator_range<std::vector<std::pair<double, double>>::iterator> histogram_type;

        acc fragment_size(tag::density::num_bins = 20, tag::density::cache_size = 10);
        acc graph_fragment_size(tag::density::num_bins = 20, tag::density::cache_size = 10);

        uint64_t problematic_fragments_linear = 0;
        uint64_t problematic_fragments_graph = 0;
        uint64_t single_read_fragments = 0;
        uint64_t paired_read_fragments = 0;
        uint64_t multi_read_fragments = 0;
        for (auto& f : fragments)
        {
            const uint64_t fsize = f->get_bam_fragment_length();
            const uint64_t graph_fsize = f->get_graph_fragment_length();

            if (fsize != std::numeric_limits<uint64_t>::max())
            {
                if (f->get_n_reads() >= 2)
                {
                    fragment_size(fsize);
                }
            }
            else
            {
                ++problematic_fragments_linear;
            }
            if (graph_fsize != std::numeric_limits<uint64_t>::max())
            {
                if (f->get_n_reads() >= 2)
                {
                    graph_fragment_size(graph_fsize);
                }
            }
            else
            {
                ++problematic_fragments_graph;
            }

            if (f->get_n_reads() == 1)
            {
                ++single_read_fragments;
            }
            else if (f->get_n_reads() == 2)
            {
                ++paired_read_fragments;
            }
            else
            {
                ++multi_read_fragments;
            }
        }

        output["fragment_statistics"] = Json::ValueType::objectValue;

        output["fragment_statistics"]["mean_linear"] = mean(fragment_size);
        output["fragment_statistics"]["mean_graph"] = mean(graph_fragment_size);
        output["fragment_statistics"]["median_linear"] = median(fragment_size);
        output["fragment_statistics"]["median_graph"] = median(graph_fragment_size);
        output["fragment_statistics"]["variance_linear"] = variance(fragment_size);
        output["fragment_statistics"]["variance_graph"] = variance(graph_fragment_size);

        output["fragment_statistics"]["single_read"] = (Json::UInt64)single_read_fragments;
        output["fragment_statistics"]["paired_read"] = (Json::UInt64)paired_read_fragments;
        output["fragment_statistics"]["multi_read"] = (Json::UInt64)multi_read_fragments;
        output["fragment_statistics"]["problematic_linear"] = (Json::UInt64)problematic_fragments_linear;
        output["fragment_statistics"]["problematic_graph"] = (Json::UInt64)problematic_fragments_graph;

        output["fragment_statistics"]["linear_histogram"] = Json::arrayValue;
        output["fragment_statistics"]["graph_histogram"] = Json::arrayValue;

        histogram_type linear_hist = density(fragment_size);
        histogram_type graph_hist = density(graph_fragment_size);

        for (auto const& bin : linear_hist)
        {
            Json::Value bin_value = Json::objectValue;
            bin_value["lb"] = bin.first;
            bin_value["value"] = bin.second;
            output["fragment_statistics"]["linear_histogram"].append(bin_value);
        }
        for (auto const& bin : graph_hist)
        {
            Json::Value bin_value = Json::objectValue;
            bin_value["lb"] = bin.first;
            bin_value["value"] = bin.second;
            output["fragment_statistics"]["graph_histogram"].append(bin_value);
        }
    }

    for (auto& f : fragments)
    {
        if (by_node)
        {
            countFragment(*f, output["read_counts_by_node"], true, false, false);
        }
        if (by_edge)
        {
            countFragment(*f, output["read_counts_by_edge"], false, true, false);
        }
        if (by_path)
        {
            if (!f->graph_sequences_supported().empty())
            {
                std::vector<std::string> seqs;
                seqs.reserve(f->graph_sequences_supported().size());
                seqs.insert(seqs.end(), f->graph_sequences_supported().begin(), f->graph_sequences_supported().end());
                std::sort(seqs.begin(), seqs.end());
                auto joined_sequence_name = boost::algorithm::join(seqs, ",");
                if (!output["read_counts_by_sequence"].isMember(joined_sequence_name))
                {
                    output["read_counts_by_sequence"][joined_sequence_name] = Json::Value();
                }
                countFragment(
                    *f, output["read_counts_by_sequence"][joined_sequence_name], by_path_detailed, by_path_detailed,
                    by_path_detailed);
            }
        }
    }
}

/**
 * Take apart CIGAR string and collect variant candidates
 * @param read read after alignment
 * @param target vector of candidate lists
 */
void updateVariantCandidateLists(
    graphs::Graph& g, Read const& read, std::unordered_map<uint64_t, variant::VariantCandidateList>& target)
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

        const auto i_nodenum = static_cast<const uint64_t>(atoll(nodenum.c_str()));
        int ref_left = 0;
        int alt_left = 0;
        std::list<variant::RefVar> vars_this_node = variant::cigarToRefVar(
            g.nodes[i_nodenum]->sequence().substr((unsigned long)pos_in_node), remaining_read, nodecigar, ref_left,
            alt_left, true);

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
            vcl_it = target.emplace(i_nodenum, variant::VariantCandidateList(g.nodes[i_nodenum]->sequence())).first;
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
    graphs::GraphCoordinates const& coordinates, vector<common::p_Read>& reads, Json::Value& output,
    int min_reads_for_variant, float min_frac_for_variant, Json::Value const& paths, bool write_variants = false,
    bool write_node_coverage = false, bool write_path_coverage = false)
{
    graphs::WalkableGraph const& graph(coordinates.getGraph());
    // collect variant candidates for every node
    typedef std::unordered_map<uint64_t, variant::VariantCandidateList> NodeCandidates;
    NodeCandidates candidates;
    std::unordered_map<std::string, NodeCandidates> candidates_by_sequence;
    for (const auto& r : reads)
    {
        try
        {
            if (write_variants || write_node_coverage)
            {
                updateVariantCandidateLists(graph, *r, candidates);
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
                    updateVariantCandidateLists(graph, *r, candidate_list->second);
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
        Json::Reader jsonReader;
        for (auto const& node_candidates : candidates)
        {
            const std::string node_name = graph.node(node_candidates.first)->name();
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

                Json::Value val;
                std::string str;
                google::protobuf::util::MessageToJsonString(*((google::protobuf::Message*)variant), &str);
                jsonReader.parse(str, val);
                output["variants"][node_name].append(val);
            }
        }
    }
    if (write_node_coverage)
    {
        output["node_coverage"] = Json::Value(Json::ValueType::objectValue);
        for (auto const& node_candidates : candidates)
        {
            const std::string node_name = graph.node(node_candidates.first)->name();
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
            for (auto const& node_id : p["nodes"])
            {
                auto const& node_id_str = node_id.asString();
                if (this_sequence_candidates == candidates_by_sequence.end())
                {
                    variant::VariantCandidateList vcl(graph.node(node_id_str)->sequence());
                    vcl.appendCoverage(coordinates, node_id_str, output["path_coverage"][path_id]);
                }
                else
                {
                    auto nc = this_sequence_candidates->second.find(graph.nodeId(node_id_str));
                    if (nc == this_sequence_candidates->second.end())
                    {
                        variant::VariantCandidateList vcl(graph.node(node_id_str)->sequence());
                        vcl.appendCoverage(coordinates, node_id_str, output["path_coverage"][path_id]);
                    }
                    else
                    {
                        nc->second.appendCoverage(coordinates, node_id_str, output["path_coverage"][path_id]);
                    }
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
Json::Value alignAndDisambiguate(Parameters& parameters, common::ReadBuffer& all_reads)
{
    auto logger = LOG();

    // Initialize the graph aligner.
    graphs::Graph graph;

    graphs::fromJson(parameters.description(), parameters.reference_path(), graph);

    graphs::WalkableGraph wgraph(graph);
    graphs::GraphCoordinates coordinates(wgraph);

    Json::Value output = parameters.description();
    output["reference"] = parameters.reference_path();
    output["bam"] = parameters.bam_path();

    common::ReadBuffer output_reads;

    size_t bad_align = 0;
    size_t nonuniq_align = 0;

    if (parameters.output_enabled(Parameters::ALIGNMENTS) || parameters.output_enabled(Parameters::FILTERED_ALIGNMENTS))
    {
        output["alignments"] = Json::Value();
    }

    const float bad_align_frac = parameters.bad_align_frac();
    std::mutex output_mutex;
    auto read_filter = [&wgraph, bad_align_frac, &bad_align, &nonuniq_align, &parameters, &output, &output_reads,
                        &output_mutex](Read& r) -> bool {
        bool result = false;
        std::string error;
        if (!r.is_graph_alignment_unique())
        {
            ++nonuniq_align;
            error = "nonuniq";
            result = true;
        }
        graphs::GraphMapping mapping(r.graph_pos(), r.graph_cigar(), r.bases(), wgraph);
        if (mapping.querySpan() < 16 || r.graph_alignment_score() < (int)(bad_align_frac * (mapping.querySpan())))
        {
            bad_align++;
            error = "bad_align";
            result = true;
        }
        if (result && parameters.output_enabled(Parameters::FILTERED_ALIGNMENTS))
        {
            Json::Value r_json = r.asJson();
            r_json["error"] = error;
            std::lock_guard<std::mutex> output_guard(output_mutex);
            output["alignments"].append(r_json);
            output_reads.emplace_back(new Read(r));
        }
        return result;
    };

    grm::alignReads(
        wgraph, parameters.description()["paths"], all_reads, read_filter, parameters.exact_sequence_matching(),
        parameters.threads());

    auto nodefilter = [&wgraph](Read& read, const std::string& node) -> bool {
        try
        {
            graphs::GraphMapping mapping(read.graph_pos(), read.graph_cigar(), read.bases(), wgraph);

            const auto node_id = wgraph.nodeId(node);
            const auto nodeinfo = wgraph.node(node);
            const bool is_short_node = nodeinfo->sequence().size() < read.bases().size() / 2;

            for (size_t i = 0; i < mapping.size(); ++i)
            {
                auto const& nodemapping = mapping[i];
                const uint64_t this_node = nodemapping.node_id();

                if (node_id == this_node)
                {
                    const size_t nonmatch = nodemapping.mismatched() + nodemapping.clipped();
                    const size_t indel = nodemapping.inserted() + nodemapping.deleted();

                    if (is_short_node && (nonmatch > 0 || indel > 0))
                    {
                        return false;
                    }
                    return nonmatch + indel <= read.bases().size() / 2;
                }
            }
        }
        catch (std::exception const&)
        {
            LOG()->warn("Invalid read mapping for {} : {}", read.fragment_id(), read.graph_cigar());
        }
        return false; // node not covered by read
    };

    auto edgefilter = [&wgraph, nodefilter](Read& read, const std::string& node1, const std::string& node2) -> bool {
        try
        {
            graphs::GraphMapping mapping(read.graph_pos(), read.graph_cigar(), read.bases(), wgraph);

            const auto node_id1 = wgraph.nodeId(node1);
            const auto node_id2 = wgraph.nodeId(node2);

            const graphs::NodeMapping* previous_mapping = nullptr;
            for (size_t i = 0; i < mapping.size(); ++i)
            {
                auto const& nodemapping = mapping[i];
                const uint64_t this_node = nodemapping.node_id();

                if (previous_mapping != nullptr && previous_mapping->node_id() == node_id1 && this_node == node_id2)
                {
                    return previous_mapping->matched() >= (unsigned)std::min(previous_mapping->referenceSpan(), 16)
                        && nodemapping.matched() >= (unsigned)std::min(nodemapping.referenceSpan(), 16);
                }

                previous_mapping = &nodemapping;
            }
        }
        catch (std::exception const&)
        {
            LOG()->warn("Invalid read mapping for {} : {}", read.fragment_id(), read.graph_cigar());
        }
        return false; // edge not covered by read
    };

    disambiguateReads(wgraph, all_reads, nodefilter, edgefilter, parameters.description()["paths"]);

    countReads(
        coordinates, all_reads, output, parameters.output_enabled(Parameters::NODE_READ_COUNTS),
        parameters.output_enabled(Parameters::EDGE_READ_COUNTS),
        parameters.output_enabled(Parameters::PATH_READ_COUNTS),
        parameters.output_enabled(Parameters::DETAILED_READ_COUNTS));

    getVariants(
        coordinates, all_reads, output, parameters.min_reads_for_variant(), parameters.min_frac_for_variant(),
        parameters.description()["paths"], parameters.output_enabled(Parameters::VARIANTS),
        parameters.output_enabled(Parameters::NODE_COVERAGE), parameters.output_enabled(Parameters::PATH_COVERAGE));

    if (parameters.output_enabled(Parameters::ALIGNMENTS))
    {
        output_reads.reserve(all_reads.size() + output_reads.size());
        for (auto& r : all_reads)
        {
            Json::Value r_json = r->asJson();
            output["alignments"].append(r_json);
            output_reads.emplace_back(std::move(r));
        }
    }
    all_reads = std::move(output_reads);
    return output;
}

/**
 * Update sequence labels in read according to nodes the read has traversed
 * @param g graph structure
 * @param reads list of aligned reads
 * @param nodefilter filter to check if a read supports a particular node
 * @param edgefilter filter to check if a read supports a particular edge
 * @param paths
 */
void disambiguateReads(
    graphs::WalkableGraph& g, std::vector<common::p_Read>& reads, ReadSupportsNode nodefilter,
    ReadSupportsEdge edgefilter, Json::Value const& paths)
{
    // this variable tells us the ordered paths by sequence
    std::unordered_map<uint64_t, std::list<std::set<uint64_t>>> paths_by_sequence;

    // go through paths and order them by sequence
    for (auto const& p : paths)
    {
        auto sequence_name = p["sequence"].asString();
        uint64_t sequence_id = g.sequenceId(sequence_name);
        // sorted node_sequence -- this works as long as nodes are in topological order
        std::set<uint64_t> node_sequence;
        uint64_t previous_node = 0;
        for (const auto& n : p["nodes"])
        {
            const uint64_t this_node = g.nodeId(n.asString());
            if (!node_sequence.empty())
            {
                // check topological sorting
                assert(this_node > previous_node);
            }

            node_sequence.insert(this_node);
            previous_node = this_node;
        }
        paths_by_sequence[sequence_id].push_back(node_sequence);
    }

    for (auto& read : reads)
    {
        read->clear_graph_sequences_supported();
        read->clear_graph_nodes_supported();
        read->clear_graph_edges_supported();
        if (read->graph_mapping_status() == reads::MappingStatus::MAPPED)
        {
            bool has_previous = false;
            uint64_t previous_node = 0;
            std::string previous_node_name;

            std::set<std::pair<uint64_t, uint64_t>> edges_supported_by_read;
            std::set<uint64_t> nodes_supported_by_read;
            std::set<uint64_t> supported_sequences;

            GraphMapping gm(read->graph_pos(), read->graph_cigar(), read->bases(), g);

            for (size_t node_index = 0; node_index < gm.size(); ++node_index)
            {
                auto const& node_name = g.nodeName(gm[node_index].node_id());
                auto node = g.nodeId(node_name);

                if (has_previous && (edgefilter == nullptr || edgefilter(*read, previous_node_name, node_name)))
                {
                    edges_supported_by_read.emplace(previous_node, node);
                    for (uint64_t s : g.edge(previous_node, node)->sequence_ids())
                    {
                        supported_sequences.emplace(s);
                    }
                }
                has_previous = true;
                previous_node = node;
                previous_node_name = node_name;

                // check if node is rejected
                if (nodefilter != nullptr && !nodefilter(*read, node_name))
                {
                    continue;
                }

                nodes_supported_by_read.emplace(node);
                for (uint64_t s : g.node(node)->sequence_ids())
                {
                    supported_sequences.emplace(s);
                }
            }

            for (auto n : nodes_supported_by_read)
            {
                read->add_graph_nodes_supported(g.nodeName(n));
            }

            for (auto const& e : edges_supported_by_read)
            {
                read->add_graph_edges_supported(g.nodeName(e.first) + "_" + g.nodeName(e.second));
            }

            // A read must support at least one node.
            // Note that technically we could have a read that supports an edge but
            // not the two nodes next to it but practically we don't deal with this
            // case atm. since it's an uncommon scenario and we may not want to trust
            // such reads anyway.
            if (nodes_supported_by_read.empty())
            {
                continue;
            }

            for (uint64_t sequence_id : supported_sequences)
            {
                int paths_supported = 0;
                int paths_broken = 0;
                int paths_unsupported = 0;

                for (auto const& path : paths_by_sequence[sequence_id])
                {
                    std::list<uint64_t> common;
                    std::set_intersection(
                        path.begin(), path.end(), nodes_supported_by_read.begin(), nodes_supported_by_read.end(),
                        std::back_inserter(common));
                    if (common.empty())
                    {
                        // read path starts / ends before
                        ++paths_unsupported;
                    }
                    else
                    {
                        const uint64_t min_common = *(std::min_element(common.begin(), common.end()));
                        const uint64_t max_common = *(std::max_element(common.begin(), common.end()));
                        auto start_path = path.find(min_common);
                        auto start_read = nodes_supported_by_read.find(min_common);

                        if (start_path != path.begin() && start_read != nodes_supported_by_read.begin())
                        {
                            // Scenario where read crosses into path:
                            //
                            //   path ... --*-->*===*=== ...
                            //   read ... --*--/
                            //
                            // in this case the read doesn't support our path.
                            ++paths_broken;
                        }
                        else
                        {
                            int matched = 0;
                            has_previous = false;
                            previous_node = 0;
                            while (start_path != path.end() && start_read != nodes_supported_by_read.end()
                                   && *start_path <= max_common && *start_read <= max_common)
                            {
                                if (*start_path == *start_read)
                                {
                                    // needs to support all edges along the shared bit with path
                                    if (has_previous
                                        && edges_supported_by_read.count(std::make_pair(previous_node, *start_path))
                                            == 0)
                                    {
                                        paths_broken += 1;
                                        break;
                                    }
                                    has_previous = true;
                                    previous_node = *start_path;
                                    ++matched;
                                    ++start_path;
                                    ++start_read;
                                }
                                else
                                {
                                    break;
                                }
                            }
                            // one by one path and read follow the same nodes
                            // also, if the read crosses over at the end we consider the path broken
                            if (matched != ((int)common.size())
                                || (start_path != path.end() && start_read != nodes_supported_by_read.end()))
                            {
                                // two cases:
                                // either the read and path diverge after the end:
                                //
                                // ... ===*===>*---*--> read
                                //              \--*--> path
                                //
                                // or we have missed or added a node along the read / path.
                                // both cases mean that the read doesn't support the path.
                                paths_broken += 1;
                            }
                            else
                            {
                                // read and path either share a suffix / prefix, match or
                                // have one fully contained in the other
                                paths_supported += 1;
                            }
                        }
                    }
                }

                if (paths_supported > 0 && paths_broken == 0)
                {
                    read->add_graph_sequences_supported(
                        ((graphs::Graph const&)g).header->sequencenames((int)sequence_id));
                }
                else if (paths_broken > 0)
                {
                    read->add_graph_sequences_broken(((graphs::Graph const&)g).header->sequencenames((int)sequence_id));
                }
            }
        }
    }
}
}
