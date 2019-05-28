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
 * \brief Counts reads/fragments supporting different elements of the graph
 */

#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string/join.hpp>

#include "common/Fragment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "paragraph/ReadCounting.hh"

// Error.hh always needs to go last
#include "common/Error.hh"

using common::Fragment;
using common::FragmentList;
using common::Read;
using common::ReadBuffer;
using common::readsToFragments;
using graphtools::Graph;
using graphtools::GraphCoordinates;

// #define FRAGMENT_STATS_HISTOGRAM

namespace paragraph
{

void addFragmentCount(Json::Value& out, std::string const& element, Fragment const& frag)
{
    if (!out.isMember(element))
    {
        out[element] = 1;
        out[element + ":READS"] = frag.get_n_reads();
        out[element + ":FWD"] = frag.get_n_graph_forward_reads();
        out[element + ":REV"] = frag.get_n_graph_reverse_reads();
    }
    else
    {
        out[element] = out[element].asUInt64() + 1;
        out[element + ":READS"] = out[element + ":READS"].asUInt64() + frag.get_n_reads();
        out[element + ":FWD"] = out[element + ":FWD"].asUInt64() + frag.get_n_graph_forward_reads();
        out[element + ":REV"] = out[element + ":REV"].asUInt64() + frag.get_n_graph_reverse_reads();
    }
}

Json::Value countNodes(FragmentList const& fragments)
{
    Json::Value out = Json::ValueType::objectValue;
    for (auto const& frag : fragments)
    {
        for (const auto& n : frag->graph_nodes_supported())
        {
            addFragmentCount(out, n, *frag);
        }
    }
    return out;
}

Json::Value countEdges(FragmentList const& fragments)
{
    Json::Value out = Json::ValueType::objectValue;
    for (auto const& frag : fragments)
    {
        for (const auto& e : frag->graph_edges_supported())
        {
            addFragmentCount(out, e, *frag);
        }
    }
    return out;
}

Json::Value countPathFamilies(FragmentList const& fragments, bool detailed)
{
    Json::Value out = Json::ValueType::objectValue;
    for (auto const& frag : fragments)
    {
        if (!frag->graph_sequences_supported().empty())
        {
            std::vector<std::string> seqs;
            seqs.reserve(frag->graph_sequences_supported().size());
            seqs.insert(seqs.end(), frag->graph_sequences_supported().begin(), frag->graph_sequences_supported().end());
            std::sort(seqs.begin(), seqs.end());
            auto joined_sequence_name = boost::algorithm::join(seqs, ",");
            if (!out.isMember(joined_sequence_name))
            {
                out[joined_sequence_name] = Json::Value();
            }
            addFragmentCount(out[joined_sequence_name], "total", *frag);
            if (detailed) // Count Nodes/Edges within this path family
            {
                for (const auto& n : frag->graph_nodes_supported())
                {
                    addFragmentCount(out[joined_sequence_name], n, *frag);
                }
                for (const auto& e : frag->graph_edges_supported())
                {
                    addFragmentCount(out[joined_sequence_name], e, *frag);
                }
            }
        }
    }
    return out;
}

Json::Value alignmentStats(FragmentList const& fragments)
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

    Json::Value stats = Json::ValueType::objectValue;
    stats["mean_linear"] = mean(fragment_size);
    stats["mean_graph"] = mean(graph_fragment_size);
    stats["median_linear"] = median(fragment_size);
    stats["median_graph"] = median(graph_fragment_size);
    stats["variance_linear"] = variance(fragment_size);
    stats["variance_graph"] = variance(graph_fragment_size);
    stats["single_read"] = (Json::UInt64)single_read_fragments;
    stats["paired_read"] = (Json::UInt64)paired_read_fragments;
    stats["multi_read"] = (Json::UInt64)multi_read_fragments;
    stats["problematic_linear"] = (Json::UInt64)problematic_fragments_linear;
    stats["problematic_graph"] = (Json::UInt64)problematic_fragments_graph;

#ifdef FRAGMENT_STATS_HISTOGRAM
    stats["linear_histogram"] = Json::arrayValue;
    stats["graph_histogram"] = Json::arrayValue;

    histogram_type linear_hist = density(fragment_size);
    histogram_type graph_hist = density(graph_fragment_size);

    for (auto const& bin : linear_hist)
    {
        Json::Value bin_value = Json::objectValue;
        bin_value["lb"] = bin.first;
        bin_value["value"] = bin.second;
        stats["linear_histogram"].append(bin_value);
    }
    for (auto const& bin : graph_hist)
    {
        Json::Value bin_value = Json::objectValue;
        bin_value["lb"] = bin.first;
        bin_value["value"] = bin.second;
        stats["graph_histogram"].append(bin_value);
    }
#endif
    return stats;
}

void countReads(
    GraphCoordinates const& coordinates, ReadBuffer const& reads, Json::Value& output, bool by_node, bool by_edge,
    bool by_pathFam, bool pathFam_detailed)
{
    FragmentList fragments;
    readsToFragments(coordinates, reads, fragments);
    output["fragment_statistics"] = alignmentStats(fragments);
    if (by_node)
    {
        output["read_counts_by_node"] = countNodes(fragments);
    }
    if (by_edge)
    {
        output["read_counts_by_edge"] = countEdges(fragments);
    }
    if (by_pathFam)
    {
        output["read_counts_by_sequence"] = countPathFamilies(fragments, pathFam_detailed);
    }
}
}
