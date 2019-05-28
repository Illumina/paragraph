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

#include "genotyping/BreakpointStatistics.hh"
#include <algorithm>
#include <boost/algorithm/string/join.hpp>

#include "common/Error.hh"

using graphtools::Graph;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace genotyping
{

/**
 * Create breakpoint statistics on graph
 * @param wgraph graph information
 * @param node_name name of node
 * @param forward true to use successor node/edges, false to use predecessors
 */
BreakpointStatistics::BreakpointStatistics(Graph const& graph, NodeId node_id, bool forward)
{
    const auto& node_name = graph.nodeName(node_id);
    const auto allele_nodes = forward ? graph.successors(node_id) : graph.predecessors(node_id);

    // note that sequences / sequence labels that have the same name as an edge
    // will not work well with this.
    assert(allele_nodes.size() > 1);

    std::map<std::string, std::set<std::string>> allele_edge_sets;
    for (auto const& an : allele_nodes)
    {
        const auto& an_name = graph.nodeName(an);
        const std::string edge_name = forward ? (node_name + "_" + an_name) : (an_name + "_" + node_name);

        // create index of counts <-> edge names
        edge_names.push_back(edge_name);
        edge_name_to_index[edge_name] = edge_names.size() - 1;

        const auto& edge_labels = forward ? graph.edgeLabels(node_id, an) : graph.edgeLabels(an, node_id);
        for (const auto& allele_name : edge_labels)
        {
            allele_edge_sets[allele_name].insert(edge_name);
            auto a_it = std::find(all_allele_names.begin(), all_allele_names.end(), allele_name);
            if (a_it == all_allele_names.end())
            {
                all_allele_names.push_back(allele_name);
            }
        }
    }

    // create canonical alleles
    std::map<std::string, std::list<std::string>> canonical_allele_to_allele;
    for (const auto& allele : allele_edge_sets)
    {
        const std::string canonical_allele_id = boost::algorithm::join(allele.second, ";");
        canonical_allele_to_allele[canonical_allele_id].push_back(allele.first);
    }

    std::vector<std::pair<std::string, std::string>> ordered_canonical_alleles;
    // choose representative allele from each equivalence class
    for (const auto& canonical_allele : canonical_allele_to_allele)
    {
        const auto ref_it = std::find(canonical_allele.second.begin(), canonical_allele.second.end(), "REF");
        auto canonical_allele_name = canonical_allele.second.front();
        if (ref_it != canonical_allele.second.end())
        {
            ordered_canonical_alleles.insert(ordered_canonical_alleles.begin(), std::make_pair("REF", "REF"));
            canonical_allele_name = "REF";
        }
        canonical_allele_names.push_back(canonical_allele_name);
        const auto this_allele_index = canonical_allele_names.size() - 1;
        for (const auto& edge : allele_edge_sets[canonical_allele_name])
        {
            edgename_to_alleles[edge].push_back(this_allele_index);
        }
        for (auto const& noncanonical_allele : canonical_allele.second)
        {
            allele_name_to_index[noncanonical_allele] = this_allele_index;
            allele_name_to_canonical_allele_name[noncanonical_allele] = canonical_allele_name;
            if (noncanonical_allele != "REF")
            {
                ordered_canonical_alleles.emplace_back(std::make_pair(canonical_allele_name, noncanonical_allele));
            }
        }
    }
}

/**
 * Add edge counts from paragraph output
 * @param paragraph_json JSON output from alignAndDisambiguate
 */
void BreakpointStatistics::addCounts(Json::Value const& paragraph_json)
{
    if (!paragraph_json.isMember("read_counts_by_edge"))
    {
        error("Cannot find key read_counts_by_edge in JSON");
    }
    for (auto const& edge_name : edge_names)
    {
        const auto e_it = edge_name_to_index.find(edge_name);
        assert(e_it != edge_name_to_index.end());
        const int this_edge_count = paragraph_json["read_counts_by_edge"].isMember(edge_name)
            ? paragraph_json["read_counts_by_edge"][edge_name].asInt()
            : 0;

        if (this_edge_count == 0)
        {
            continue;
        }

        if (edge_counts.size() <= e_it->second)
        {
            edge_counts.resize(edge_names.size(), 0);
            assert(edge_counts.size() > e_it->second);
        }

        edge_counts[e_it->second] += this_edge_count;

        // add counts for canonical alleles also
        for (const auto& allele : edgename_to_alleles[edge_name])
        {
            if (allele_counts.size() <= allele)
            {
                allele_counts.resize(canonical_allele_names.size(), 0);
                assert(allele_counts.size() > allele);
            }
            allele_counts[allele] += this_edge_count;
        }
    }
}

int32_t BreakpointStatistics::getCount(std::string const& edge_or_allele_name) const
{
    const auto e_it = edge_name_to_index.find(edge_or_allele_name);
    const auto a_it = allele_name_to_index.find(edge_or_allele_name);
    int32_t result = 0;
    if (e_it != edge_name_to_index.end() && a_it != allele_name_to_index.end())
    {
        error("Allele / sequence name %s is ambiguous with an edge name.", edge_or_allele_name.c_str());
    }
    else if (e_it != edge_name_to_index.end())
    {
        result = e_it->second >= edge_counts.size() ? 0 : edge_counts[e_it->second];
    }
    else if (a_it != allele_name_to_index.end())
    {
        result = a_it->second >= allele_counts.size() ? 0 : allele_counts[a_it->second];
    }
    else
    {
        // unknown edge or allele -- this is allowed since for a complex breakpoint set
        // we may not see all alleles at all breakpoints
        return 0;
    }
    return result;
}
}
