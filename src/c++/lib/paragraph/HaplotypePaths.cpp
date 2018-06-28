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
 * \summary Find haplotypes from mapped reads
 *
 * \file HaplotypePaths.cpp
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#include <numeric>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphCoordinates.hh"
#include "graphcore/PathFamily.hh"
#include "graphcore/PathFamilyOperations.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/PairHashing.hh"
#include "paragraph/HaplotypePaths.hh"

#include "common/Error.hh"

using common::ReadBuffer;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::GraphCoordinates;
using graphtools::NodeId;
using graphtools::Path;
using graphtools::checkPathPrefixSuffixOverlap;
using graphtools::decodeGraphAlignment;
using graphtools::exhaustiveMerge;
using graphtools::greedyMerge;
using graphtools::mergePaths;
using std::list;
using std::string;
using std::unordered_map;

namespace paragraph
{

// Returns a string with all nodeIDs in a path
string pathIds(Path const& path)
{
    auto const& nodes = path.nodeIds();
    return std::accumulate(
        std::next(nodes.begin()), nodes.end(), std::to_string(nodes[0]),
        [](std::string const& a, int b) { return a + '-' + std::to_string(b); });
}

/**
 * Return the alignment path(s) for each aligned fragment
 * Merges overlapping mates
 * @param graph The graph
 * @param reads Aligned reads
 * @returns map fragmentID -> List of alignment paths for that fragment
 */
std::map<std::string, std::list<Path>> getFragmentPaths(Graph const& graph, ReadBuffer const& reads)
{
    std::map<std::string, std::list<Path>> path_map;
    // list fragments. every fragment contributes a path. Paths
    // are then extended to the start and end of the start / end node
    for (auto const& read : reads)
    {
        GraphAlignment mapping = decodeGraphAlignment(read->graph_pos(), read->graph_cigar(), &graph);
        if (mapping.path().numNodes() > 0)
        {
            path_map[read->fragment_id()].push_back(mapping.path());
        }
    }
    for (auto& fragment_paths : path_map)
    {
        greedyMerge(fragment_paths.second);
    }
    return path_map;
}

// Creates a path family with all edges from all paths
graphtools::PathFamily pathsToFamily(Graph* const graph, std::list<Path> paths)
{
    graphtools::PathFamily family(graph);
    for (auto const& path : paths)
    {
        for (auto start = path.begin(); start != std::prev(path.end()); ++start)
        {
            auto end = std::next(start);
            family.addEdge(*start, *end);
        }
    }
    return family;
}

typedef std::vector<graphtools::NodeIdPair> Edges;
// Sorted vector with all edges in the path family
Edges edgesVector(graphtools::PathFamily const& fam)
{
    Edges sortedEdges(fam.edges().begin(), fam.edges().end());
    std::sort(sortedEdges.begin(), sortedEdges.end());
    return sortedEdges;
}

std::vector<PhasingFamily> getPhasingFamilies(Graph* const graph, ReadBuffer const& reads)
{
    // TODO Could add a hash function for PathFamily instead of map with Edges
    std::map<Edges, PhasingFamily> phasingFams;
    GraphCoordinates coordinates(graph);
    for (auto const& pair : getFragmentPaths(*graph, reads))
    {
        graphtools::PathFamily family(pathsToFamily(graph, pair.second));

        // don't add empty families
        if (family.edges().empty())
        {
            continue;
        }

        const auto edges = edgesVector(family);

        graphtools::NodeId prev = 0;
        bool has_prev = false;
        bool is_linear = true;
        for (const auto& edge : edges)
        {
            if (has_prev
                && coordinates.distance(
                       coordinates.canonicalPos(graph->nodeName(prev), 0),
                       coordinates.canonicalPos(graph->nodeName(edge.first), 0))
                    == std::numeric_limits<uint64_t>::max())
            {
                is_linear = false;
                LOG()->trace("Family {} is not linear.", family);
            }
            prev = edge.second;
            has_prev = true;
        }

        if (!is_linear)
        {
            continue;
        }

        auto fam = phasingFams.find(edges);
        if (fam == phasingFams.end())
        {
            phasingFams.emplace(edgesVector(family), PhasingFamily{ family, 1 });
        }
        else
        {
            fam->second.second++;
        }
    }

    std::vector<PhasingFamily> result;
    result.reserve(phasingFams.size());
    for (auto const& kvp : phasingFams)
    {
        result.push_back(kvp.second);
    }
    return result;
}

void addHaplotypePaths(
    common::ReadBuffer const& reads, graphtools::Graph& graph, Json::Value& paths, Json::Value& output)
{
    // Compute and output phasing families
    Json::Value phasing = Json::arrayValue;
    auto families = getPhasingFamilies(&graph, reads);
    graphtools::PathFamily uber_family(&graph);
    for (auto& family : families)
    {
        Json::Value json_fam = Json::objectValue;
        json_fam["edges"] = Json::arrayValue;
        for (const auto& edge : family.first.edges())
        {
            Json::Value json_edge = Json::objectValue;
            json_edge["from"] = graph.nodeName(edge.first);
            json_edge["to"] = graph.nodeName(edge.second);
            json_fam["edges"].append(json_edge);
            uber_family.addEdge(edge.first, edge.second);
        }
        json_fam["count"] = family.second;
        phasing.append(json_fam);
    }
    output["phasing"] = phasing;

    struct HaplotypeGroup
    {
        std::list<size_t> paths;
        NodeId start = 0;
        NodeId end = 0;
    };

    // each path segment is part of one haplotype group
    std::vector<HaplotypeGroup> groups;
    std::vector<graphtools::Path> path_segments;
    {
        const std::list<graphtools::Path> path_segment_list = graphtools::getPathSegmentsForFamily(uber_family);
        path_segments = { path_segment_list.begin(), path_segment_list.end() };
        std::sort(path_segments.begin(), path_segments.end(), [](Path const& p1, Path const& p2) -> bool {
            return p1.nodeIds().front() < p2.nodeIds().front();
        });

        std::map<NodeId, std::list<size_t>> starts;
        for (const auto ps : path_segments | boost::adaptors::indexed(0))
        {
            starts[ps.value().nodeIds().front()].push_back(ps.index());
        }

        bool has_group = false;
        size_t current_group = 0;
        const std::function<NodeId(size_t)> end_node
            = [&path_segments](size_t ix) -> NodeId { return path_segments[ix].nodeIds().back(); };
        for (const auto& s : starts)
        {
            if (has_group && groups[current_group].end <= s.first)
            {
                has_group = false;
                ++current_group;
            }
            if (!has_group)
            {
                groups.resize(current_group + 1);
                auto& current = groups.back();
                current.start = s.first;
                current.end = *(boost::range::max_element(s.second | boost::adaptors::transformed(end_node)));
                has_group = true;
            }
            auto& current = groups.back();
            current.end = std::max(
                current.end, *(boost::range::max_element(s.second | boost::adaptors::transformed(end_node))));
            current.paths.insert(current.paths.end(), s.second.begin(), s.second.end());
        }
    }

    auto this_hap_group = groups.begin();
    auto next_hap_group = std::next(this_hap_group);
    while (this_hap_group != groups.end() && next_hap_group != groups.end())
    {
        bool has_merged = false;
        std::list<Path> group_merge_paths;
        for (const auto p1 : this_hap_group->paths)
        {
            bool can_merge = true;
            std::list<Path> p1_merge_paths;
            for (const auto p2 : next_hap_group->paths)
            {
                Path const& pp1 = path_segments[p1];
                Path const& pp2 = path_segments[p2];
                if (pp1.nodeIds().back() == pp2.nodeIds().front())
                {
                    const Path merged = mergePaths(pp1, pp2);
                    for (const auto& f : families)
                    {
                        if (f.first.containsPath(pp1) && f.first.containsPath(pp2) && f.first.containsPath(merged))
                        {
                            p1_merge_paths.push_back(merged);
                            break;
                        }
                    }
                }
                else
                {
                    can_merge = false;
                    break;
                }
            }
            if (!can_merge)
            {
                has_merged = false;
                break;
            }
            if (!p1_merge_paths.empty()
                && (this_hap_group->paths.size() == 1 || next_hap_group->paths.size() == 1
                    || p1_merge_paths.size() == 1))
            {
                group_merge_paths.insert(group_merge_paths.end(), p1_merge_paths.begin(), p1_merge_paths.end());
                has_merged = true;
            }
            else
            {
                has_merged = false;
                break;
            }
        }
        if (has_merged)
        {
            // replace group paths by merged paths -- note each path is in only one group
            // this is always positive because we merge all paths in the next group with at least one path in this group
            const size_t path_count_difference
                = this_hap_group->paths.size() + next_hap_group->paths.size() - group_merge_paths.size();
            const auto first_deleted_path = std::min(
                *boost::range::min_element(this_hap_group->paths), *boost::range::min_element(next_hap_group->paths));
            auto delete_begin = std::next(path_segments.begin(), first_deleted_path);
            auto delete_end = std::next(
                path_segments.begin(),
                std::max(
                    *boost::range::max_element(this_hap_group->paths),
                    *boost::range::max_element(next_hap_group->paths))
                    + 1);
            path_segments.erase(delete_begin, delete_end);
            path_segments.insert(
                std::next(path_segments.begin(), first_deleted_path), group_merge_paths.begin(),
                group_merge_paths.end());

            HaplotypeGroup new_hg;
            const std::function<NodeId(Path const&)> start_node
                = [](Path const& p) -> NodeId { return p.nodeIds().front(); };
            const std::function<NodeId(Path const&)> end_node
                = [](Path const& p) -> NodeId { return p.nodeIds().back(); };
            new_hg.start = *(boost::range::min_element(group_merge_paths | boost::adaptors::transformed(start_node)));
            new_hg.end = *(boost::range::max_element(group_merge_paths | boost::adaptors::transformed(end_node)));
            for (size_t j = first_deleted_path; j < first_deleted_path + group_merge_paths.size(); ++j)
            {
                new_hg.paths.push_back(j);
            }

            auto hg_pos = std::distance(groups.begin(), this_hap_group);
            ++next_hap_group;
            groups.erase(this_hap_group, next_hap_group);
            this_hap_group = std::next(groups.begin(), hg_pos);
            groups.insert(this_hap_group, new_hg);
            next_hap_group = std::next(this_hap_group);
            while (next_hap_group != groups.end())
            {
                for (auto& p : next_hap_group->paths)
                {
                    p -= path_count_difference;
                }
                ++next_hap_group;
            }
            this_hap_group = groups.begin();
            next_hap_group = std::next(this_hap_group);
        }
        else
        {
            ++this_hap_group;
            ++next_hap_group;
        }
    }

    std::vector<int> path_ix(path_segments.size(), 0);
    for (const auto ps : path_segments | boost::adaptors::indexed(0))
    {
        const auto& path_segment = ps.value();

        Json::Value haplo_json = Json::objectValue;
        haplo_json["path_length"] = (Json::UInt64)path_segment.length();
        haplo_json["path_start"] = path_segment.startPosition();
        haplo_json["path_end"] = path_segment.endPosition();
        haplo_json["path_encoding"] = path_segment.encode();
        haplo_json["nodes"] = Json::arrayValue;
        string path_id;
        for (const auto& node : path_segment.nodeIds())
        {
            haplo_json["nodes"].append(graph.nodeName(node));
            if (!path_id.empty())
            {
                path_id += "_";
            }
            path_id += std::to_string(node);
        }
        haplo_json["path_id"] = path_id;
        path_ix[ps.index()] = (int)paths.size();
        paths.append(haplo_json);
    }
    output["paths"] = paths;

    Json::Value hap_groups = Json::arrayValue;
    for (const auto& hg : groups)
    {
        Json::Value group = Json::objectValue;
        group["start_node"] = hg.start;
        group["end_node"] = hg.end;
        group["paths"] = Json::arrayValue;
        for (const auto& p : hg.paths)
        {
            group["paths"].append(paths[path_ix[p]]["path_id"]);
        }
        hap_groups.append(group);
    }

    output["phased_path_groups"] = hap_groups;
}
}
