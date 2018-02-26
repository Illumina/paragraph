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

#include <boost/algorithm/string/join.hpp>
#include <numeric>

#include "common/HashHelper.hh"
#include "graphs/GraphCoordinates.hh"
#include "graphs/GraphMapping.hh"
#include "graphs/GraphMappingOperations.hh"
#include "graphs/GraphPathOperations.hh"
#include "paragraph/HaplotypePaths.hh"

using common::ReadBuffer;
using graphs::GraphPath;
using graphs::WalkableGraph;
using graphs::checkPathPrefixSuffixOverlap;
using graphs::decodeFromString;
using graphs::greedyMerge;
using graphs::mergePaths;
using std::list;
using std::string;
using std::unordered_map;

namespace paragraph
{

list<GraphPath> findHaplotypes(WalkableGraph const& wgraph, ReadBuffer const& reads)
{
    auto graph = std::make_shared<WalkableGraph>(wgraph);
    graphs::GraphCoordinates coordinates(wgraph);

    unordered_map<string, list<const common::Read*>> fragments;

    for (auto& read : reads)
    {
        const auto fragment_id = read->fragment_id();
        fragments[fragment_id].emplace_back(read.get());
    }

    // list fragments. every fragment contributes a path. Paths
    // are then extended to the start and end of the start / end node
    std::map<std::string, std::vector<int32_t>> path_map;
    for (const auto& fragment : fragments)
    {
        list<GraphPath> fragment_paths;
        for (const auto& read : fragment.second)
        {
            graphs::GraphMapping mapping
                = decodeFromString(read->graph_pos(), read->graph_cigar(), read->bases(), graph->graph());

            std::vector<int32_t> nodes;
            nodes.reserve(mapping.size());
            for (const auto& m : mapping)
            {
                nodes.push_back((int32_t)m.node_id);
            }
            GraphPath read_path{ graph, (int32_t)read->graph_pos(), nodes,
                                 (int32_t)(mapping.back().mapping.referenceSpan()) };

            fragment_paths.push_back(read_path);
        }
        greedyMerge(fragment_paths);

        for (const auto& path : fragment_paths)
        {
            if (path.num_nodes() == 0)
            {
                continue;
            }
            const std::string path_id = std::accumulate(
                std::next(path.node_ids().begin()), path.node_ids().end(), std::to_string(path.node_ids()[0]),
                [](const std::string& a, int b) { return a + '-' + std::to_string(b); });
            auto path_it = path_map.find(path_id);
            if (path_it == path_map.end())
            {
                path_map.emplace(std::make_pair(path_id, path.node_ids()));
            }
        }
    }

    list<GraphPath> result;
    for (const auto& path : path_map)
    {
        result.emplace_back(
            graph, 0, path.second,
            static_cast<int32_t>(wgraph.node(static_cast<uint64_t>(path.second.back()))->sequence().size() - 1));
    }

    graphs::exhaustiveMerge(result);

    return result;
}
}
