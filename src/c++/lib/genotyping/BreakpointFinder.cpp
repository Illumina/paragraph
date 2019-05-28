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
 * Functions to find breakpoints in a graph
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/BreakpointFinder.hh"

#include <string>
#include <vector>

#include "common/Error.hh"

using std::make_pair;
using std::map;
using std::string;
using std::vector;

namespace genotyping
{

using graphtools::NodeId;

/**
 * Create a map of all breakpoints in a graph
 * @param wgraph the graph
 * @return a breakpoint map for the graph
 */
BreakpointMap createBreakpointMap(graphtools::Graph const& wgraph)
{
    BreakpointMap breakpoint_map;
    const auto source_node = static_cast<NodeId>(0);
    const auto sink_node = static_cast<NodeId>(wgraph.numNodes() - 1);
    const bool has_source_and_sink = wgraph.nodeName(source_node) == "source" && wgraph.nodeName(sink_node) == "sink";

    for (auto node = source_node; node <= sink_node; ++node)
    {
        if (has_source_and_sink && (node == source_node || node == sink_node))
        {
            continue;
        }

        const string& node_name = wgraph.nodeName(node);
        if (wgraph.successors(node).size() > 1)
        {
            breakpoint_map.emplace(node_name + "_", BreakpointStatistics(wgraph, node, true));
        }

        if (wgraph.predecessors(node).size() > 1)
        {
            breakpoint_map.emplace(string("_") + node_name, BreakpointStatistics(wgraph, node, false));
        }
    }
    return breakpoint_map;
}
}
