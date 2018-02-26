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

/**
 * Create a map of all breakpoints in a graph
 * @param wgraph the graph
 * @return a breakpoint map for the graph
 */
BreakpointMap createBreakpointMap(graphs::WalkableGraph const& wgraph)
{
    BreakpointMap breakpoint_map;
    const bool has_source_and_sink
        = wgraph.nodeName(wgraph.source()) == "source" && wgraph.nodeName(wgraph.sink()) == "sink";

    for (auto node : wgraph.allNodes())
    {
        if (has_source_and_sink && (node == wgraph.source() || node == wgraph.sink()))
        {
            continue;
        }

        const string node_name = wgraph.nodeName(node);
        try
        {
            breakpoint_map.emplace(node_name + "_", BreakpointStatistics(wgraph, node_name, true));
        }
        catch (std::exception const&)
        {
            // ignore exception, only create breakpoints for nodes with multiple out-edges
        }
        try
        {
            breakpoint_map.emplace(string("_") + node_name, BreakpointStatistics(wgraph, node_name, false));
        }
        catch (std::exception const&)
        {
            // ignore exception, only create breakpoints for nodes with multiple out-edges
        }
    }
    return breakpoint_map;
}
}
