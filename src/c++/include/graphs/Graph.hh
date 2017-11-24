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
 * \brief Graph class to collect all information for a single graph, protobuf helpers, ...
 *
 * \file Graph.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "graph.pb.h"

#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

#include "json/json.h"

namespace graphs
{
typedef std::unique_ptr<graphs::GraphHeader> p_Header;
typedef std::shared_ptr<graphs::Node> p_Node;
typedef std::shared_ptr<graphs::Edge> p_Edge;

struct Graph
{
    p_Header header;
    std::map<uint64_t, p_Node> nodes;
    std::vector<p_Edge> edges;
};

/**
 * Initialize graph from JSON
 * @param in Input JSON node
 * @param reference reference fasta name
 * @param out Output graph to initialize
 * @param store_ref_sequence store sequence for reference nodes
 */
void fromJson(Json::Value const& in, std::string const& reference, Graph& out, bool store_ref_sequence = true);

/**
 * Reverse a DAG
 * @param g input graph
 * @param g_rev reversed graph
 */
void reverse(const Graph& g, Graph& g_rev);
};
