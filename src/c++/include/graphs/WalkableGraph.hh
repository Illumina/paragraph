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
 * \brief WalkableGraph class -- collect indexes necessary to walk a graph
 *
 * \file WalkableGraph.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "Graph.hh"

namespace graphs
{

typedef std::set<uint64_t> NeighborSet;
typedef std::unordered_map<uint64_t, NeighborSet> AdjacencyList;

/**
 * Graph class with adjacency lists
 */
class WalkableGraph
{
public:
    WalkableGraph();
    explicit WalkableGraph(Graph const& g);
    WalkableGraph(WalkableGraph const& g);
    WalkableGraph(WalkableGraph&& wg) noexcept;
    ~WalkableGraph();

    WalkableGraph& operator=(Graph const& g);
    WalkableGraph& operator=(WalkableGraph const& g);
    WalkableGraph& operator=(WalkableGraph&& g) noexcept;

    operator Graph&() const;

    /**
     * Nodes in topological order
     */
    std::list<uint64_t> const& allNodes() const;

    /**
     * Get predecessor of node
     */
    std::list<uint64_t> pred(uint64_t) const;

    /**
     * Get successor of node
     */
    std::list<uint64_t> succ(uint64_t) const;

    /**
     * get source and sink node.
     */
    uint64_t source() const;
    uint64_t sink() const;

    graphs::GraphHeader const& header() const;

    uint64_t sequenceId(const std::string& sequence_name) const;
    std::string sequenceName(uint64_t sequence_name) const;

    /**
     * by-name index
     */
    graphs::Node* node(uint64_t node_id);

    graphs::Node* node(const std::string& node_name);

    uint64_t nodeId(const std::string& node_name) const;
    std::string nodeName(uint64_t node_id) const;

    graphs::Node const* node(const std::string& node_name) const;

    graphs::Node const* node(uint64_t node_id) const;

    graphs::Edge* edge(uint64_t from, uint64_t to);

    graphs::Edge const* edge(uint64_t from, uint64_t to) const;

    bool hasEdge(uint64_t from, uint64_t to) const;

    graphs::Edge* edge(const std::string& node1, const std::string& node2);

    graphs::Edge const* edge(const std::string& node1, const std::string& node2) const;

private:
    struct WalkableGraphImpl;
    std::unique_ptr<WalkableGraphImpl> _impl;
};

/**
 * Create an adjacency list
 * @param g Graph structure
 * @return the adjacency list
 */
AdjacencyList makeAdjacencyList(Graph const& g);

/**
 * Create inverted adjacency list
 * @param g Graph structure
 * @return Adjacency list for graph with edges reversed
 */
AdjacencyList makeReverseAdjacencyList(Graph const& g);
};
