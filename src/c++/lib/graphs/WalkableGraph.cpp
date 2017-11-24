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

#include "graphs/WalkableGraph.hh"

#include "common/HashHelper.hh"

#include <algorithm>

namespace graphs
{

struct WalkableGraph::WalkableGraphImpl
{
    explicit WalkableGraphImpl(Graph const& g)
    {
        if (g.header)
        {
            graph.header = std::unique_ptr<GraphHeader>(new GraphHeader(*g.header));
            for (int s_id = 0; s_id < graph.header->sequencenames().size(); ++s_id)
            {
                sequence_name_to_id[graph.header->sequencenames(s_id)] = static_cast<uint64_t>(s_id);
            }
        }

        std::vector<uint64_t> node_ids;
        for (auto const& n : g.nodes)
        {
            graph.nodes.emplace(n.first, p_Node(new Node(*n.second)));
            node_ids.push_back(n.first);
            name_to_id[n.second->name()] = n.first;
        }
        std::sort(node_ids.begin(), node_ids.end());
        nodes_in_topological_order.insert(nodes_in_topological_order.end(), node_ids.begin(), node_ids.end());

        for (auto const& e : g.edges)
        {
            graph.edges.emplace_back(p_Edge(new Edge(*e)));
            edge_index[std::make_pair(e->from(), e->to())] = graph.edges.back().get();
        }

        adj = makeAdjacencyList(g);
        reverse_adj = makeReverseAdjacencyList(g);
    }

    Graph graph;
    AdjacencyList adj;
    AdjacencyList reverse_adj;

    std::unordered_map<std::string, uint64_t> name_to_id;
    std::unordered_map<std::pair<uint64_t, uint64_t>, graphs::Edge*> edge_index;
    std::unordered_map<std::string, uint64_t> sequence_name_to_id;
    std::list<uint64_t> nodes_in_topological_order;
};

WalkableGraph::WalkableGraph() = default;

WalkableGraph::WalkableGraph(Graph const& g)
    : _impl(new WalkableGraphImpl(g))
{
}

WalkableGraph::WalkableGraph(WalkableGraph const& g)
    : _impl(new WalkableGraphImpl(g))
{
}

WalkableGraph::WalkableGraph(WalkableGraph&& g) noexcept
    : _impl(std::move(g._impl))
{
}

WalkableGraph::~WalkableGraph() = default;

WalkableGraph& WalkableGraph::operator=(WalkableGraph&& g) noexcept
{
    _impl = std::move(g._impl);
    return *this;
}

WalkableGraph& WalkableGraph::operator=(WalkableGraph const& g)
{
    if (&g == this)
    {
        return *this;
    }
    _impl.reset(new WalkableGraphImpl(g));
    return *this;
}

WalkableGraph& WalkableGraph::operator=(Graph const& g)
{
    _impl.reset(new WalkableGraphImpl(g));
    return *this;
}

WalkableGraph::operator Graph&() const { return _impl->graph; }

graphs::GraphHeader const& WalkableGraph::header() const { return *_impl->graph.header; }

uint64_t WalkableGraph::sequenceId(const std::string& sequence_name) const
{
    auto pos = _impl->sequence_name_to_id.find(sequence_name);
    assert(pos != _impl->sequence_name_to_id.end());
    return pos->second;
}

std::string WalkableGraph::sequenceName(uint64_t sequence_name) const
{
    return _impl->graph.header->sequencenames((int)sequence_name);
}

/**
 * Nodes in topological order
 */
std::list<uint64_t> const& WalkableGraph::allNodes() const { return _impl->nodes_in_topological_order; }

uint64_t WalkableGraph::source() const { return _impl->nodes_in_topological_order.front(); }

uint64_t WalkableGraph::sink() const { return _impl->nodes_in_topological_order.back(); }

/**
 * Get predecessor of node
 */
std::list<uint64_t> WalkableGraph::pred(uint64_t node) const
{
    std::list<uint64_t> result;
    auto it = _impl->reverse_adj.find(node);
    if (it != _impl->reverse_adj.end())
    {
        for (auto p : it->second)
        {
            result.push_back(p);
        }
    }
    return result;
}

/**
 * Get successor of node
 */
std::list<uint64_t> WalkableGraph::succ(uint64_t node) const
{
    std::list<uint64_t> result;
    auto it = _impl->adj.find(node);
    if (it != _impl->adj.end())
    {
        for (auto p : it->second)
        {
            result.push_back(p);
        }
    }
    return result;
}

graphs::Node* WalkableGraph::node(uint64_t node_id)
{
    assert(_impl->graph.nodes.count(node_id));
    return _impl->graph.nodes[node_id].get();
}

graphs::Node const* WalkableGraph::node(uint64_t node_id) const
{
    assert(_impl->graph.nodes.count(node_id));
    auto n_it = _impl->graph.nodes.find(node_id);
    return n_it->second.get();
}

graphs::Node* WalkableGraph::node(const std::string& node_name)
{
    assert(_impl->name_to_id.count(node_name));
    return node(_impl->name_to_id[node_name]);
}

uint64_t WalkableGraph::nodeId(const std::string& node_name) const
{
    assert(_impl->name_to_id.count(node_name));
    return _impl->name_to_id[node_name];
}

std::string WalkableGraph::nodeName(uint64_t node_id) const
{
    assert(node_id < _impl->graph.nodes.size());
    return _impl->graph.nodes[node_id]->name();
}

graphs::Node const* WalkableGraph::node(const std::string& node_name) const
{
    assert(_impl->name_to_id.count(node_name));
    return node(_impl->name_to_id[node_name]);
}

graphs::Edge* WalkableGraph::edge(uint64_t from, uint64_t to)
{
    const auto key = std::make_pair(from, to);
    assert(_impl->edge_index.count(key));
    return _impl->edge_index[key];
}

graphs::Edge const* WalkableGraph::edge(uint64_t from, uint64_t to) const
{
    const auto key = std::make_pair(from, to);
    assert(_impl->edge_index.count(key));
    return _impl->edge_index[key];
}

bool WalkableGraph::hasEdge(uint64_t from, uint64_t to) const
{
    const auto key = std::make_pair(from, to);
    return static_cast<bool>(_impl->edge_index.count(key));
}

graphs::Edge* WalkableGraph::edge(const std::string& node1, const std::string& node2)
{
    assert(_impl->name_to_id.count(node1));
    assert(_impl->name_to_id.count(node2));
    const auto from = _impl->name_to_id[node1];
    const auto to = _impl->name_to_id[node2];
    return edge(from, to);
}

graphs::Edge const* WalkableGraph::edge(const std::string& node1, const std::string& node2) const
{
    assert(_impl->name_to_id.count(node1));
    assert(_impl->name_to_id.count(node2));
    const auto from = _impl->name_to_id[node1];
    const auto to = _impl->name_to_id[node2];
    return edge(from, to);
}

AdjacencyList makeAdjacencyList(Graph const& g)
{
    AdjacencyList al;
    for (auto const& e : g.edges)
    {
        auto node_it = al.find(e->from());
        if (node_it == al.end())
        {
            al.emplace(e->from(), NeighborSet{ e->to() });
        }
        else
        {
            node_it->second.insert(e->to());
        }
    }
    return al;
}

AdjacencyList makeReverseAdjacencyList(Graph const& g)
{
    AdjacencyList al;
    for (auto const& e : g.edges)
    {
        auto node_it = al.find(e->to());
        if (node_it == al.end())
        {
            al.emplace(e->to(), NeighborSet{ e->from() });
        }
        else
        {
            node_it->second.insert(e->from());
        }
    }
    return al;
}
}
