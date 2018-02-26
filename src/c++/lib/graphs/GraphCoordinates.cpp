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

#include "graphs/GraphCoordinates.hh"
#include "common/HashHelper.hh"

namespace graphs
{

struct GraphCoordinates::GraphCoordinatesImpl
{
    explicit GraphCoordinatesImpl(WalkableGraph const& g)
        : graph(g)
    {
        auto const& nodes_to_walk = g.allNodes();
        uint64_t canonical_offset = 0;
        for (auto n_id : nodes_to_walk)
        {
            auto const& n_name = g.nodeName(n_id);
            canonical_offsets[n_name] = canonical_offset;
            node_starts[canonical_offset] = n_name;
            canonical_offset += std::max((size_t)1, g.node(n_id)->sequence().size());

            // nodes are sorted in topological order, so we can compute distances as min over all predecessors
            for (auto n_source : nodes_to_walk)
            {
                // distance = zero in these cases
                if (n_id == n_source || g.hasEdge(n_source, n_id))
                {
                    continue;
                }

                size_t min_dist = std::numeric_limits<size_t>::max();
                for (auto pred : g.pred(n_id))
                {
                    auto pred_distance_it = node_end_to_start_distance.find(std::make_pair(n_source, pred));
                    if (pred_distance_it != node_end_to_start_distance.end())
                    {
                        // minimal distance via that predecessor
                        min_dist = std::min(min_dist, pred_distance_it->second + g.node(pred)->sequence().size());
                    }
                    else if (g.hasEdge(n_source, pred))
                    {
                        min_dist = std::min(min_dist, g.node(pred)->sequence().size());
                    }
                }

                if (min_dist != std::numeric_limits<size_t>::max())
                {
                    node_end_to_start_distance[std::make_pair(n_source, n_id)] = min_dist;
                }
            }
        }
    }

    WalkableGraph graph;

    std::unordered_map<std::string, uint64_t> canonical_offsets;
    std::map<uint64_t, std::string> node_starts;

    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> node_end_to_start_distance;
};

GraphCoordinates::GraphCoordinates(WalkableGraph const& g)
    : _impl(new GraphCoordinatesImpl(g))
{
}
GraphCoordinates::~GraphCoordinates() = default;

GraphCoordinates::GraphCoordinates(GraphCoordinates&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

GraphCoordinates& GraphCoordinates::operator=(GraphCoordinates&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

/**
 * Get a "canonical" / linearized position -- every base on the graph has such a position
 * @param node node name
 * @param offset offset relative to start of node
 * @return canonical position
 */
uint64_t GraphCoordinates::canonicalPos(std::string const& node, uint64_t offset) const
{
    auto ioffset = _impl->canonical_offsets.find(node);
    assert(ioffset != _impl->canonical_offsets.end());
    return ioffset->second + offset;
}

/**
 * Calculated canonical start and end positions for a graph mapping
 * @param mapping
 * @return start and end
 */
std::pair<uint64_t, uint64_t> GraphCoordinates::canonicalStartAndEnd(graphs::GraphMapping const& mapping) const
{
    std::pair<uint64_t, uint64_t> start_end{ -1, -1 };

    start_end.first = canonicalPos(
        _impl->graph.nodeName(mapping.front().node_id),
        static_cast<uint64_t>(mapping.front().mapping.reference_start()));

    auto end_offset = static_cast<uint64_t>(mapping.back().mapping.referenceSpan());
    if (mapping.size() > 0 && end_offset > 0)
    {
        if (mapping.size() == 1)
        {
            end_offset += mapping.front().mapping.reference_start();
        }
        start_end.second = canonicalPos(_impl->graph.nodeName(mapping.back().node_id), end_offset - 1);
    }

    if (start_end.first > start_end.second)
    {
        std::swap(start_end.first, start_end.second);
    }

    return start_end;
}

/**
 * Reverse lookup : get node and offset from a canonical pos
 * @param canonical_pos canonical position
 * @param node output variable for node name
 * @param offset output variable for offset
 */
void GraphCoordinates::nodeAndOffset(uint64_t canonical_pos, std::string& node, uint64_t& offset) const
{
    auto lb = _impl->node_starts.lower_bound(canonical_pos);
    if (lb != _impl->node_starts.end())
    {
        if (lb != _impl->node_starts.begin() && canonical_pos < lb->first)
        {
            lb = std::prev(lb);
        }
        node = lb->second;
        offset = canonical_pos - lb->first;
    }
    else
    {
        node = _impl->node_starts.rbegin()->second;
        offset = canonical_pos - _impl->node_starts.rbegin()->first;
    }
}

/**
 * Calculate the minimum distance in bp between two canonical positions
 * @param pos1 start pos
 * @param pos2 end pos
 * @return basepairs between pos1 and pos2
 */
uint64_t GraphCoordinates::distance(uint64_t pos1, uint64_t pos2) const
{
    if (pos1 == pos2)
    {
        return 0;
    }
    if (pos2 < pos1)
    {
        std::swap(pos1, pos2);
    }

    std::string n1, n2;
    uint64_t offset1, offset2;
    nodeAndOffset(pos1, n1, offset1);
    nodeAndOffset(pos2, n2, offset2);

    // on on the same node-> can compute distance directly
    if (n1 == n2)
    {
        return pos2 - pos1;
    }

    const uint64_t n1_id = _impl->graph.nodeId(n1);
    const uint64_t n2_id = _impl->graph.nodeId(n2);
    const size_t n1_length = _impl->graph.node(n1_id)->sequence().size();

    uint64_t result = std::numeric_limits<uint64_t>::max();
    if (_impl->graph.hasEdge(n1_id, n2_id))
    {
        result = n1_length - offset1 + offset2;
    }
    else
    {
        auto dist_it = _impl->node_end_to_start_distance.find(std::make_pair(n1_id, n2_id));

        if (dist_it != _impl->node_end_to_start_distance.end())
        {
            result = n1_length - offset1 + offset2 + dist_it->second;
        }
    }
    return result;
}

/**
 * @return the graph for these coordinates
 */
WalkableGraph const& GraphCoordinates::getGraph() const { return _impl->graph; }
}