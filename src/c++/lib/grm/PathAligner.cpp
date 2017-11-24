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

#include "grm/PathAligner.hh"
#include "common/Genetics.hh"

#include "common/Error.hh"

namespace grm
{

struct PathAligner::PathAlignerImpl
{
    struct Path
    {
        std::string path_id;
        std::string sequence_id;
        std::string path_sequence;
        std::map<size_t, size_t> starts;
        std::vector<std::string> node_names;
        std::vector<uint64_t> node_ids;
    };

    std::list<Path> paths;
    std::unordered_map<uint64_t, size_t> node_lengths;
};

PathAligner::PathAligner()
    : _impl(new PathAlignerImpl())
{
}

PathAligner::~PathAligner() = default;

PathAligner::PathAligner(PathAligner&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

PathAligner& PathAligner::operator=(PathAligner&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

/**
 * Set the graph to align to
 * @param g a graph
 * @param paths list of paths
 */
void PathAligner::setGraph(graphs::Graph const& g, Json::Value const& paths)
{
    std::unordered_map<std::string, uint64_t> node_id_map;

    for (const auto& n : g.nodes)
    {
        assert(node_id_map.count(n.second->name()) == 0);
        node_id_map[n.second->name()] = n.first;
        _impl->node_lengths[n.first] = n.second->sequence().size();
    }

    for (const auto& p : paths)
    {
        assert(p.isMember("nodes"));
        assert(p.isMember("path_id"));
        assert(p.isMember("sequence"));

        PathAlignerImpl::Path path;
        path.path_id = p["path_id"].asString();
        path.sequence_id = p["sequence"].asString();
        for (const auto& n : p["nodes"])
        {
            const auto node_name = n.asString();
            const auto node_id = node_id_map[node_name];
            path.starts.emplace(path.path_sequence.size(), path.node_ids.size());
            auto node = g.nodes.find(node_id);
            assert(node != g.nodes.cend());
            path.path_sequence += node->second->sequence();
            path.node_ids.push_back(node_id);
            path.node_names.push_back(node_name);
        }
        _impl->paths.emplace_back(path);
    }
}

/**
 * Align a read to the graph and update the graph_* fields.
 *
 * We will align both the forward and the reverse strand and return
 * the alignment which is unique, or with the better score, or default
 * to the forward strand.
 *
 * @param read read structure
 */
void PathAligner::alignRead(common::Read& read) const
{
    bool is_mapped = false;
    for (const auto& path : _impl->paths)
    {
        auto match_pos_1 = path.path_sequence.find(read.bases());

        bool is_reverse_match = false;
        const auto rev_bases = common::reverseComplement(read.bases());
        if (match_pos_1 == std::string::npos)
        {
            match_pos_1 = path.path_sequence.find(rev_bases);
            if (match_pos_1 != std::string::npos)
            {
                is_reverse_match = true;
            }
        }

        if (match_pos_1 == std::string::npos)
        {
            // no exact match
            continue;
        }

        auto start_node = path.starts.lower_bound(match_pos_1);
        if (start_node == path.starts.end())
        {
            assert(path.starts.size() >= 1);
            start_node = std::prev(path.starts.end());
        }
        else if (start_node->first > match_pos_1)
        {
            start_node = std::prev(start_node);
        }
        assert(start_node != path.starts.end());

        const auto start = static_cast<google::protobuf::int32>(match_pos_1 - start_node->first);
        auto length_left = read.bases().size();
        auto this_start = (size_t)start;
        std::string cigar;
        int alignment_score = 0;
        std::list<uint64_t> graph_nodes_traversed;
        while (start_node != path.starts.end() && length_left > 0)
        {
            auto this_length = length_left;
            auto next_step = std::next(start_node);
            if (next_step != path.starts.end())
            {
                this_length = std::min(length_left, next_step->first - start_node->first - this_start);
            }

            if (this_length > 0)
            {
                const auto this_node_length = _impl->node_lengths[path.node_ids[start_node->second]];
                assert(this_start + this_length <= this_node_length);
                graph_nodes_traversed.push_back(start_node->second);
                cigar += std::to_string(path.node_ids[start_node->second]) + "[" + std::to_string(this_length) + "M]";
            }
            alignment_score += this_length;

            length_left -= this_length;
            start_node = next_step;
            this_start = 0;
        }

        if (is_mapped)
        {
            // check if this path maps our read somewhere else
            if ((is_reverse_match != read.is_reverse_strand()) != read.is_graph_reverse_strand()
                || start != read.graph_pos() || cigar != read.graph_cigar())
            {
                read.set_graph_mapq(0);
                read.set_is_graph_alignment_unique(false);
                break;
            }
            continue;
        }

        read.set_graph_pos(start);
        size_t reverse_strand_match_pos = std::string::npos;
        if (is_reverse_match)
        {
            // we don't need to update reverse_strand_match_pos: the match
            // already is on the reverse strand because we didn't find one forward
            // and second matches to the reverse strand are checked below.
            read.set_bases(rev_bases);
            read.set_is_graph_reverse_strand(!read.is_reverse_strand());
        }
        else
        {
            // check if read also matches on the reverse strand
            reverse_strand_match_pos = path.path_sequence.find(rev_bases);
            read.set_is_graph_reverse_strand(read.is_reverse_strand());
        }
        read.set_graph_cigar(cigar);
        read.set_graph_alignment_score(alignment_score);

        const auto match_pos_2 = path.path_sequence.find(read.bases(), match_pos_1 + 1);
        if (match_pos_2 != std::string::npos || reverse_strand_match_pos != std::string::npos)
        {
            read.set_graph_mapq(0);
            read.set_is_graph_alignment_unique(false);
        }
        else
        {
            read.set_graph_mapq(60);
            read.set_is_graph_alignment_unique(true);
        }
        read.set_graph_mapping_status(reads::MAPPED);
        is_mapped = true;
    }
}
}
