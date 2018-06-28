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
 * \brief Filter that detects whether reads were mapped uniquely using a graph kmer index
 *
 * \file KmerFilter.hh
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <unordered_set>

#include "KmerFilter.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/KmerIndexOperations.hh"

#include "common/Error.hh"

namespace paragraph
{
namespace readfilters
{
    using graphtools::Graph;
    using graphtools::GraphAlignment;
    using graphtools::NodeId;
    using graphtools::decodeGraphAlignment;
    struct KmerFilter::Impl
    {
        Impl(Graph const* g, int32_t kmer_len_)
            : graph(g)
            , graph_index(*graph, kmer_len_)
            , kmer_len(kmer_len_)
        {
        }
        Graph const* graph;
        graphtools::KmerIndex graph_index;
        int32_t kmer_len;
    };

    KmerFilter::KmerFilter(Graph const* graph, int32_t kmer_len)
    {
        if (kmer_len < 0)
        {
            kmer_len = graphtools::findMinCoveringKmerLength(
                graph, static_cast<size_t>(-kmer_len), static_cast<size_t>(-kmer_len));
            LOG()->info("Auto-detected kmer length is {}.", kmer_len);
        }
        _impl.reset(new Impl(graph, kmer_len));
    }

    KmerFilter::~KmerFilter() = default;

    std::pair<bool, std::string> KmerFilter::filterRead(common::Read const& r)
    {
        const GraphAlignment alignment = decodeGraphAlignment(r.graph_pos(), r.graph_cigar(), _impl->graph);
        if (alignment.size() < 1)
        {
            return { true, "kmer_nomapping" };
        }
        const auto sc_left = alignment[0].numClipped();
        const auto sc_right = alignment[alignment.size() - 1].numClipped();
        const auto& bases = r.bases();
        if ((signed)(bases.size() - sc_left - sc_right) < _impl->kmer_len)
        {
            return { true, "kmer_tooshort" };
        }

        std::unordered_set<std::string> kmers;
        for (size_t pos = sc_left; pos <= bases.size() - sc_right - _impl->kmer_len; ++pos)
        {
            kmers.insert(bases.substr(pos, static_cast<unsigned long>(_impl->kmer_len)));
        }

        std::unordered_set<NodeId> nodes_not_covered;
        std::list<NodeId> nodes_supported;
        for (int32_t node_index = 0; node_index != (int32_t)alignment.size(); ++node_index)
        {
            const auto node_id = static_cast<const NodeId>(alignment.getNodeIdByIndex(node_index));
            // only check the nodes which actually have unique overlapping kmers
            if (_impl->graph_index.numUniqueKmersOverlappingNode(node_id) > 0)
            {
                nodes_not_covered.insert(node_id);
                nodes_supported.push_back(node_id);
            }
        }

        for (auto const& kmer : kmers)
        {
            if (_impl->graph_index.numPaths(kmer) == 1)
            {
                auto paths = _impl->graph_index.getPaths(kmer);
                for (const auto& node_id : paths.front().nodeIds())
                {
                    nodes_not_covered.erase(node_id);
                    if (nodes_not_covered.empty())
                    {
                        return { false, "" };
                    }
                }
            }
        }

        std::string result_msg = "kmer_uncov";
        for (auto const& node : nodes_supported)
        {
            if (nodes_not_covered.count(node) != 0)
            {
                result_msg += "_";
                result_msg += std::to_string(node);
            }
        }

        return { true, result_msg };
    }
}
}
