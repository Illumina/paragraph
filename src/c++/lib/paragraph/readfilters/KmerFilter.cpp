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
#include "graphs/GraphMapping.hh"
#include "graphs/GraphMappingOperations.hh"
#include "graphs/KmerIndex.hh"
#include "graphs/KmerIndexOperations.hh"

#include "common/Error.hh"

namespace paragraph
{
namespace readfilters
{
    struct KmerFilter::Impl
    {
        Impl(std::shared_ptr<graphs::WalkableGraph> g, int32_t kmer_len_)
            : graph(std::move(g))
            , graph_index(graph, kmer_len_)
            , kmer_len(kmer_len_)
        {
        }
        std::shared_ptr<graphs::WalkableGraph> graph;
        graphs::KmerIndex graph_index;
        int32_t kmer_len;
    };

    KmerFilter::KmerFilter(graphs::WalkableGraph const* graph, int32_t kmer_len)
    {
        auto wgraph = std::make_shared<graphs::WalkableGraph>(*graph);
        if (kmer_len < 0)
        {
            kmer_len = graphs::findMinCoveringKmerLength(
                wgraph, static_cast<size_t>(-kmer_len), static_cast<size_t>(-kmer_len));
            LOG()->info("Auto-detected kmer length is {}.", kmer_len);
        }
        _impl.reset(new Impl(wgraph, kmer_len));
    }

    KmerFilter::~KmerFilter() = default;

    std::pair<bool, std::string> KmerFilter::filterRead(common::Read const& r)
    {
        const graphs::GraphMapping mapping
            = graphs::decodeFromString(r.graph_pos(), r.graph_cigar(), r.bases(), *_impl->graph);
        if (mapping.size() < 1)
        {
            return { true, "kmer_nomapping" };
        }
        const auto sc_left = mapping[0].mapping.clipped();
        const auto sc_right = mapping[mapping.size() - 1].mapping.clipped();
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

        std::unordered_set<uint64_t> nodes_not_covered;
        std::list<uint64_t> nodes_supported;
        for (const auto& nodemapping : mapping)
        {
            const auto id = static_cast<const uint64_t>(nodemapping.node_id);
            // only check the nodes which actually have unique overlapping kmers
            if (_impl->graph_index.numUniqueKmersOverlappingNode(static_cast<uint32_t>(id)) > 0)
            {
                nodes_not_covered.insert(id);
                nodes_supported.push_back(id);
            }
        }

        for (auto const& kmer : kmers)
        {
            if (_impl->graph_index.numPaths(kmer) == 1ull)
            {
                auto paths = _impl->graph_index.getPaths(kmer);
                for (const auto& node_id : paths.front().node_ids())
                {
                    nodes_not_covered.erase(static_cast<uint64_t>(node_id));
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
