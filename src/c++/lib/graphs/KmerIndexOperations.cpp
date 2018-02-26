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
 * \brief Kmer index for graphs
 *
 * \file KmerIndexOperations.cpp
 * \author Egor Dolzhenko & Peter Krusche
 * \email edolzhenko@illumina.com, pkrusche@illumina.com
 *
 */

#include "graphs/KmerIndexOperations.hh"

namespace graphs
{

/**
 * Find minimum kmer length that covers each node with a unique kmer
 * @param graph  a graph
 * @param min_unique_kmers_per_edge  min number of unique kmers to cover each edge
 * @param min_unique_kmers_per_node  min number of unique kmers to cover each node
 * @return
 */
int findMinCoveringKmerLength(
    std::shared_ptr<WalkableGraph>& graph, size_t min_unique_kmers_per_edge, size_t min_unique_kmers_per_node)
{
    for (int k = 10; k < 64; ++k)
    {
        KmerIndex index(graph, k);

        bool any_below = false;
        for (const auto& n : (static_cast<Graph&>(*(graph.get()))).nodes)
        {
            if (index.numUniqueKmersOverlappingNode((uint32_t)n.first) < min_unique_kmers_per_node)
            {
                any_below = true;
                break;
            }
        }
        if (any_below)
        {
            continue;
        }

        for (const auto& e : (static_cast<Graph&>(*(graph.get()))).edges)
        {
            if (index.numUniqueKmersOverlappingEdge((uint32_t)e->from(), (uint32_t)e->to()) < min_unique_kmers_per_edge)
            {
                any_below = true;
                break;
            }
        }
        if (any_below)
        {
            continue;
        }

        return k;
    }
    return -1;
}
}
