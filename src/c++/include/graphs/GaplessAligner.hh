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

#pragma once

#include <list>
#include <string>

#include "graphs/GraphMapping.hh"
#include "graphs/GraphPath.hh"
#include "graphs/KmerIndex.hh"

namespace graphs
{
/**
 * @brief Gapless graph aligner
 *
 * Enables gapless alignment of any sequence to a graph.
 */
class GaplessAligner
{
public:
    GaplessAligner(std::shared_ptr<WalkableGraph> wgraph_ptr, int32_t kmer_len)
        : _kmer_len(kmer_len)
        , _kmer_index(wgraph_ptr, kmer_len)
    {
    }
    GraphMapping getBestAlignment(const std::string& sequence) const;

private:
    int32_t _kmer_len;
    KmerIndex _kmer_index;
};

/**
 * @brief Computes a top-scoring gapless alignment of a sequence to the graph that goes through the path starting at the
 * given position on the sequence
 *
 * @param path: Any path shorter than the sequence
 * @param start_pos: Position on the sequence corrsponding to the start of the path
 * @param sequence: Any sequence
 * @return Best gapless alignment with the above properety
 */
GraphMapping getBestAlignmentToShortPath(const GraphPath& path, int32_t start_pos, const std::string& sequence);

/**
 * @brief Aligns a sequence to a path of the same length
 *
 * @param path: Any graph path
 * @param sequence: Any sequence that has the same length as the path
 * @return Result of the alignment
 */
GraphMapping alignWithoutGaps(const GraphPath& path, const std::string& sequence);

/**
 * @brief Aligns query sequence to the reference sequence starting at the given position
 *
 * @param query: Any sequence
 * @param ref_start: Position of the start of the alignment on the reference
 * @param reference: Any sequence
 * @return Result of the alignment
 */
Mapping alignWithoutGaps(const std::string& query, int32_t ref_start, const std::string& reference);

/**
 * @brief Extracts kmers starting at each position
 *
 * @param sequence: Any sequence
 * @param kmer_len: Kmer length
 * @return List of kmers indexed by start position in the original sequence
 */
std::list<std::string> extractKmersFromAllPositions(const std::string& sequence, int32_t kmer_len);
}