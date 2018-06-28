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

#include <memory>
#include <string>
#include <vector>

#include "common/Read.hh"

#include "graphcore/Graph.hh"

namespace grm
{

/**
 * Wrapper for GSSW graph alignment code
 */
class GraphAligner
{
public:
    GraphAligner();

    virtual ~GraphAligner();

    GraphAligner(GraphAligner&& rhs) noexcept;

    GraphAligner& operator=(GraphAligner&& rhs) noexcept;

    /**
     * Set the graph to align to
     * @param g a graph
     */
    void setGraph(graphtools::Graph const* g);

    /**
     * Smith-Waterman align a string to the graph and return a cigar string and
     * mapping score which can be either 0 or 60. Score of 60 means
     * that all top-scoring alignments end as the same node; otherwise
     * the score is 0.
     *
     * @param read read sequence
     * @param mapq mapping score
     * @param position start position in first node
     * @param score alignment score
     * @return cigar string
     */
    std::string align(const std::string& read, int& mapq, int& position, int& score) const;

    /** alignment flags */
    static const unsigned int AF_CIGAR = 0x01;
    static const unsigned int AF_BOTH_STRANDS = 0x02;
    static const unsigned int AF_REVERSE_GRAPH = 0x04;
    static const unsigned int AF_ALL = (unsigned int)-1;

    /**
     * Align a read to the graph and update the graph_* fields.
     *
     * We will align both the forward and the reverse strand and return
     * the alignment which is unique, or with the better score, or default
     * to the forward strand.
     *
     * @param read read structure
     * @param alignment_flags flags for alignment (see above)
     */
    void alignRead(common::Read& read, unsigned int alignment_flags = AF_ALL) const;

private:
    struct GraphAlignerImpl;
    std::unique_ptr<GraphAlignerImpl> _impl;
};
}