// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Paragraph
// Copyright (c) 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//		http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
// See the License for the specific language governing permissions and limitations
//
//

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