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

#include "Read.hh"

#include <list>
#include <memory>
#include <string>
#include <unordered_set>

#include "graphcore/GraphCoordinates.hh"

namespace common
{

class Fragment
{
public:
    Fragment() = default;
    Fragment(Fragment const&) = default;
    Fragment& operator=(Fragment const&) = default;

    Fragment(Fragment&&) = delete;
    Fragment& operator=(Fragment&&) = delete;

    virtual ~Fragment() = default;

    /**
     * Add read and update length estimate
     * @param coordinates coordinates and graph for length calculation
     * @param read read information
     */
    void addRead(graphtools::GraphCoordinates const& coordinates, Read const& read);

    std::string get_fragment_id() const { return fragment_id; }
    int get_n_reads() const { return n_reads; }
    int get_n_graph_forward_reads() const { return n_graph_forward_reads; }
    int get_n_graph_reverse_reads() const { return n_graph_reverse_reads; }
    uint64_t get_bam_fragment_length() const { return bam_fragment_length_; }
    uint64_t get_graph_fragment_length() const { return fragment_length_; }

    std::unordered_set<std::string> const& graph_nodes_supported() const { return graph_nodes_supported_; }
    std::unordered_set<std::string> const& graph_edges_supported() const { return graph_edges_supported_; }
    std::unordered_set<std::string> const& graph_sequences_supported() const { return graph_sequences_supported_; }
    std::unordered_set<std::string> const& graph_sequences_broken() const { return graph_sequences_broken_; }

private:
    std::string fragment_id;
    int n_reads = 0;
    int n_graph_forward_reads = 0;
    int n_graph_reverse_reads = 0;

    std::unordered_set<std::string> graph_nodes_supported_ = {};
    std::unordered_set<std::string> graph_edges_supported_ = {};
    std::unordered_set<std::string> graph_sequences_supported_ = {};
    std::unordered_set<std::string> graph_sequences_broken_ = {};

    uint64_t fragment_length_ = 0;
    std::list<std::pair<uint64_t, uint64_t>> read_positions_ = {};
    std::list<uint64_t> read_lengths_ = {};
    uint64_t bam_fragment_length_ = 0;
};

typedef std::list<std::unique_ptr<common::Fragment>> FragmentList;

/**
 * Extract fragments from reads
 * @param coordinates graph coordinates structure for length calculation
 * @param reads read buffer
 * @param output_list output list for fragments
 */
void readsToFragments(
    graphtools::GraphCoordinates const& coordinates, common::ReadBuffer const& reads,
    std::list<std::unique_ptr<Fragment>>& output_list);
}
