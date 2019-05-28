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
