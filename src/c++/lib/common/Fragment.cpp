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

#include "common/Fragment.hh"

#include "graphalign/GraphAlignmentOperations.hh"

#include "common/Error.hh"

namespace common
{

/**
 * add read and update length estimate
 * @param coordinates coordinates and graph for length calculation
 * @param read read information
 */
void Fragment::addRead(graphtools::GraphCoordinates const& coordinates, Read const& read)
{
    if (fragment_id.empty())
    {
        fragment_id = read.fragment_id();
    }
    else
    {
        assert(read.fragment_id() == fragment_id);
    }
    ++n_reads;

    const bool is_proper_pair = read.is_mapped() && read.is_mate_mapped()
        && ((read.is_reverse_strand() && (!read.is_mate_reverse_strand()))
            || ((!read.is_reverse_strand()) && read.is_mate_reverse_strand()))
        && read.mate_chrom_id() == read.chrom_id();
    if (!is_proper_pair || n_reads > 2)
    {
        bam_fragment_length_ = std::numeric_limits<uint64_t>::max();
    }
    else
    {
        bam_fragment_length_ = abs(read.mate_pos() - read.pos()) + read.bases().size();
    }

    if (read.graph_mapping_status() == common::Read::MAPPED)
    {
        if (read.is_graph_reverse_strand())
        {
            ++n_graph_reverse_reads;
        }
        else
        {
            ++n_graph_forward_reads;
        }

        graphtools::GraphAlignment mapping
            = graphtools::decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), &coordinates.getGraph());
        read_positions_.emplace_back(coordinates.canonicalStartAndEnd(mapping.path()));
        read_lengths_.emplace_back(mapping.queryLength());
        if (read_positions_.size() == 1)
        {
            fragment_length_ = read_lengths_.front();
        }
        else if (read_positions_.size() == 2)
        {
            const uint64_t r1_start = read_positions_.front().first;
            const uint64_t r1_end = read_positions_.front().second;
            const uint64_t r2_start = read_positions_.back().first;
            const uint64_t r2_end = read_positions_.back().second;
            const uint64_t d1 = coordinates.distance(r1_end, r2_start);
            const uint64_t d2 = coordinates.distance(r2_end, r1_start);
            const uint64_t distance = std::min(d1, d2);
            if (distance == std::numeric_limits<uint64_t>::max())
            {
                fragment_length_ = (uint64_t)-1;
            }
            else
            {
                fragment_length_ = read_lengths_.front() + read_lengths_.back() + distance;
            }
        }
        else
        {
            bool has_previous = false;
            uint64_t previous = 0;
            uint64_t length = 0;
            struct ByStart
            {
                bool operator()(std::pair<uint64_t, uint64_t> const& p1, std::pair<uint64_t, uint64_t> const& p2)
                {
                    return p1.first < p2.first;
                }
            };
            read_positions_.sort(ByStart());
            for (auto const& reads : read_positions_)
            {
                if (has_previous)
                {
                    const uint64_t d2p = coordinates.distance(previous, reads.first);
                    if (d2p != std::numeric_limits<uint64_t>::max())
                    {
                        length += d2p;
                    }
                    else
                    {
                        length = std::numeric_limits<uint64_t>::max();
                        break;
                    }
                }
                const uint64_t d2e = coordinates.distance(previous, reads.first);
                if (d2e != std::numeric_limits<uint64_t>::max())
                {
                    length += d2e;
                }
                else
                {
                    length = std::numeric_limits<uint64_t>::max();
                    break;
                }
                previous = reads.first;
                has_previous = true;
            }
            fragment_length_ = length;
        }
    }

    for (auto const& n : read.graph_nodes_supported())
    {
        graph_nodes_supported_.insert(n);
    }
    for (auto const& n : read.graph_edges_supported())
    {
        graph_edges_supported_.insert(n);
    }
    for (auto const& n : read.graph_sequences_supported())
    {
        graph_sequences_supported_.insert(n);
    }
    for (auto const& n : read.graph_sequences_broken())
    {
        graph_sequences_broken_.insert(n);
    }
}

/**
 * Extract fragments from reads
 * @param coordinates graph coordinates structure for length calculation
 * @param reads read buffer
 * @param output_list output list for fragments
 */
void readsToFragments(
    graphtools::GraphCoordinates const& coordinates, common::ReadBuffer const& reads,
    std::list<std::unique_ptr<Fragment>>& output_list)
{
    std::unordered_map<std::string, Fragment*> fragment_map;
    for (auto& read : reads)
    {
        const std::string fragment_id = read->fragment_id();
        auto frag_it = fragment_map.find(fragment_id);
        if (frag_it == fragment_map.end())
        {
            auto fragment = output_list.emplace(output_list.end(), new Fragment());
            frag_it = fragment_map.insert(std::make_pair(fragment_id, fragment->get())).first;
        }
        frag_it->second->addRead(coordinates, *read);
    }
}
}
