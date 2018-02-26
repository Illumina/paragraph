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

#include "graphs/GraphPathOperations.hh"

#include "common/Error.hh"

using std::list;
using std::string;
using std::vector;

namespace graphs
{
vector<string> splitByPath(const GraphPath& path, const std::string& sequence)
{
    if (path.length() != sequence.length())
    {
        error("Split operation requires that %s and %s have same length", path.encode().c_str(), sequence.c_str());
    }

    vector<string> split_seq;
    const vector<int32_t> path_node_ids = path.node_ids();

    size_t cur_position = 0;
    for (int32_t node_id : path_node_ids)
    {
        const size_t length_on_node = path.lengthOnNode(node_id);
        split_seq.push_back(sequence.substr(cur_position, length_on_node));
        cur_position += length_on_node;
    }
    return split_seq;
}

/**
 * Return true if two paths overlap either prefix - suffix or suffix-prefix
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths overlap
 */
bool checkPathPrefixSuffixOverlap(GraphPath const& p1, GraphPath const& p2)
{
    // technically we'd want to check that the two graphs match also

    if (p1.node_ids().empty() || p2.node_ids().empty())
    {
        return false;
    }
    if (p1.node_ids().back() < p2.node_ids().front() || // p1 ends before p2
        p1.node_ids().front() > p2.node_ids().back()) // p1 starts after p2
    {
        return false;
    }

    auto p1_it = p1.node_ids().begin();
    auto p2_it = p2.node_ids().begin();

    int shared_nodes = 0;
    while (p1_it != p1.node_ids().end() && p2_it != p2.node_ids().end())
    {
        if (*p1_it < *p2_it)
        {
            if (p2_it != p2.node_ids().begin())
            {
                // paths diverged
                return false;
            }
            // --> ignore non-matching prefix of p1 until paths meet
            ++p1_it;
        }
        else if (*p1_it > *p2_it)
        {
            if (p1_it != p1.node_ids().begin())
            {
                // paths diverged
                return false;
            }
            // --> ignore non-matching prefix of p2 until paths meet
            ++p2_it;
        }
        else
        { // *p1_it == *p2_it
            // paths have met. They must now match until one of them ends
            ++shared_nodes;
            ++p1_it;
            ++p2_it;
        }
    }

    if (shared_nodes == 0)
    {
        return false;
    }

    // if we only share one node, the paths may not overlap on that node
    if (shared_nodes == 1)
    {
        if (p1_it == p1.node_ids().end() && p2_it == p2.node_ids().end())
        {
            if (p1.num_nodes() > 1 && p2.num_nodes() > 1)
            {
                // if they both have > 1 nodes, they should also share more than one of them;
                // otherwise they cannot both end here
                assert(false);
            }
            else if (p1.num_nodes() == 1 && p2.num_nodes() > 1)
            {
                // p1 starts here, p2 ends here
                if (p2.end_position() < p1.start_position())
                {
                    return false;
                }
            }
            else if (p1.num_nodes() > 1 && p2.num_nodes() == 1)
            {
                // p2 starts here, p1 ends here
                if (p1.end_position() < p2.start_position())
                {
                    return false;
                }
            }
            else if (p1.num_nodes() == 1 && p2.num_nodes() == 1)
            {
                // both paths on same node, check if they overlap there
                return p1.end_position() >= p2.start_position() && p2.end_position() >= p1.start_position();
            }
        }
        else if (p1_it != p1.node_ids().end() && p2_it == p2.node_ids().end())
        {
            // p2 starts+ends on this node p1 starts -- check that p1 starts before p2 ends
            if (p2.end_position() < p1.start_position())
            {
                return false;
            }
        }
        else if (p1_it == p1.node_ids().end() && p2_it != p2.node_ids().end())
        {
            // p1 starts+ends on this node p2 starts -- check that p2 starts before p1 ends
            if (p1.end_position() < p2.start_position())
            {
                return false;
            }
        }
        else
        {
            // this shouldn't happen. we iterate until one of them reaches end() above
            assert(false);
        }
    }

    return true;
}

/**
 * Paths can be merged if they overlap prefix-suffix / suffix-prefix.
 *
 * @param p1 first path
 * @param p2 second path
 * @return merged path
 */
GraphPath mergePaths(GraphPath const& p1, GraphPath const& p2)
{
    assert(checkPathPrefixSuffixOverlap(p1, p2));

    int32_t start = -1;
    int32_t end = -1;
    std::vector<int32_t> nodes;
    auto p1_it = p1.node_ids().begin();
    auto p2_it = p2.node_ids().begin();
    while ((p1_it != p1.node_ids().end()) && (p2_it != p2.node_ids().end()))
    {
        if (*p1_it < *p2_it)
        {
            if (start < 0)
            {
                start = p1.start_position();
            }
            nodes.push_back(*p1_it);
            ++p1_it;
        }
        else if (*p1_it > *p2_it)
        {
            if (start < 0)
            {
                start = p2.start_position();
            }
            nodes.push_back(*p2_it);
            ++p2_it;
        }
        else
        { // *p1_it == *p2_it
            if (start < 0)
            {
                start = std::min(p1.start_position(), p2.start_position());
            }
            nodes.push_back(*p1_it);
            ++p1_it;
            ++p2_it;
        }
    }
    if (p1_it == p1.node_ids().end() && p2_it == p2.node_ids().end())
    {
        end = std::max(p1.end_position(), p2.end_position());
    }
    else if (p1_it != p1.node_ids().end() && p2_it == p2.node_ids().end())
    {
        nodes.insert(nodes.end(), p1_it, p1.node_ids().end());
        end = p1.end_position();
    }
    else if (p1_it == p1.node_ids().end() && p2_it != p2.node_ids().end())
    {
        nodes.insert(nodes.end(), p2_it, p2.node_ids().end());
        end = p2.end_position();
    }
    assert(start >= 0 && end >= 0);
    return GraphPath(p1.wgraph_ptr(), start, nodes, end);
}

/**
 * Merge a set of paths
 *
 * This will merge paths until none of the resulting paths overlap
 *
 * @param paths a list of paths
 */
void greedyMerge(std::list<GraphPath>& paths)
{
    bool has_merged = true;
    while (has_merged && paths.size() > 1)
    {
        auto path_a = paths.begin();
        has_merged = false;
        while (path_a != paths.end())
        {
            auto path_b = std::next(path_a);
            while (path_b != paths.end())
            {
                if (checkPathPrefixSuffixOverlap(*path_a, *path_b))
                {
                    const GraphPath merged_a_b{ mergePaths(*path_a, *path_b) };
                    paths.erase(path_a);
                    paths.erase(path_b);
                    paths.push_back(merged_a_b);
                    has_merged = true;
                    break;
                }
                ++path_b;
            }
            if (has_merged)
            {
                break;
            }
            ++path_a;
        }
    }
}

/**
 * Merge a set of paths
 *
 * This will merge paths exhaustively, each path is merged with all
 * paths it overlaps until we cannot merge anymore
 *
 * @param paths a list of paths
 */
void exhaustiveMerge(std::list<GraphPath>& paths)
{
    bool has_merged = true;
    while (has_merged && paths.size() > 1)
    {
        auto path_a = paths.begin();
        has_merged = false;

        list<GraphPath> new_paths;
        while (path_a != paths.end())
        {
            auto path_b = paths.begin();
            while (path_b != paths.end())
            {
                if (path_a == path_b)
                {
                    ++path_b;
                    continue;
                }
                if (checkPathPrefixSuffixOverlap(*path_a, *path_b))
                {
                    const GraphPath merged_a_b{ mergePaths(*path_a, *path_b) };
                    const bool a_contained_in_b = merged_a_b.encode() == path_b->encode();
                    const bool b_contained_in_a = merged_a_b.encode() == path_a->encode();
                    const bool a_eq_b = a_contained_in_b && b_contained_in_a;

                    if (a_eq_b)
                    {
                        // keep only one of them
                        new_paths.push_back(*path_b);
                    }
                    else if (a_contained_in_b || b_contained_in_a)
                    {
                        new_paths.push_back(merged_a_b);
                    }
                    else
                    {
                        new_paths.push_back(merged_a_b);
                        new_paths.push_back(*path_a);
                        new_paths.push_back(*path_b);
                    }
                    ++path_b;
                    has_merged = true;
                }
                else
                {
                    new_paths.push_back(*path_b);
                    ++path_b;
                }
            }
            if (has_merged)
            {
                break;
            }
            new_paths.push_back(*path_a);
            ++path_a;
        }
        if (has_merged)
        {
            paths = new_paths;
        }
    }
}
}