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

#include "graphs/GraphSearch.hh"
#include <graphs/GraphMapping.hh>

namespace graphs
{

size_t singleNodePrefixMatch(WalkableGraph const& wg, uint64_t node, int pos, std::string const& s)
{
    auto const& node_string = wg.node(node)->sequence();
    size_t match_len = 0;
    size_t s_pos = 0;
    while (pos < (int)node_string.size() && s_pos < s.size())
    {
        if (node_string[pos] != s[s_pos])
        {
            break;
        }
        ++s_pos;
        ++pos;
        ++match_len;
    }
    return match_len;
}

/**
 * Match prefix of string to graph, starting at a particular position, extending
 * the match until the first mismatch is encountered. When multiple paths can
 * be matched, this will return the longest match.
 * @param wg a graph
 * @param node start node ID
 * @param pos start position
 * @param s string to match
 * @return graph CIGAR string of the form node[xM]...[zMkS] where k
 *         is the number of unmatched characters in s
 */
std::string prefixMatch(WalkableGraph const& wg, uint64_t node, int pos, std::string const& s)
{
    const size_t this_node_len = wg.node(node)->sequence().size();
    const size_t matched = singleNodePrefixMatch(wg, node, pos, s);
    assert(!s.empty());
    std::string cigar = std::to_string(node) + "[";
    // This handles N's, which graphMapping doesn't like represented as 'M'
    if (matched > 0)
    {
        char op = 'M';
        size_t count = 0;
        for (size_t p = 0; p < matched; ++p)
        {
            if (s[p] == 'N' && op == 'M')
            {
                if (count > 0)
                {
                    cigar += std::to_string(count) + std::string(1, op);
                }
                count = 1;
                op = 'N';
            }
            else if (s[p] != 'N' && op == 'N')
            {
                if (count > 0)
                {
                    cigar += std::to_string(count) + std::string(1, op);
                }
                count = 1;
                op = 'M';
            }
            else
            {
                ++count;
            }
        }
        if (count > 0)
        {
            cigar += std::to_string(count) + std::string(1, op);
        }
    }
    auto succs = wg.succ(node);
    if (matched < s.size() && matched < this_node_len - pos && succs.empty())
    {
        cigar += std::to_string(s.size()) + "S";
    }
    cigar += "]";

    const std::string suffix_left = s.substr(matched);
    // have we matched everything to the end of the node?
    if (!suffix_left.empty() && matched == this_node_len - pos)
    {
        size_t current_best_match_len = 0;
        std::string current_best;
        for (auto const& successor : succs)
        {
            auto succ_cigar = prefixMatch(wg, successor, 0, suffix_left);
            GraphMapping mapping(0, succ_cigar, suffix_left, wg);
            size_t succ_matches = 0;
            for (size_t i = 0; i < mapping.size(); ++i)
            {
                succ_matches += mapping[i].matched();
            }
            if (succ_matches >= current_best_match_len)
            {
                current_best_match_len = succ_matches;
                current_best = succ_cigar;
            }
        }
        cigar += current_best;
    }

    return cigar;
}
}
