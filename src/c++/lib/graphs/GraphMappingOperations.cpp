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

#include "graphs/GraphMappingOperations.hh"

#include "common/Error.hh"

using std::string;
using std::vector;

namespace graphs
{

GraphMapping
decodeFromString(int32_t first_node_start, const string& graph_cigar, const string& query, const Graph& graph)
{
    vector<int32_t> node_ids;
    vector<Mapping> node_mappings;
    int32_t query_pos = 0;
    string node_cigar;
    for (size_t index = 0; index != graph_cigar.length(); ++index)
    {
        node_cigar += graph_cigar[index];
        if (node_cigar.back() == ']')
        {
            string query_piece = query.substr((size_t)query_pos);
            int32_t ref_pos = node_mappings.empty() ? first_node_start : 0;

            string cigar;
            int32_t node_id;
            splitNodeCigar(node_cigar, cigar, node_id);
            node_ids.push_back(node_id);
            const string& node_seq = graph.nodes.at(node_id)->sequence();
            Mapping node_mapping(ref_pos, cigar, query_piece, node_seq);
            node_mappings.push_back(node_mapping);
            query_pos += node_mapping.querySpan();
            node_cigar.clear();
        }
    }
    return GraphMapping(node_ids, node_mappings);
}

void splitNodeCigar(const string& node_cigar, string& cigar, int32_t& node_id)
{
    node_id = -1;
    string nodeid_encoding;
    for (size_t index = 0; index != node_cigar.length(); ++index)
    {
        if (node_cigar[index] == '[')
        {
            node_id = std::stoull(nodeid_encoding);
            cigar = node_cigar.substr(index + 1);
            cigar.pop_back();
            return;
        }
        if (isdigit(node_cigar[index]) == 0)
        {
            error("Error: %s is a malformed node CIGAR", node_cigar.c_str());
        }
        nodeid_encoding += node_cigar[index];
    }
}
}