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

#include "graphs/GaplessAligner.hh"

#include <iostream>
#include <vector>

#include "common/Error.hh"
#include "graphs/GraphPathOperations.hh"

using std::list;
using std::string;
using std::vector;

namespace graphs
{
GraphMapping GaplessAligner::getBestAlignment(const std::string& sequence) const
{
    const list<string> kmers = extractKmersFromAllPositions(sequence, _kmer_len);

    int32_t pos = 0;
    for (const string& kmer : kmers)
    {
        // Initiate alignment from a unique kmer.
        if (_kmer_index.numPaths(kmer) == 1)
        {
            GraphPath kmer_path = _kmer_index.getPaths(kmer).front();
            return getBestAlignmentToShortPath(kmer_path, pos, sequence);
        }
        ++pos;
    }
    return GraphMapping();
}

GraphMapping getBestAlignmentToShortPath(const GraphPath& path, int32_t start_pos, const string& sequence)
{
    const int32_t start_extension = start_pos;
    const int32_t end_extension = sequence.length() - start_pos - path.length();
    const list<GraphPath> full_paths = path.extendBy(start_extension, end_extension);

    GraphMapping best_mapping;
    int32_t max_matches = -1;

    for (const GraphPath& full_path : full_paths)
    {
        GraphMapping mapping = alignWithoutGaps(full_path, sequence);
        if (mapping.numMatches() > max_matches)
        {
            max_matches = mapping.numMatches();
            best_mapping = mapping;
        }
    }

    return best_mapping;
}

GraphMapping alignWithoutGaps(const GraphPath& path, const string& sequence)
{
    vector<string> sequence_pieces = splitByPath(path, sequence);
    vector<Mapping> node_mappings;

    std::shared_ptr<WalkableGraph> wgraph_ptr = path.wgraph_ptr();
    size_t index = 0;
    for (int32_t node_id : path.node_ids())
    {
        const string node_seq = wgraph_ptr->node(node_id)->sequence();
        const string sequence_piece = sequence_pieces[index];
        const int32_t ref_start = index == 0 ? path.start_position() : 0;
        node_mappings.push_back(alignWithoutGaps(sequence_piece, ref_start, node_seq));
        ++index;
    }

    return GraphMapping(path.node_ids(), node_mappings);
}

Mapping alignWithoutGaps(const std::string& query, int32_t ref_start, const std::string& reference)
{
    if (reference.length() < ref_start + query.length())
    {
        error(
            "Gapless alignment requires that read %s is shorter than reference %s.", query.c_str(), reference.c_str());
    }

    if (query.empty() || reference.empty())
    {
        error("Cannot align empty sequences");
    }

    vector<Operation> operations;
    int32_t previous_run_end = 0;
    int32_t run_length = 0;
    char run_operation = '\0';
    for (size_t index = 0; index != query.length(); ++index)
    {
        char cur_operation = 'X';
        if (query[index] == reference[ref_start + index])
        {
            cur_operation = 'M';
        }

        if (cur_operation == run_operation)
        {
            ++run_length;
        }
        else
        {
            if (run_operation != '\0')
            {
                const string query_piece = query.substr(previous_run_end, run_length);
                const string reference_piece = reference.substr(ref_start + previous_run_end, run_length);
                operations.push_back(Operation(run_operation, run_length, query_piece, reference_piece));
            }
            previous_run_end += run_length;
            run_length = 1;
            run_operation = cur_operation;
        }
    }

    const string query_piece = query.substr(previous_run_end, run_length);
    const string reference_piece = reference.substr(ref_start + previous_run_end, run_length);
    operations.push_back(Operation(run_operation, run_length, query_piece, reference_piece));

    return Mapping(ref_start, operations);
}

list<string> extractKmersFromAllPositions(const std::string& sequence, int32_t kmer_len)
{
    list<string> kmers;
    for (size_t pos = 0; pos + kmer_len <= sequence.length(); ++pos)
    {
        kmers.push_back(sequence.substr(pos, kmer_len));
    }
    return kmers;
}
}