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

#include "grm/Align.hh"
#include <graphs/GraphMapping.hh>

#include "grm/GraphAligner.hh"
#include "grm/PathAligner.hh"

#include "common/Error.hh"

using namespace grm;

/**
 * Sequential helper to produce read alignments
 * @param graph graph to align to
 * @param paths list of paths through graph for exact matching; pass NO_PATHS for none
 * @param reads vector of reads that will be updated with graph alignment information
 * @param filter filter function to discard reads if alignment isn't good
 * @param exact_path_matching enable / disable exact matching step (=force Smith Waterman)
 */
static void sequentialAlignReads(
    const graphs::Graph& graph, Json::Value const& paths, std::vector<common::p_Read>& reads, ReadFilter filter,
    bool exact_path_matching)
{
    auto logger = LOG();
    logger->info("[Aligning {} reads]", reads.size());

    grm::GraphAligner graph_aligner;
    graph_aligner.setGraph(graph);

    grm::PathAligner path_aligner;
    if (exact_path_matching && !paths.empty())
    {
        path_aligner.setGraph(graph, paths);
    }

    int c = 0;
    int f = 0;
    std::vector<common::p_Read> filtered_reads;
    int mapped_exactly = 0;
    int mapped_sw = 0;
    for (auto& read : reads)
    {
        if (read->bases().empty())
        {
            continue;
        }
        read->set_graph_mapping_status(reads::UNMAPPED);
        if (exact_path_matching)
        {
            path_aligner.alignRead(*read);
            if (read->graph_mapping_status() == reads::MAPPED)
            {
#ifdef _DEBUG
                // check a valid alignment was produced
                graphs::GraphMapping mapping(read->graph_pos(), read->graph_cigar(), read->bases(), graph);
#endif
                ++mapped_exactly;
            }
        }
        if (read->graph_mapping_status() != reads::MAPPED)
        {
            graph_aligner.alignRead(*read);
            if (read->graph_mapping_status() == reads::MAPPED)
            {
#ifdef _DEBUG
                // check a valid alignment was produced
                graphs::GraphMapping mapping(read->graph_pos(), read->graph_cigar(), read->bases(), graph);
#endif
                ++mapped_sw;
            }
        }
        if (filter != nullptr && filter(*read))
        {
            read->set_graph_mapping_status(reads::MappingStatus::BAD_ALIGN);
            ++f;
        }
        else
        {
            read->set_graph_mapping_status(reads::MappingStatus::MAPPED);
            filtered_reads.emplace_back(std::move(read));
        }
        logger->debug("    [aligned {} reads]", ++c);
    }
    reads.clear();
    for (auto& read : filtered_reads)
    {
        reads.emplace_back(std::move(read));
    }
    logger->info(
        "[Done with alignment step {} total aligned (exact: {} / sw: {}) ; {} were filtered]", c, mapped_exactly,
        mapped_sw, f);
}

/**
 * Wrapper / helper to produce read alignments
 * @param graph graph to align to
 * @param paths list of paths through graph for exact matching; pass NO_PATHS for none
 * @param reads vector of reads that will be updated with graph alignment information
 * @param filter filter function to discard reads if alignment isn't good
 * @param exact_path_matching enable / disable exact matching step (=force Smith Waterman)
 * @param threads number of threads to use for parallel execution
 */
void grm::alignReads(
    const graphs::Graph& graph, Json::Value const& paths, std::vector<common::p_Read>& reads, ReadFilter const& filter,
    bool exact_path_matching, int threads)
{
    if (threads > 1)
    {
        const size_t chunksize = std::min(200ull, (long long unsigned int)reads.size() / threads);
        size_t first_read = 0;
        std::vector<common::p_Read> output_reads;
        std::mutex sync1;
        std::mutex sync2;
        auto worker = [&graph, &paths, &reads, &filter, exact_path_matching, &output_reads, &sync1, &sync2, &first_read,
                       chunksize]() {
            while (first_read < reads.size())
            {
                size_t my_first_read;
                std::vector<common::p_Read> input_reads;
                {
                    std::lock_guard<std::mutex> lock(sync1);
                    my_first_read = first_read;
                    first_read += chunksize;
                    input_reads.resize(chunksize);
                    size_t reads_added = 0;
                    for (size_t i = my_first_read; i < my_first_read + chunksize; ++i)
                    {
                        if (i >= reads.size())
                        {
                            break;
                        }
                        input_reads[i - my_first_read] = std::move(reads[i]);
                        reads_added++;
                    }
                    input_reads.resize(reads_added);
                }
                sequentialAlignReads(graph, paths, input_reads, filter, exact_path_matching);
                {
                    std::lock_guard<std::mutex> lock(sync2);
                    // there probably is a smarter way to do this and keep the same order of the
                    // reads by preallocating output_reads and just moving the unique_ptrs back
                    for (auto& i : input_reads)
                    {
                        output_reads.emplace_back(std::move(i));
                    }
                }
            }
        };

        std::vector<std::thread> workers;
        for (int t = 0; t < threads; ++t)
        {
            workers.emplace_back(worker);
        }
        for (auto& w : workers)
        {
            w.join();
        }
        reads.clear();
        reads.reserve(output_reads.size());
        for (auto& r : output_reads)
        {
            reads.emplace_back(std::move(r));
        }
    }
    else
    {
        sequentialAlignReads(graph, paths, reads, filter, exact_path_matching);
    }
}
