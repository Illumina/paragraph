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

#include <boost/range.hpp>

#include "common/Error.hh"
#include "common/Threads.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "grm/Align.hh"
#include "grm/CompositeAligner.hh"
#include "grm/ValidationAligner.hh"

using namespace grm;

using common::Read;

void logAlignerStats(const CompositeAligner& aligner)
{
    LOG()->info(
        "[Done with alignment step {} total aligned (path: {} [{} anchored] kmers: {} / ksw: {} / gssw: {}) ; {} were "
        "filtered]",
        aligner.attempted(), aligner.mappedPath(), aligner.anchoredPath(), aligner.mappedKlib(), aligner.mappedKmers(),
        aligner.mappedSw(), aligner.filtered());
}

template <typename AlignerT> void logAlignerStats(const ValidationAligner<AlignerT>& aligner)
{
    logAlignerStats(aligner.base());

    LOG()->info("[VALIDATION]\tMAPQ\tEmpMAPQ\tWrong\tTotal");
    LOG()->info("[VALIDATION]\tunalgnd\t0\t0\t{}", (aligner.total() - aligner.aligned() - aligner.repeats()));
    LOG()->info("[VALIDATION]\trepeat\t0\t0\t{}", aligner.repeats());
    LOG()->info(
        "[VALIDATION]\t60\t{}\t{}\t{}",
        (!aligner.mismapped()
             ? 60
             : (aligner.aligned()) ? -10 * log10(double(aligner.mismapped()) / double(aligner.aligned())) : 0),
        aligner.mismapped(), aligner.aligned());
}

/**
 * Sequential helper to produce read alignments
 * @param graph graph to align to
 * @param paths list of paths through graph for exact matching; pass NO_PATHS for none
 * @param reads vector of reads that will be updated with graph alignment information
 * @param filter filter function to discard reads if alignment isn't good
 */
template <typename IteratorT, typename AlignerT>
static void sequentialAlignReads(
    const IteratorT begin, IteratorT end, const graphtools::Graph* graph, std::list<graphtools::Path> const& paths,
    ReadFilter filter, std::vector<common::p_Read>& filtered_reads, AlignerT& aligner)
{
    auto logger = LOG();
    logger->info("[Aligning {} reads]", std::distance(begin, end));

    for (auto& read : boost::make_iterator_range(begin, end))
    {
        if (read->bases().empty())
        {
            continue;
        }
        read->set_graph_mapping_status(Read::UNMAPPED);
        aligner.alignRead(*read, filter);

        if (Read::MAPPED == read->graph_mapping_status())
        {
            filtered_reads.emplace_back(std::move(read));
        }
    }

    logAlignerStats(aligner);
}

template <typename IteratorT>
static void sequentialAlignReads(
    const IteratorT begin, IteratorT end, const graphtools::Graph* graph, std::list<graphtools::Path> const& paths,
    ReadFilter filter, bool path_sequence_matching, bool graph_sequence_matching, bool klib_sequence_matching,
    bool kmer_sequence_matching, bool validate_alignments, std::vector<common::p_Read>& filtered_reads)
{
    if (validate_alignments)
    {
        grm::ValidationAligner<grm::CompositeAligner> aligner(
            grm::CompositeAligner(
                path_sequence_matching, graph_sequence_matching, klib_sequence_matching, kmer_sequence_matching),
            graph, paths);
        aligner.setGraph(graph, paths);
        sequentialAlignReads(begin, end, graph, paths, filter, filtered_reads, aligner);
    }
    else
    {
        grm::CompositeAligner aligner(
            path_sequence_matching, graph_sequence_matching, klib_sequence_matching, kmer_sequence_matching);
        aligner.setGraph(graph, paths);
        sequentialAlignReads(begin, end, graph, paths, filter, filtered_reads, aligner);
    }
}

void grm::alignReads(
    const graphtools::Graph* graph, std::list<graphtools::Path> const& paths, std::vector<common::p_Read>& reads,
    ReadFilter const& filter, bool path_sequence_matching, bool graph_sequence_matching, bool klib_sequence_matching,
    bool kmer_sequence_matching, bool validate_alignments, uint32_t threads)
{
    auto next = reads.begin();
    const std::size_t step = std::max((reads.size() + threads - 1) / threads, std::size_t(1));
    std::mutex m;
    std::vector<common::p_Read> allFilteredReads;
    bool terminate = false;
    common::CPU_THREADS(threads).execute(
        [&]() {
            std::unique_lock<std::mutex> lock(m);
            while (reads.end() != next)
            {
                const std::size_t ourStep
                    = std::min<std::size_t>(step, static_cast<const size_t>(std::distance(next, reads.end())));
                if (ourStep)
                {
                    auto begin = next;
                    next += ourStep;
                    auto end = next;
                    std::vector<common::p_Read> filteredReads;
                    ASYNC_BLOCK_WITH_CLEANUP([&](bool failure) { terminate |= failure; })
                    {
                        if (terminate)
                        {
                            LOG()->warn("terminating");
                            break;
                        }
                        common::unlock_guard<std::unique_lock<std::mutex>> unlock(lock);
                        sequentialAlignReads(
                            begin, end, graph, paths, filter, path_sequence_matching, graph_sequence_matching,
                            klib_sequence_matching, kmer_sequence_matching, validate_alignments, filteredReads);
                    }
                    std::move(filteredReads.begin(), filteredReads.end(), std::back_inserter(allFilteredReads));
                }
            }
        },
        std::max(reads.size() / step, std::size_t(1)));

    reads.swap(allFilteredReads);
}
