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

/**
 * \brief CompositeAligner implementation
 *
 * \file CompositeAligner.cpp
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#include "grm/CompositeAligner.hh"

#ifdef _DEBUG
#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#endif

using namespace grm;

CompositeAligner::CompositeAligner(
    bool pathMatching, bool graphMatching, bool klibMatching, bool kmerMatching, unsigned grapAlignmentflags)
    : pathMatching_(pathMatching)
    , graphMatching_(graphMatching)
    , klibMatching_(klibMatching)
    , kmerMatching_(kmerMatching)
    , grapAlignmentflags_(grapAlignmentflags)
{
}

CompositeAligner::~CompositeAligner() = default;

CompositeAligner::CompositeAligner(CompositeAligner&& rhs) noexcept = default;

void CompositeAligner::setGraph(graphtools::Graph const* graph, std::list<graphtools::Path> const& paths)
{
    if (pathMatching_)
    {
        pathAligner_.setGraph(graph, paths);
    }

    if (graphMatching_)
    {
        graphAligner_.setGraph(graph);
    }

    if (klibMatching_)
    {
        klibAligner_.setGraph(graph, paths);
    }

    if (kmerMatching_)
    {
        kmerAligner_.setGraph(graph, paths);
    }
#ifdef _DEBUG
    graph_ = graph;
#endif
}

void CompositeAligner::alignRead(common::Read& read, ReadFilter filter)
{
    ++attempted_;

    if (pathMatching_)
    {
        pathAligner_.alignRead(read);
        if (read.graph_mapping_status() == common::Read::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphtools::GraphAlignment mapping
                = graphtools::decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), graph_);
#endif
            ++mappedPath_;
        }
        anchoredPath_ = pathAligner_.anchored();
    }

    // Filter here if filter is set. This allows second-chance alignment with kmer + graph aligner
    if (read.graph_mapping_status() == common::Read::MAPPED && filter && filter(read))
    {
        read.set_graph_mapping_status(common::Read::BAD_ALIGN);
        // increment filtered count if we are not using graph aligner
        filtered_ += !kmerMatching_ && !klibMatching_ && !graphMatching_;
    }

    if (read.graph_mapping_status() != common::Read::MAPPED && kmerMatching_)
    {
        kmerAligner_.alignRead(read);
        if (read.graph_mapping_status() == common::Read::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphtools::GraphAlignment mapping
                = graphtools::decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), graph_);
#endif
            if (filter && filter(read))
            {
                read.set_graph_mapping_status(common::Read::BAD_ALIGN);
                // increment filtered count if we are not trying to align further
                filtered_ += !klibMatching_ && !graphMatching_;
            }
            else
            {
                ++mappedKmers_;
            }
        }
    }

    if (read.graph_mapping_status() != common::Read::MAPPED && klibMatching_)
    {
        klibAligner_.alignRead(read);
        // Filter here if filter is set. This allows second-chance alignment with graph aligner
        if (read.graph_mapping_status() == common::Read::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphtools::GraphAlignment mapping
                = graphtools::decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), graph_);
#endif
            if (filter && filter(read))
            {
                read.set_graph_mapping_status(common::Read::BAD_ALIGN);
                // increment filtered count if we are not using graph aligner
                filtered_ += !graphMatching_;
            }
            else
            {
                ++mappedKlib_;
            }
        }
    }

    if (read.graph_mapping_status() != common::Read::MAPPED && graphMatching_)
    {
        graphAligner_.alignRead(read);
        // graph aligner always produces a mapping, It just does not set the status for some reason
        read.set_graph_mapping_status(common::Read::MAPPED);

        if (read.graph_mapping_status() == common::Read::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphtools::GraphAlignment mapping
                = graphtools::decodeGraphAlignment(read.graph_pos(), read.graph_cigar(), graph_);
#endif
            if (filter && filter(read))
            {
                read.set_graph_mapping_status(common::Read::BAD_ALIGN);
                ++filtered_;
            }
            else
            {
                ++mappedSw_;
            }
        }
    }
}
