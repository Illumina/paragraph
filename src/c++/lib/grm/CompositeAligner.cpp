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
#include "graphs/GraphMapping.hh"
#include "graphs/GraphMappingOperations.hh"
using namespace grm;

CompositeAligner::CompositeAligner(
    bool exactPathMatching, bool graphMatching, bool kmerMatching, unsigned grapAlignmentflags)
    : exactPathMatching_(exactPathMatching)
    , graphMatching_(graphMatching)
    , kmerMatching_(kmerMatching)
    , grapAlignmentflags_(grapAlignmentflags)
{
}

CompositeAligner::~CompositeAligner() = default;

CompositeAligner::CompositeAligner(CompositeAligner&& rhs) noexcept
    : exactPathMatching_(rhs.exactPathMatching_)
    , graphMatching_(rhs.graphMatching_)
    , kmerMatching_(rhs.kmerMatching_)
    , grapAlignmentflags_(rhs.grapAlignmentflags_)
    , graphAligner_(std::move(rhs.graphAligner_))
    , kmerAligner_(std::move(rhs.kmerAligner_))
    , pathAligner_(std::move(rhs.pathAligner_))
{
}

void CompositeAligner::setGraph(graphs::Graph const& graph, Json::Value const& paths)
{
    if (exactPathMatching_)
    {
        pathAligner_.setGraph(graph, paths);
    }
    if (graphMatching_)
    {
        graphAligner_.setGraph(graph);
    }
    if (kmerMatching_)
    {
        kmerAligner_.setGraph(graph, paths);
    }
#ifdef _DEBUG
    graph_ = std::unique_ptr<graphs::WalkableGraph>(new graphs::WalkableGraph(graph));
#endif
}

void CompositeAligner::alignRead(common::Read& read, ReadFilter filter)
{
    ++attempted_;
    if (exactPathMatching_)
    {
        pathAligner_.alignRead(read);
        if (read.graph_mapping_status() == reads::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphs::GraphMapping mapping
                = graphs::decodeFromString(read.graph_pos(), read.graph_cigar(), read.bases(), *graph_);
#endif
            ++mappedExactly_;
        }
    }

    if (read.graph_mapping_status() != reads::MAPPED && kmerMatching_)
    {
        kmerAligner_.alignRead(read);
        if (read.graph_mapping_status() == reads::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphs::GraphMapping mapping
                = graphs::decodeFromString(read.graph_pos(), read.graph_cigar(), read.bases(), *graph_);
#endif
            ++mappedKmers_;
        }
    }

    // Filter here if filter is set. This allows second-chance alignment with graph aligner
    if (read.graph_mapping_status() == reads::MAPPED && filter && filter(read))
    {
        read.set_graph_mapping_status(reads::MappingStatus::BAD_ALIGN);
        // increment filtered count if we are not using graph aligner
        filtered_ += !graphMatching_;
    }

    if (read.graph_mapping_status() != reads::MAPPED && graphMatching_)
    {
        graphAligner_.alignRead(read);
        // graph aligner always produces a mapping, It just does not set the status for some reason
        read.set_graph_mapping_status(reads::MappingStatus::MAPPED);

        if (read.graph_mapping_status() == reads::MAPPED)
        {
#ifdef _DEBUG
            // check a valid alignment was produced
            graphs::GraphMapping mapping
                = graphs::decodeFromString(read.graph_pos(), read.graph_cigar(), read.bases(), *graph_);
#endif
            ++mappedSw_;
        }

        if (filter && filter(read))
        {
            read.set_graph_mapping_status(reads::MappingStatus::BAD_ALIGN);
            ++filtered_;
        }
    }
}
