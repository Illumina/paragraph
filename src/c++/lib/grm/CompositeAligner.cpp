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
