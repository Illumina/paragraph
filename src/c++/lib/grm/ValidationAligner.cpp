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
 * \brief ValidationAligner implementation
 *
 * \file ValidationAligner.cpp
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#include <algorithm>
#include <atomic>
#include <iostream>

#include "common/Error.hh"
#include "grm/CompositeAligner.hh"
#include "grm/ValidationAligner.hh"

namespace grm
{

using graphtools::NodeId;

template <typename AlignerT>
ValidationAligner<AlignerT>::ValidationAligner(
    AlignerT&& aligner, const graphtools::Graph* graph, std::list<graphtools::Path> const& paths)
    : AlignerT(std::move(aligner))
{
    for (const auto& p : paths)
    {
        std::string pathId = p.encode();

        pathNodes_[pathId] = "";
        for (const auto& n : p.nodeIds())
        {
            pathNodes_[pathId] += (pathNodes_[pathId].empty() ? "" : "->") + std::to_string(n);
        }
        LOG()->trace("Validation path: {} : {}", pathId, pathNodes_[pathId]);
    }
}

template <typename AlignerT> void ValidationAligner<AlignerT>::alignRead(common::Read& read, ReadFilter filter)
{
    ++total_;
    AlignerT::alignRead(read, filter);

    if (read.graph_mapping_status() == common::Read::MAPPED)
    {
        ++aligned_;
    }

    if (read.graph_mapping_status() == common::Read::MAPPED)
    {
        const std::string simulatedPathId = getSimulatedPathId(read);
        const std::string cigarNodes = getNodes(read.graph_cigar());
        const std::string simulatedPathNodes = pathNodes_[simulatedPathId];
        bool supportsPath = std::string::npos != simulatedPathNodes.find(cigarNodes);
        mismapped_ += !supportsPath;
        if (!supportsPath)
        {
            LOG()->debug(
                "misplaced:{}:pn'{}' cn'{}'{}", simulatedPathId, pathNodes_[simulatedPathId], cigarNodes,
                read.toJson().toStyledString());
        }
    }
    else if (read.graph_mapping_status() == common::Read::BAD_ALIGN && !read.is_graph_alignment_unique())
    {
        ++repeats_;
        LOG()->debug("repeat:{}", read.toJson().toStyledString());
    }
}

template <typename AlignerT> std::atomic<unsigned> ValidationAligner<AlignerT>::mismapped_(0);
template <typename AlignerT> std::atomic<unsigned> ValidationAligner<AlignerT>::repeats_(0);
template <typename AlignerT> std::atomic<unsigned> ValidationAligner<AlignerT>::aligned_(0);
template <typename AlignerT> std::atomic<unsigned> ValidationAligner<AlignerT>::total_(0);

template <typename AlignerT> std::string ValidationAligner<AlignerT>::getNodes(const std::string& cigar)
{
    std::string cigarNodes;
    bool inCigar = false;
    for (const char c : cigar)
    {
        if ('[' == c)
        {
            inCigar = true;
        }
        else if (']' == c)
        {
            inCigar = false;
        }
        else if (!inCigar)
        {
            if (!cigarNodes.empty())
            {
                cigarNodes += "->";
            }
            cigarNodes += c;
        }
    }
    return cigarNodes;
}

template <typename AlignerT> std::string ValidationAligner<AlignerT>::getSimulatedPathId(common::Read& read)
{
    return read.fragment_id().substr(0, read.fragment_id().find('_'));
}

template class ValidationAligner<CompositeAligner>;
} // namespace
