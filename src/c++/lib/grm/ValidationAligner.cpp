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

template <typename AlignerT>
ValidationAligner<AlignerT>::ValidationAligner(AlignerT&& aligner, const graphs::Graph& graph, Json::Value const& paths)
    : AlignerT(std::move(aligner))
{
    std::unordered_map<std::string, uint64_t> nodeNameToId;

    for (const auto& n : graph.nodes)
    {
        assert(nodeNameToId.count(n.second->name()) == 0);
        nodeNameToId[n.second->name()] = n.first;
    }

    for (const auto& p : paths)
    {
        std::string pathId = p["path_id"].asString();
        std::replace(pathId.begin(), pathId.end(), '|', '-');

        pathNodes_[pathId] = "";
        for (const auto& n : p["nodes"])
        {
            const auto nodeName = n.asString();
            pathNodes_[pathId] += (pathNodes_[pathId].empty() ? "" : "->") + std::to_string(nodeNameToId[nodeName]);
        }
        LOG()->debug("Validation path: {} : {}", pathId, pathNodes_[pathId]);
    }
}

template <typename AlignerT> void ValidationAligner<AlignerT>::alignRead(common::Read& read)
{
    ++total_;
    AlignerT::alignRead(read);

    if (read.graph_mapping_status() == reads::MAPPED)
    {
        ++aligned_;
    }

    if (read.graph_mapping_status() == reads::MAPPED)
    {
        const std::string simulatedPathId = getSimulatedPathId(read);
        const std::string cigarNodes = getNodes(read.graph_cigar());
        const std::string simulatedPathNodes = pathNodes_[simulatedPathId];
        bool supportsPath = std::string::npos != simulatedPathNodes.find(cigarNodes);
        mismapped_ += !supportsPath;
        if (!supportsPath)
        {
            LOG()->debug(
                "misplaced:{}:{} {}{}", simulatedPathId, pathNodes_[simulatedPathId], cigarNodes,
                read.asJson().toStyledString());
        }
    }
    else if (read.graph_mapping_status() == reads::BAD_ALIGN && !read.is_graph_alignment_unique())
    {
        ++repeats_;
        LOG()->debug("repeat:{}", read.asJson().toStyledString());
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
