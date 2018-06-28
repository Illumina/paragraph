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
 * \brief Aligner which seeds and extends exact unique matches
 *
 * \file PathAligner.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Read.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include <memory>

namespace grm
{

class PathAligner
{
public:
    explicit PathAligner(int32_t kmer_size = 32);
    virtual ~PathAligner();

    PathAligner(PathAligner&& rhs) noexcept;
    PathAligner& operator=(PathAligner&& rhs) noexcept;

    /**
     * Set the graph to align to
     * @param g a graph
     * @param paths list of paths
     */
    void setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths);

    /**
     * Align a read to the graph and update the graph_* fields.
     *
     * We will align both the forward and the reverse strand and return
     * the alignment which is unique, or with the better score, or default
     * to the forward strand.
     *
     * @param read read structure
     */
    void alignRead(common::Read& read);

    unsigned attempted() const { return attempted_; }
    unsigned anchored() const { return anchored_; }
    unsigned mapped() const { return mapped_; }

private:
    unsigned attempted_ = 0;
    unsigned anchored_ = 0;
    unsigned mapped_ = 0;

    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}