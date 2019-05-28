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
 * \brief CompositeAligner interface
 *
 * \file CompositeAligner.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include "PathAligner.hh"
#include "grm/Filter.hh"
#include "grm/GraphAligner.hh"
#include "grm/KlibAligner.hh"
#include "grm/KmerAligner.hh"
#include "grm/PathAligner.hh"

namespace grm
{

/**
 * Wrapper for a combination of KmerAligner and GraphAligner
 */
class CompositeAligner
{
public:
    CompositeAligner(
        bool pathMatching, bool graphMatching, bool klibMatching, bool kmerMatching,
        unsigned grapAlignmentflags = GraphAligner::AF_ALL);

    virtual ~CompositeAligner();

    CompositeAligner(CompositeAligner&& rhs) noexcept;

    CompositeAligner& operator=(CompositeAligner&& rhs) noexcept = delete;

    void setGraph(graphtools::Graph const* graph, std::list<graphtools::Path> const& paths);
    void alignRead(common::Read& read, ReadFilter filter);

    unsigned attempted() const { return attempted_; }
    unsigned filtered() const { return filtered_; }
    unsigned mappedKlib() const { return mappedKlib_; }
    unsigned mappedPath() const { return mappedPath_; }
    unsigned anchoredPath() const { return anchoredPath_; }
    unsigned mappedKmers() const { return mappedKmers_; }
    unsigned mappedSw() const { return mappedSw_; }

private:
    const bool pathMatching_;
    const bool graphMatching_;
    const bool klibMatching_;
    const bool kmerMatching_;
    const unsigned int grapAlignmentflags_;

    grm::PathAligner pathAligner_;
    grm::GraphAligner graphAligner_;
    grm::KlibAligner klibAligner_;
    grm::KmerAligner<16> kmerAligner_;
    // grm::KmerAligner<32> kmerAligner_;

    unsigned attempted_ = 0;
    unsigned filtered_ = 0;
    unsigned mappedKlib_ = 0;
    unsigned mappedPath_ = 0;
    unsigned anchoredPath_ = 0;
    unsigned mappedKmers_ = 0;
    unsigned mappedSw_ = 0;
#ifdef _DEBUG
    graphtools::Graph const* graph_;
#endif
};
}
