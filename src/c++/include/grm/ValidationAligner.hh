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
 * \brief ValiationAligner interface
 *
 * \file ValiationAligner.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <list>
#include <map>
#include <string>

#include "common/Read.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace grm
{

/**
 * Uses supplied aligner to perform the alignments. Compares the alignment results
 * with information stored in read name and collects statistics
 */
template <typename AlignerT> class ValidationAligner : private AlignerT
{
public:
    ValidationAligner(AlignerT&& aligner, const graphtools::Graph* graph, std::list<graphtools::Path> const& paths);

    virtual ~ValidationAligner() { ; }

    using AlignerT::setGraph;

    void alignRead(common::Read& read, ReadFilter filter);
    const AlignerT& base() const { return *this; }
    static unsigned mismapped() { return mismapped_; }
    static unsigned repeats() { return repeats_; }
    static unsigned aligned() { return aligned_; }
    static unsigned total() { return total_; }

private:
    std::unordered_map<std::string, std::string> pathNodes_;
    static std::atomic<unsigned> mismapped_;
    static std::atomic<unsigned> repeats_;
    static std::atomic<unsigned> aligned_;
    static std::atomic<unsigned> total_;

    static std::string getNodes(const std::string& cigar);
    static std::string getSimulatedPathId(common::Read& read);
};

extern template class ValidationAligner<CompositeAligner>;
}
