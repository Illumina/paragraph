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
 * \brief ValiationAligner interface
 *
 * \file ValiationAligner.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <map>
#include <string>

#include "common/Read.hh"
#include "graphs/Graph.hh"

namespace grm
{

/**
 * Uses supplied aligner to perform the alignments. Compares the alignment results
 * with information stored in read name and collects statistics
 */
template <typename AlignerT> class ValidationAligner : private AlignerT
{
public:
    ValidationAligner(AlignerT&& aligner, const graphs::Graph& graph, Json::Value const& paths);

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
