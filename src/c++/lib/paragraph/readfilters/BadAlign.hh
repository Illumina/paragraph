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
 * \summary Read filter based on alignment quality
 *
 * \file BadAlign.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#include <algorithm>

#include "graphs/GraphMapping.hh"
#include "graphs/GraphMappingOperations.hh"
#include "paragraph/ReadFilter.hh"

#include "common/Error.hh"

namespace paragraph
{
namespace readfilters
{
    class BadAlign : public ReadFilter
    {
    public:
        explicit BadAlign(graphs::WalkableGraph const* graph, double bad_align_frac)
            : graph_(graph)
            , bad_align_frac_(bad_align_frac)
        {
            assert(graph_ != nullptr);
        }

        std::pair<bool, std::string> filterRead(common::Read const& r) override
        {
            const graphs::GraphMapping mapping
                = graphs::decodeFromString(r.graph_pos(), r.graph_cigar(), r.bases(), *graph_);
            const auto query_aligned = mapping.querySpan() - mapping.queryClipped();
            const bool is_bad = query_aligned < (int)round(bad_align_frac_ * (mapping.querySpan()));
            return { is_bad, is_bad ? "bad_align" : "" };
        }

    private:
        graphs::WalkableGraph const* graph_;
        double bad_align_frac_;
    };
}
}