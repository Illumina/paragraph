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
 * \summary Read filtering interface
 *
 * \file ReadFilter.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Read.hh"
#include "graphcore/Graph.hh"

#include <map>
#include <memory>
#include <string>

namespace paragraph
{

class ReadFilter
{
public:
    virtual ~ReadFilter() = default;
    /**
     * Create a read filter
     * @param read the read
     * @return { true if filtered, reason for filtering or "" }
     */
    virtual std::pair<bool, std::string> filterRead(common::Read const& read) = 0;
};

/**
 * Create parameterized read filter
 * @param graph the graph the read was aligned to
 * @param remove_nonuniq remove reads that don't have a unique best alignment
 * @param bad_align_frac fraction of read that must be aligned in order to pass
 * @param kmer_len length of kmers for uniqueness
 * @return a read filter
 */
std::unique_ptr<ReadFilter>
createReadFilter(graphtools::Graph const* graph, bool remove_nonuniq, double bad_align_frac, int32_t kmer_len = 0);
}
