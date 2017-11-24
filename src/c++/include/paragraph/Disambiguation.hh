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
 * \brief Wrapper for alignment and disambiguation of reads
 *
 * \file Disambiguation.hh
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "Parameters.hh"
#include "common/Read.hh"
#include "variant/Variant.hh"
#include "json/json.h"

#include <graphs/WalkableGraph.hh>
#include <unordered_map>
#include <vector>

namespace paragraph
{

/**
 * Align reads from single BAM file to graph and disambiguate reads
 * to produce counts.
 *
 * @param parameters Graph alignment parameters
 * @param all_reads pass read bufer with all reads to be aligned and disambiguated
 * @return results as JSON value
 */
Json::Value alignAndDisambiguate(Parameters& parameters, common::ReadBuffer& all_reads);

/**
 * Internal functions (exposed for testing)
 */

/**
 * Take apart CIGAR string and collect variant candidates
 * @param read read after alignment
 * @param target vector of candidate lists for nodes
 */
void updateVariantCandidateLists(
    graphs::Graph& g, common::Read const& read, std::unordered_map<uint64_t, variant::VariantCandidateList>& target);

/**
 * Node and edge filters / return True to indicate a node or edge is supported by a read
 */
typedef std::function<bool(common::Read&, std::string const& node)> ReadSupportsNode;
typedef std::function<bool(common::Read&, std::string const& node1, std::string const& node2)> ReadSupportsEdge;

/**
 * Update sequence labels in read according to nodes the read has traversed
 * @param g graph structure
 * @param reads list of aligned reads
 * @param nodefilter filter to check if a read supports a particular node
 * @param edgefilter filter to check if a read supports a particular edge
 * @param paths JSON array with information about paths.
 *              Each path must give a "sequence" and a list of nodes.
 */
void disambiguateReads(
    graphs::WalkableGraph& g, std::vector<common::p_Read>& reads, ReadSupportsNode nodefilter = nullptr,
    ReadSupportsEdge edgefilter = nullptr, Json::Value const& paths = Json::Value(Json::objectValue));
}
