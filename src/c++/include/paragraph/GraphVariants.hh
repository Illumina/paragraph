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

#pragma once

#include "Parameters.hh"
#include "common/Read.hh"
#include "graphcore/Graph.hh"
#include "variant/Variant.hh"
#include "json/json.h"

namespace paragraph
{

typedef std::unordered_map<graphtools::NodeId, variant::VariantCandidateList> NodeCandidates;
/**
 * Take apart CIGAR string and collect variant candidates
 * @param g a graph
 * @param read read after alignment
 * @param target vector of candidate lists for nodes
 */
void updateVariantCandidateLists(graphtools::Graph const* g, common::Read const& read, NodeCandidates& target);

/**
 * Extract on-graph variants
 * @param coordinates graph coordinates and graph information
 * @param reads list of reads
 * @param output  output JSON node
 * @param min_reads_for_variant minumum number of reads that must support a variant
 * @param min_frac_for_variant minimum fraction of reads that must support a variant
 * @param paths set of paths to compute coverage over
 * @param write_variants output variants
 * @param write_node_coverage output coverage for nodes
 * @param write_node_coverage output coverage for paths
 */
void getVariants(
    graphtools::GraphCoordinates const& coordinates, common::ReadBuffer const& reads, Json::Value& output,
    int min_reads_for_variant, float min_frac_for_variant, Json::Value const& paths, bool write_variants = false,
    bool write_node_coverage = false, bool write_path_coverage = false);
}
