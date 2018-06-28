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

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "common/ReadExtraction.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"
#include "grm/Filter.hh"
#include "json/json.h"

namespace grm
{
static const Json::Value NO_PATHS = Json::objectValue;

/**
 * Wrapper / helper to produce read alignments
 * @param graph graph to align to
 * @param paths list of paths through graph for exact matching; pass NO_PATHS for none
 * @param reads vector of reads that will be updated with graph alignment information
 * @param filter filter function to discard reads if alignment isn't good
 * @param graph_sequence_matching enable smith waterman graph sequence matching
 * @param kmer_sequence_matching enable kmer sequence matching
 * @param validate_alignments enable validation using read ids
 * @param threads number of threads to use for parallel execution
 */
void alignReads(
    const graphtools::Graph* graph, std::list<graphtools::Path> const& paths, std::vector<common::p_Read>& reads,
    ReadFilter const& filter, bool path_sequence_matching, bool graph_sequence_matching, bool klib_sequence_matching,
    bool kmer_sequence_matching, bool validate_alignments, uint32_t threads = 1);
}
