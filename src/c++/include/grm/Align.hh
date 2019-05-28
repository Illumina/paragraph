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
