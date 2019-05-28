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
