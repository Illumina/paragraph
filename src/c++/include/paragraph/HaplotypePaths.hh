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
 * \summary Find haplotypes from mapped reads
 *
 * \file HaplotypePaths.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Read.hh"
#include "graphcore/Path.hh"
#include "graphcore/PathFamily.hh"
#include <memory>

namespace paragraph
{

typedef std::pair<graphtools::PathFamily, int> PhasingFamily;

/**
 * Summarize phasing evidence from read/pair alignments
 * @param graph the graph
 * @param reads aligned reads
 * @return Number of fragments directly phasing together each path family
 */
std::vector<PhasingFamily> getPhasingFamilies(graphtools::Graph* graph, common::ReadBuffer const& reads);

/**
 * Add paths based on read-supported haplotypes to graph and output JSON
 * @param reads
 * @param[out] graph
 * @param[out] paths
 * @param[out] output
 */
void addHaplotypePaths(
    common::ReadBuffer const& reads, graphtools::Graph& graph, Json::Value& paths, Json::Value& output);
}
