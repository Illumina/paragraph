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
