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
#include "graphcore/Graph.hh"
#include "graphcore/PathFamily.hh"

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
Json::Value alignAndDisambiguate(const Parameters& parameters, common::ReadBuffer& all_reads);

/**
 * Node and edge filters / return True to indicate a node or edge is supported by a read
 */
typedef std::function<bool(common::Read&, std::string const& node)> ReadSupportsNode;
typedef std::function<bool(common::Read&, std::string const& node1, std::string const& node2)> ReadSupportsEdge;

/**
 * Update sequence labels in read according to nodes the read has traversed
 * @param g graph, passed by reference
 * @param reads list of aligned reads
 * @param nodefilter filter to check if a read supports a particular node
 * @param edgefilter filter to check if a read supports a particular edge
 * @param paths JSON array with information about paths.
 *              Each path must give a "sequence" and a list of nodes.
 */
void disambiguateReads(
    graphtools::Graph* g, std::vector<common::p_Read>& reads, ReadSupportsNode nodefilter = nullptr,
    ReadSupportsEdge edgefilter = nullptr);
}
