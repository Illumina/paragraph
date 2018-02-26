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
 * Data structure for one breakpoint info across samples
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/SampleInfo.hh"
#include "graphs/WalkableGraph.hh"
#include "json/json.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace genotyping
{

/**
 * Counts corresponding to breakpoint alleles
 */
typedef std::vector<int32_t> AlleleCounts;

/**
 * Class to collect single-breakpoint statistics / edge counts
 *
 * Future versions may incorporate additional alignment information (mapping qualities, strands, ...),
 * currently, this only stores the read count per allele and per edge.
 *
 */
class BreakpointStatistics
{
public:
    /**
     * Create breakpoint statistics on graph
     * @param wgraph graph information
     * @param node_name name of node
     * @param forward true to use successor node/edges, false to use predecessors
     */
    BreakpointStatistics(graphs::WalkableGraph const& wgraph, std::string const& node_name, bool forward);

    /**
     * Add edge counts from paragraph output
     * @param paragraph_json JSON output from alignAndDisambiguate
     */
    void addCounts(Json::Value const& paragraph_json);
    int32_t getCount(std::string const& edge_or_allele_name) const;

    /**
     * @return vector of edge names
     */
    std::vector<std::string> const& edgeNames() const { return edge_names; }

    /**
     * @return vector of allele names
     */
    std::vector<std::string> const& alleleNames() const { return allele_names; }

private:
    // edge to allele index and allele names
    std::unordered_map<std::string, std::vector<size_t>> edgename_to_alleles;

    // order of edge and allele names for AlleleCounts arrays below
    std::vector<std::string> edge_names;
    std::map<std::string, size_t> edge_name_to_index;
    std::vector<std::string> allele_names;
    std::map<std::string, size_t> allele_name_to_index;

    // breakpoint count and genotype information for each sample
    AlleleCounts edge_counts; //  count by edge
    AlleleCounts allele_counts; // count by allele
};
};