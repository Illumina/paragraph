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
 * Data structure for one breakpoint info across samples
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/SampleInfo.hh"
#include "graphcore/Graph.hh"
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
    BreakpointStatistics(graphtools::Graph const& graph, graphtools::NodeId node_name, bool forward);

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
     * @return vector of canonical allele names (i.e. alleles which represent different
     *         subsets of edges on the breakpoint)
     */
    std::vector<std::string> const& canonicalAlleleNames() const { return canonical_allele_names; }

    /**
     * @return vector of all allele names -- some alleles may be duplicated
     */
    std::vector<std::string> const& allAlleleNames() const { return all_allele_names; }

    /**
     * Get the canonical allele name for an allele
     * @param alleleName name of allele on the graph
     * @return the corresponding canonical allele name
     */
    std::string getCanonicalAlleleName(std::string const& alleleName) const
    {
        return allele_name_to_canonical_allele_name.at(alleleName);
    }

private:
    // edge to allele index and allele names
    std::unordered_map<std::string, std::vector<size_t>> edgename_to_alleles;

    // order of edge and allele names for AlleleCounts arrays below
    std::vector<std::string> edge_names;
    std::map<std::string, size_t> edge_name_to_index;
    std::vector<std::string> canonical_allele_names;
    std::map<std::string, size_t> allele_name_to_index;
    std::vector<std::string> all_allele_names;
    std::map<std::string, std::string> allele_name_to_canonical_allele_name;

    // breakpoint count and genotype information for each sample
    AlleleCounts edge_counts; //  count by edge
    AlleleCounts allele_counts; // count by allele
};
};