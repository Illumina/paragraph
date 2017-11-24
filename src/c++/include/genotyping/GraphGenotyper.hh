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
 * \Genotyper for the graph-represented whole variant
 *  \Most likely genotype is voted from all breakpoint genotypes
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/BreakpointStat.hh"
#include "genotyping/Genotype.hh"
#include "genotyping/Idxdepth.hh"
#include "genotyping/VariantGenotype.hh"
#include "graphs/WalkableGraph.hh"
#include "json/json.h"
#include <map>
#include <memory>
#include <vector>

namespace genotyping
{

class GraphGenotyper
{
public:
    GraphGenotyper(
        const double genotype_error_rate, const int min_overlap_bases, const int max_read_times,
        const std::string reference_allele_name = "REF");

    ~GraphGenotyper();

    /**
     * main function for getting whole
     */
    void genotypeGraph(
        const std::string& paragraph_input_path, const std::string& reference_path, const std::string& manifest_path,
        bool use_em);

    // output functions (different format)
    void toCsv(std::ostream* out);

    void toJson(std::ostream* out);

private:
    /**
     * load graph description, sequence index and read counts from input Json
     */
    void loadGraphAndCounts(const std::string& paragraph_input_path, const std::string& reference_path);

    /**
     *  load sequence names and put reference sequence name at the beginning.
     */
    void loadSeqNames(Json::Value& paragraph_json);

    /**
     *  load breakpoints, sequence index and edge names
     */
    void loadBreakpointInfo(graphs::WalkableGraph& wgraph, Json::Value& paragraph_json);

    /**
     * load the map edge_name->sequence index
     * with the map, breakpoint genotype will be represented by sequence index
     */
    std::map<std::string, std::vector<uint64_t>>
    generateEgeNameToSeqIndexMap(Json::Value& paragraph_json, bool source_exist);

    /**
     * special handling of graph with source & sink
     */
    std::vector<uint64_t> removeSourceSink(const std::vector<uint64_t>& node_vec, uint64_t source, uint64_t sink);

    /**
     *  load counts of all BP and samples from paragraph json
     */
    void loadEdgeCounts(Json::Value& paragraph_json);

    /**
     * generate breakpoint specific genotypes from loaded counts
     */
    void computeBreakpointGenotypes(bool use_em);

    /**
     * perform genotyping for the whole variant
     */
    void computeVariantGenotypes();

    /**
     *  append pooled stats as HWE/AF/CR
     */
    void updateVariantStats();

    std::vector<std::string> getUniqueEdgeNames();

    struct GraphGenotyperImpl;
    std::unique_ptr<GraphGenotyperImpl> _impl;
};
};
