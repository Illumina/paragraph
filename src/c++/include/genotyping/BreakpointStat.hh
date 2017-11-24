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
 * \data structure for one breakpoint info & genotype
 *
 * \author Sai Chen
 * \email schen6@illumina.com
 *
 */

#pragma once

#include "genotyping/EMbreakpointGenotyper.hh"
#include "genotyping/Idxdepth.hh"
#include "json/json.h"
#include <string>
#include <vector>

namespace genotyping
{
class BreakpointStat
{
public:
    BreakpointStat(
        std::string& node_name, bool from_node_fixed, std::vector<std::string>& neighbor_names,
        std::map<std::string, std::vector<uint64_t>>& edge_to_seq_indexes);

    void addBreakpointEdgeCountForOneSample(Json::Value& paragraph_json, const std::string& sample_name);

    void genotype(
        double genotype_error_rate, int sample_size, int max_read_times, int min_overlap_bases, Idxdepth& idxdepth,
        bool use_em);

    Json::Value descriptionsToJson();

    Json::Value sampleGenotypeToJson(int sample, std::vector<std::string>& allele_names);

    Genotype getGenotype(int sample) { return genotypes[sample]; }

    int32_t getEdgeCount(int sample, std::string& edge_name);

    std::vector<std::string> exportEdgeNames() { return edge_names; }

private:
    void recodeGenotypesWithSequenceIndex();

    // basic information for genotyping
    std::vector<std::vector<uint64_t>> haplotype_indexes; // edge -> seq indexes (coule be multiple)
    std::vector<std::string> edge_names;

    // load from Json
    std::vector<EdgeCounts> edge_counts; // sample -> count by haplotype

    // breakpoint genotyes
    std::vector<Genotype> genotypes;

    // pooled stats
    int num_iterations;
    double call_rate;
    double hwe_p;
};
};