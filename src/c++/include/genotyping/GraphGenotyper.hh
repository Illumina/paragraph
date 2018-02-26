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
 * Genotyper for the graph-represented whole variants across many samples
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/GenotypingParameters.hh"
#include "genotyping/SampleInfo.hh"
#include "graphs/WalkableGraph.hh"
#include "json/json.h"

#include <map>
#include <memory>
#include <vector>

namespace genotyping
{

/**
 * Graph genotyper interface and basic functionality
 */
class GraphGenotyper
{
public:
    GraphGenotyper();
    GraphGenotyper(GraphGenotyper const&) = delete;
    GraphGenotyper(GraphGenotyper&&) = delete;
    virtual ~GraphGenotyper();
    GraphGenotyper& operator=(GraphGenotyper const&) = delete;
    GraphGenotyper& operator=(GraphGenotyper&&) = delete;

    /**
     * Set the graph we genotype on
     * @param graph our graph to genotype
     */
    void reset(std::shared_ptr<graphs::WalkableGraph> graph);

    /**
     * @return the graph (asserts if no graph is set)
     */
    graphs::WalkableGraph const& getGraph() const;

    /**
     * Add a alignment and depth information for a sample
     */
    void addAlignment(SampleInfo const& sampleinfo);

    /**
     * This runs runGenotyping.
     * @return JSON encoded set of genotypes for all alignments that were added
     */
    Json::Value getGenotypes();

    /**
     * Function to set parameter values in derived classes
     * @param parameters set of parameters
     */
    virtual void setParameters(const std::string& genotyping_parameter_path) {}

protected:
    /**
     * Implemented by derived classes which implement genotyping
     */
    virtual void runGenotyping() = 0;

    /**
     * Set the genotype for a particular sample
     *
     * @param samplename sample name
     * @param breakpointname name of breakpoint ("" for combined GT)
     * @param genotype breakpoint genotype
     */
    void setGenotype(const std::string& samplename, const std::string& breakpointname, Genotype genotype);

    /**
     * Get the genotype for a particular sample
     *
     * @param samplename sample name
     * @param breakpointname name of breakpoint ("" for combined GT)
     * @return the genotype
     */
    Genotype getGenotype(const std::string& samplename, const std::string& breakpointname) const;

    /**
     * @return an ordered list of sample names
     */
    std::vector<std::string> const& sampleNames() const;

    /**
     * @return a list of breakpoint names
     */
    std::list<std::string> const& breakpointNames() const;

    /**
     * @return a list of allele names
     */
    std::vector<std::string> const& alleleNames() const;

    /**
     * Get the alignment read counts
     * @param sample_index index of sample (name is in sampleNames[sample_index])
     * @param breakpoint the name of the breakpoint
     * @param edge_or_allele_name name of edge or allele
     * @return alignment result for sample
     */
    int32_t getCount(size_t sample_index, std::string const& breakpoint, std::string const& edge_or_allele_name) const;

    /**
     * Get the depth data for a sample
     * @param sample_index index of sample (name is in sampleNames[sample_index])
     * @return pair of expected mean depth and read length
     */
    std::pair<double, int> const& getDepthAndReadlength(size_t sample_index) const;

private:
    /**
     * internal data structure
     */
    struct GraphGenotyperImpl;
    std::unique_ptr<GraphGenotyperImpl> _impl;
};
};
