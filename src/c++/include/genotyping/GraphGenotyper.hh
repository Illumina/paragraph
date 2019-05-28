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
#include "graphcore/Graph.hh"
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
    void reset(graphtools::Graph const* graph);

    /**
     * @return the graph (asserts if no graph is set)
     */
    graphtools::Graph const& getGraph() const;

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

    /**
     * Get depth standard deviation for a sample
     * @param sample_index
     * @return std deviation of average depth. if not available, return 0
     */
    double getDepthSD(size_t sample_index) const;

    /**
     * Get sex integer for a sample
     */
    SampleInfo::Sex getSampleSex(size_t sample_index) const;

private:
    /**
     * internal data structure
     */
    struct GraphGenotyperImpl;
    std::unique_ptr<GraphGenotyperImpl> _impl;
};
};
