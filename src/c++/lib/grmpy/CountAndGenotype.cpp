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
 * Functions to run counting and genotyping in grmpy
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "grmpy/CountAndGenotype.hh"
#include "common/JsonHelpers.hh"
#include "genotyping/GraphBreakpointGenotyper.hh"
#include "graphcore/Graph.hh"
#include "grm/GraphInput.hh"

#include <sstream>

namespace grmpy
{

/**
 * main function for alignment and genotyping
 *
 * @param graphPath               If empty the alignment data of the first sample is used as graph
 * @param genotypingParameterPath path to genotyper settings
 * @param samples                 Collection of samples to genotype. Cannot be empty
 */
Json::Value countAndGenotype(
    const std::string& graphPath, const std::string& referencePath, const std::string& genotypingParameterPath,
    const genotyping::Samples& samples)
{
    LOG()->info("Running genotyper");
    // Initialize walkable graph
    Json::Value root = graphPath.empty() ? samples.front().get_alignment_data() : common::getJSON(graphPath);
    graphtools::Graph graph = grm::graphFromJson(root, referencePath, true);

    unsigned int male_ploidy = 2;
    unsigned int female_ploidy = 2;

    for (auto t_region : root["target_regions"])
    {
        std::stringstream ss(t_region.asString());
        std::string chrom;
        getline(ss, chrom, ':');
        if (chrom == "chrX" || chrom == "X")
        {
            male_ploidy = 1;
        }
        else if (chrom == "chrY" || chrom == "Y")
        {
            male_ploidy = 1;
            female_ploidy = 1;
        }
    }

    genotyping::GraphBreakpointGenotyper graph_genotyper(male_ploidy, female_ploidy);
    graph_genotyper.reset(&graph);

    graph_genotyper.setParameters(genotypingParameterPath);

    for (const genotyping::SampleInfo& sample_info : samples)
    {
        graph_genotyper.addAlignment(sample_info);
    }

    Json::Value ret = graph_genotyper.getGenotypes();

    LOG()->info("Done running genotyper");
    return ret;
}

} // namespace grmpy
