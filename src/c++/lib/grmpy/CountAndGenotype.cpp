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
