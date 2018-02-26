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

#include <fstream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/Error.hh"
#include "common/JsonHelpers.hh"
#include "genotyping/GraphBreakpointGenotyper.hh"
#include "graphs/WalkableGraph.hh"
#include "grmpy/CountAndGenotype.hh"

namespace po = boost::program_options;
using namespace common;

namespace grmpy
{

/**
 * main function for alignment and genotyping
 *
 * @param parameters Input parameters
 * @param out_stream pointer to the output stream (for testing purpose)
 */
void countAndGenotype(const Parameters& parameters, std::ostream* out_stream)
{
    auto logger = LOG();
    std::unique_ptr<std::ofstream> file_out;
    const std::string& output_path = parameters.output_path();
    if (!output_path.empty())
    {
        logger->info("Output path: {}", output_path);
        file_out.reset(new std::ofstream(output_path));
        out_stream = (std::ostream*)file_out.get();
        logger->info("Done initializing output");
    }
    else
    {
        logger->info("Empty output path. Re-direct to std output or specified outstream.");
    }

    // Initialize walkable graph
    Json::Value root = common::getJSON(parameters.graph_path());
    graphs::Graph graph;
    graphs::fromJson(root, parameters.reference_path(), graph);

    auto wgraph_ptr = std::make_shared<graphs::WalkableGraph>(graph);

    genotyping::GraphBreakpointGenotyper graph_genotyper;
    graph_genotyper.reset(wgraph_ptr);
    graph_genotyper.setParameters(parameters.genotypingParameterPath());

    for (auto& sample_info : parameters.getSamples())
    {
        graph_genotyper.addAlignment(sample_info);
    }

    Json::StyledStreamWriter writer;
    writer.write(*out_stream, graph_genotyper.getGenotypes());

    logger->info("Output data written.");
}
}
