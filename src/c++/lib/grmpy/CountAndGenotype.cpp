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

#include <fstream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/Error.hh"
#include "genotyping/GraphGenotyper.hh"
#include "grm/CountAndGenotype.hh"

namespace po = boost::program_options;
using namespace common;

namespace grm
{

void countAndGenotype(const Parameters& parameters, std::ostream* out_stream)
{
    auto logger = LOG();
    std::unique_ptr<std::ofstream> file_out;
    const std::string output_path = parameters.output_path();
    if (!output_path.empty())
    {
        logger->info("Output path: {}", output_path);
        file_out.reset(new std::ofstream(output_path));
        out_stream = (std::ostream*)file_out.get();
        logger->info("Done initializing output");
    }
    else
    {
        logger->info("Empty output path. Re-direct to std output");
    }

    genotyping::GraphGenotyper graph_genotypes(
        parameters.genotype_error_rate(), parameters.min_overlap_bases(), parameters.max_read_times());
    graph_genotypes.genotypeGraph(
        parameters.input_path(), parameters.reference_path(), parameters.manifest_path(), parameters.use_em());
    logger->info("Genotyping completed");
    if (parameters.output_as_csv())
    {
        graph_genotypes.toCsv(out_stream);
    }
    else
    {
        graph_genotypes.toJson(out_stream);
    }
    logger->info("Finished analysis and printed output");
}
}
