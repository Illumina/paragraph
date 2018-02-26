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
 * graph genotyper for graph models
 *
 * \author Sai Chen & Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "spdlog/spdlog.h"

#include "grmpy/AlignSamples.hh"
#include "grmpy/CountAndGenotype.hh"
#include "grmpy/Parameters.hh"

#include "common/Error.hh"

using std::string;
namespace po = boost::program_options;

using namespace grmpy;

int main(int argc, const char* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message.")(
        "reference,r", po::value<string>(),
        "Reference genome fasta file.")("graph-spec,g", po::value<string>(), "JSON file describing the graph")(
        "genotyping-parameters,G", po::value<string>(), "JSON file with genotyping model parameters")(
        "manifest,m", po::value<string>(), "Manifest of samples with path and bam stats.")(
        "output,o", po::value<string>(), "Output file name. Will output tabular format to stdout if omitted.")(
        "max-reads-per-event,M", po::value<int>()->default_value(10000),
        "Maximum number of reads to process for a single event.")(
        "bad-align-frac", po::value<float>()->default_value(0.8f),
        "Fraction of read that needs to be mapped in order for it to be used.")(
        "log-level", po::value<string>()->default_value("info"), "Set log level (error, warning, info).")(
        "sample-threads,t", po::value<int>()->default_value(std::thread::hardware_concurrency()),
        "Number of threads for parallel sample processing.")(
        "alignment-threads", po::value<int>()->default_value(1), "Number of threads for parallel read alignment.")(
        "log-file", po::value<string>()->default_value(""), "Log to a file instead of stderr.")(
        "log-async", po::value<bool>()->default_value(true), "Enable / disable async logging.");

    po::variables_map vm;
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.empty() || (vm.count("help") != 0u))
        {
            std::cerr << desc << std::endl;
            return 1;
        }

        initLogging(
            "Genotyper", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<string>().c_str());
        logger = LOG();

        string reference_path;
        if (vm.count("reference"))
        {
            reference_path = vm["reference"].as<string>();
            logger->info("Reference path: {}", reference_path);
            assertFileExists(reference_path);
        }
        else
        {
            error("Error: Reference genome path is missing.");
        }

        string graph_path;
        if (vm.count("graph-spec"))
        {
            graph_path = vm["graph-spec"].as<string>();
            logger->info("Graph path: {}", graph_path);
            assertFileExists(graph_path);
        }
        else
        {
            error("Error: Graph spec path is missing.");
        }

        string manifest_path;
        if (vm.count("manifest"))
        {
            manifest_path = vm["manifest"].as<string>();
            logger->info("Manifest path: {}", manifest_path);
            assertFileExists(manifest_path);
        }
        else
        {
            error("Error: Manifest file is missing.");
        }
        string genotyping_parameter_path;
        if (vm.count("genotyping-parameters"))
        {
            genotyping_parameter_path = vm["genotyping-parameters"].as<string>();
        }

        const string output_path = vm["output"].as<string>();

        logger->info("Loading parameters");
        Parameters parameters(
            vm["sample-threads"].as<int>(), vm["alignment-threads"].as<int>(), vm["max-reads-per-event"].as<int>(),
            vm["bad-align-frac"].as<float>());
        parameters.load(graph_path, reference_path, manifest_path, output_path, genotyping_parameter_path);
        logger->info("Done loading parameters");

        logger->info("Running alignments");
        alignSamples(parameters);
        logger->info("Done running alignments");

        logger->info("Running genotyper");
        countAndGenotype(parameters);
        logger->info("Done running genotyper");
    }
    catch (const std::exception& e)
    {
        if (logger)
        {
            logger->critical(e.what());
        }
        else
        {
            std::cerr << e.what() << std::endl;
        }
        return 1;
    }

    return 0;
}
