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

#include "common/Error.hh"

#include "grm/CountAndGenotype.hh"
#include "grm/Parameters.hh"

using std::string;
namespace po = boost::program_options;

using namespace grm;

int main(int argc, char* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message.")(
        "paragraph,p", po::value<string>(), "JSON file output produced by Paragraph tool.")(
        "reference,r", po::value<string>(), "Reference genome fasta file.")(
        "manifest,m", po::value<string>(), "Manifest of samples with path and bam stats.")(
        "output,o", po::value<string>(), "Output file name. Will output tabular format to stdout if omitted.")(
        "print-csv", po::value<bool>()->default_value(false), "print output as csv instead of json")(
        "log-level", po::value<string>()->default_value("info"), "Set log level (error, warning, info).")(
        "genotype-error-rate", po::value<double>()->default_value(0.01),
        "Fixed genotype error rate for breakpoint genotyping")(
        "min-overlap-bases", po::value<int>()->default_value(16),
        "Minimum overlap bases used in estimating poisson model parameters.")(
        "max-read-times", po::value<int>()->default_value(40),
        "Max times of total reads in one sample for a breakpoint. Multiplied by depth.")(
        "useEM", po::value<bool>()->default_value(false),
        "use Expectation-Maximization for genotyping. Will be more accurate but take more time")(
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

        string paragraph_input_path;
        if (vm.count("paragraph") != 0u)
        {
            paragraph_input_path = vm["paragraph"].as<string>();
            logger->info("Paragraph result as input: {}", paragraph_input_path);
            assertFileExists(paragraph_input_path);
        }
        else
        {
            error("Error: Paragraph output is missing.");
        }

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

        string output_path = vm["output"].as<string>();
        Parameters parameters(
            vm["print-csv"].as<bool>(), vm["genotype-error-rate"].as<double>(), vm["min-overlap-bases"].as<int>(),
            vm["max-read-times"].as<int>(), vm["useEM"].as<bool>());
        logger->info("Loading parameters");
        parameters.load(paragraph_input_path, reference_path, manifest_path, output_path);
        logger->info("Done loading parameters");
        countAndGenotype(parameters, &std::cout);
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
