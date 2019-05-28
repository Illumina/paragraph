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
 * \brief Pam depth-based CNV/SV joint caller
 *
 * \file pam.cpp
 * \author Mitch Bekritsky & Peter Krusche & Egor Dolzhenko
 * \email mbekritsky@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "pam/Parameters.hh"

// Error.hh always needs to go last
#include "common/Error.hh"

namespace po = boost::program_options;

using namespace pam;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "bam,b", po::value<string>(), "Input BAM file for read extraction")(
        "depth-spec,g", po::value<string>(), "JSON file describing the regions for depth extraction")(
        "output,o", po::value<string>(), "Output file name. Will output to stdout if omitted.")(
        "target-regions,T", po::value<string>(),
        "Comma-separated list of target regions, e.g. chr1:1-20,chr2:2-40. "
        "This overrides the target regions in the depth spec.")(
        "output-node-coverage", po::value<bool>()->default_value(false), "Output coverage for nodes")(
        "output-alignments,a", po::value<bool>()->default_value(false), "Output alignments for every read (large).")(
        "output-filtered-alignments,A", po::value<bool>()->default_value(false),
        "Output alignments for every read even when it was filtered (larger).")(
        "output-everything,E", po::value<bool>()->default_value(false),
        "Write all information we have into JSON. (=enable all --output-* above)")(
        "max-reads-per-event,M", po::value<int>()->default_value(10000),
        "Maximum number of reads to process for a single event.")(
        "reference,r", po::value<string>(), "FASTA with reference genome")(
        "log-level", po::value<string>()->default_value("info"), "Set log level (error, warning, info).")(
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
            cerr << desc << endl;
            return 1;
        }

        initLogging(
            "pam", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<string>().c_str());
        logger = LOG();

        string bam_path;
        if (vm.count("bam") != 0u)
        {
            bam_path = vm["bam"].as<string>();
            logger->info("Input BAM: {}", bam_path);
            assertFileExists(bam_path);
        }
        else
        {
            error("ERROR: BAM file is missing.");
        }

        string depth_spec_path;
        if (vm.count("depth-spec") != 0u)
        {
            depth_spec_path = vm["depth-spec"].as<string>();
            logger->info("Depth spec: {}", depth_spec_path);
            assertFileExists(depth_spec_path);
        }
        else
        {
            error("ERROR: File with depth specification is missing.");
        }

        string output_path;
        if (vm.count("output") != 0u)
        {
            output_path = vm["output"].as<string>();
            logger->info("Output path: {}", output_path);
        }

        string reference_path;
        if (vm.count("reference") != 0u)
        {
            reference_path = vm["reference"].as<string>();
            logger->info("Reference: {}", reference_path);
            assertFileExists(reference_path);
        }
        else
        {
            error("ERROR: Reference genome is missing.");
        }

        string target_regions;
        if (vm.count("target-regions") != 0u)
        {
            target_regions = vm["target-regions"].as<string>();
            logger->info("Overriding target regions: {}", target_regions);
        }

        int output_options = Parameters::output_options::NODE_READ_COUNTS;
        if (vm["output-alignments"].as<bool>())
        {
            output_options |= Parameters::output_options::ALIGNMENTS;
        }
        if (vm["output-filtered-alignments"].as<bool>())
        {
            output_options |= Parameters::output_options::FILTERED_ALIGNMENTS;
        }
        if (vm["output-node-coverage"].as<bool>())
        {
            output_options |= Parameters::output_options::NODE_COVERAGE;
        }
        if (vm["output-everything"].as<bool>())
        {
            output_options |= Parameters::output_options::ALL;
        }

        Parameters parameters(output_options);

        logger->info("Loading parameters");
        parameters.load(bam_path, depth_spec_path, reference_path, target_regions);
        logger->info("Done loading parameters");

        // TODO actual logic for loading event spec, extracting reads, and dumping to JSON
    }
    catch (const std::exception& e)
    {
        if (logger)
        {
            logger->critical(e.what());
        }
        else
        {
            cerr << e.what() << endl;
        }
        return 1;
    }

    return 0;
}
