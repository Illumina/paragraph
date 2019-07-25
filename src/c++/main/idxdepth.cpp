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
 * \brief Index-based depth estimation for BAM and CRAM
 *
 * \file idxdepth.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "json/json.h"
#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/JsonHelpers.hh"
#include "idxdepth/DepthEstimation.hh"
#include "idxdepth/IndexBinning.hh"
#include "idxdepth/Parameters.hh"

#include "common/Error.hh"

namespace po = boost::program_options;

using idxdepth::Parameters;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
            ("help,h", "produce help message")
            ("bam,b", po::value<string>(), "BAM / CRAM input file")
            ("bam-index", po::value<string>()->default_value(""), "BAM / CRAM index file when not at default location.")
            ("output,o", po::value<string>(), "Output file name. Will output to stdout if omitted.")
            ("output-bins,O", po::value<string>(), "Output binned coverage in tsv format.")
            ("reference,r", po::value<string>(), "FASTA with reference genome")
            ("altcontig",po::value<bool>()->default_value(false), "Include ALT contigs in estimation")
            ("include-regex,I", po::value<string>()->default_value(""), "Regex to identify contigs to include")
            ("autosome-regex", po::value<string>()->default_value("(chr)?[1-9][0-9]?"),
             "Regex to identify autosome chromosome names (default: '(chr)?[1-9][0-9]?'")
            ("sex-chromosome-regex", po::value<string>()->default_value("(chr)?[XY]?"),
             "Regex to identify sex chromosome names (default: '(chr)?[XY]?'")
            ("threads", po::value<int>()->default_value(std::thread::hardware_concurrency()),
             "Number of threads to use for parallel estimation.")
            ("log-level", po::value<string>()->default_value("info"), "Set log level (error, warning, info).")
            ("log-file", po::value<string>()->default_value(""), "Log to a file instead of stderr.")
            ("log-async", po::value<bool>()->default_value(true), "Enable / disable async logging.");
    // clang-format on

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
            "idxdepth", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<string>().c_str());
        logger = LOG();

        string bam_path;
        if (vm.count("bam") != 0u)
        {
            bam_path = vm["bam"].as<string>();
            logger->info("BAM: {}", bam_path);
            assertFileExists(bam_path);
        }
        else
        {
            error("ERROR: File with variant specification is missing.");
        }
        const string bam_index_path = vm["bam-index"].as<string>();

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

        string output_path;
        if (vm.count("output") != 0u)
        {
            output_path = vm["output"].as<string>();
            logger->info("Output path: {}", output_path);
        }

        string output_path_bins;
        if (vm.count("output-bins") != 0u)
        {
            output_path_bins = vm["output-bins"].as<string>();
            logger->info("Output path for binned coverage: {}", output_path_bins);
        }

        Parameters parameters(bam_path, bam_index_path, reference_path, vm["altcontig"].as<bool>());
        parameters.set_include_regex(vm["include-regex"].as<string>());
        parameters.set_autosome_regex(vm["autosome-regex"].as<string>());
        parameters.set_sex_chromosome_regex(vm["sex-chromosome-regex"].as<string>());
        parameters.set_threads(vm["threads"].as<int>());

        Json::Value output = idxdepth::estimateDepths(parameters);

        if (output_path.empty())
        {
            std::cout << common::writeJson(output);
        }
        else
        {
            std::ofstream output_stream(output_path);
            output_stream << common::writeJson(output, false);
        }

        if (!output_path_bins.empty())
        {
            std::vector<idxdepth::IndexBin> bins;
            idxdepth::getIndexBins(bam_path, bins);

            std::ofstream bin_output(output_path_bins);

            std::list<std::string> header{
                "id", "CHROM", "POS", "END", "SLICES", "OVERLAPPING_BYTES", "ADJUSTED_BYTES", "NORMALIZED_DEPTH"
            };

            bin_output << boost::algorithm::join(header, "\t") << endl;

            for (auto const& bin : bins)
            {
                bin_output << bin.bin_id << "\t" << bin.chrom << "\t" << bin.start + 1 << "\t" << bin.end + 1 << "\t"
                           << bin.slices << "\t" << bin.overlapping_bytes << "\t" << bin.adjusted_bytes << "\t"
                           << bin.normalized_depth << endl;
            }
        }
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
