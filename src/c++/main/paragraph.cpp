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
 * \brief Graph read aligner
 *
 * \file paragraph.cpp
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/Error.hh"
#include "common/ReadExtraction.hh"

#include "paragraph/Disambiguation.hh"
#include "paragraph/Parameters.hh"

namespace po = boost::program_options;

using namespace paragraph;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "bam,b", po::value<string>(),
        "Input BAM file for read extraction")("graph-spec,g", po::value<string>(), "JSON file describing the graph")(
        "output,o", po::value<string>(), "Output file name. Will output to stdout if omitted.")(
        "target-regions,T", po::value<string>(),
        "Comma-separated list of target regions, e.g. chr1:1-20,chr2:2-40. "
        "This overrides the target regions in the graph spec.")(
        "exact-sequence-matching", po::value<bool>()->default_value(true),
        "Switch this off to always use Smith Waterman, don't try to find exact sequence matches first.")(
        "output-detailed-read-counts", po::value<bool>()->default_value(false),
        "Output detailed read counts not just for paths but also for each node/edge on the paths.")(
        "output-variants,v", po::value<bool>()->default_value(false), "Output variants not present in the graph.")(
        "output-path-coverage", po::value<bool>()->default_value(false), "Output coverage for paths")(
        "output-node-coverage", po::value<bool>()->default_value(false), "Output coverage for nodes")(
        "output-alignments,a", po::value<bool>()->default_value(false), "Output alignments for every read (large).")(
        "output-filtered-alignments,A", po::value<bool>()->default_value(false),
        "Output alignments for every read even when it was filtered (larger).")(
        "output-everything,E", po::value<bool>()->default_value(false),
        "Write all information we have into JSON. (=enable all --output-* above)")(
        "max-reads-per-event,M", po::value<int>()->default_value(10000),
        "Maximum number of reads to process for a single event.")(
        "variant-min-reads", po::value<int>()->default_value(3),
        "Minimum number of reads required to report a variant.")(
        "variant-min-frac", po::value<float>()->default_value(0.01f),
        "Minimum fraction of reads required to report a variant.")(
        "bad-align-frac", po::value<float>()->default_value(0.8f),
        "Fraction of read that needs to be mapped in order for it to be used.")(
        "reference,r", po::value<string>(), "FASTA with reference genome")(
        "threads", po::value<int>()->default_value(std::thread::hardware_concurrency()),
        "Number of threads to use for parallel alignment.")(
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
            "paragraph", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
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

        string graph_spec_path;
        if (vm.count("graph-spec") != 0u)
        {
            graph_spec_path = vm["graph-spec"].as<string>();
            logger->info("Graph spec: {}", graph_spec_path);
            assertFileExists(graph_spec_path);
        }
        else
        {
            error("ERROR: File with variant specification is missing.");
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

        int output_options = Parameters::output_options::NODE_READ_COUNTS | Parameters::output_options::EDGE_READ_COUNTS
            | Parameters::output_options::PATH_READ_COUNTS;
        if (vm["output-alignments"].as<bool>())
        {
            output_options |= Parameters::output_options::ALIGNMENTS;
        }
        if (vm["output-filtered-alignments"].as<bool>())
        {
            output_options |= Parameters::output_options::FILTERED_ALIGNMENTS;
        }
        if (vm["output-variants"].as<bool>())
        {
            output_options |= Parameters::output_options::VARIANTS;
        }
        if (vm["output-detailed-read-counts"].as<bool>())
        {
            output_options |= Parameters::output_options::DETAILED_READ_COUNTS;
        }
        if (vm["output-path-coverage"].as<bool>())
        {
            output_options |= Parameters::output_options::PATH_COVERAGE;
        }
        if (vm["output-node-coverage"].as<bool>())
        {
            output_options |= Parameters::output_options::NODE_COVERAGE;
        }
        if (vm["output-everything"].as<bool>())
        {
            output_options |= Parameters::output_options::ALL;
        }

        Parameters parameters(
            vm["max-reads-per-event"].as<int>(), vm["variant-min-reads"].as<int>(), vm["variant-min-frac"].as<float>(),
            vm["bad-align-frac"].as<float>(), output_options, vm["exact-sequence-matching"].as<bool>());

        parameters.set_threads(vm["threads"].as<int>());

        logger->info("Loading parameters");
        parameters.load(bam_path, graph_spec_path, reference_path, target_regions);
        logger->info("Done loading parameters");

        common::ReadBuffer all_reads;
        common::extractReads(
            parameters.bam_path(), parameters.reference_path(), parameters.target_regions(),
            (int)parameters.max_reads(), all_reads);
        const auto output = alignAndDisambiguate(parameters, all_reads);

        Json::FastWriter fastWriter;
        std::ofstream output_stream(output_path);
        output_stream << fastWriter.write(output);
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
