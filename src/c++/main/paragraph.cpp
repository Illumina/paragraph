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

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "common/Error.hh"
#include "common/StringUtil.hh"

#include "paragraph/Parameters.hh"
#include "paragraph/Workflow.hh"

namespace po = boost::program_options;

using namespace paragraph;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "bam,b", po::value<std::vector<string>>()->multitoken(),
        "Input BAM file(s) for read extraction. We align all reads to all graphs.")(
        "graph-spec,g", po::value<std::vector<string>>()->multitoken(), "JSON file(s) describing the graph(s)")(
        "output-file,o", po::value<string>(),
        "Output file name. Will output to stdout if '-' or neither of output-file or output-folder provided.")(
        "output-folder,O", po::value<string>(),
        "Output folder path. paragraph will attempt to create "
        "the folder but not the entire path. Will output to stdout if neither of output-file or "
        "output-folder provided. If specified, paragraph will produce one output file for each "
        "input file bearing the same name.")(
        "target-regions,T", po::value<string>(),
        "Comma-separated list of target regions, e.g. chr1:1-20,chr2:2-40. "
        "This overrides the target regions in the graph spec.")(
        "exact-sequence-matching", po::value<bool>()->default_value(true)->implicit_value(true),
        "Use exact sequence match to find the best candidate among all available paths")(
        "graph-sequence-matching", po::value<bool>()->default_value(true)->implicit_value(true),
        "Enables smith waterman graph alignment")(
        "kmer-sequence-matching", po::value<bool>()->default_value(false)->implicit_value(true),
        "Use k-mers to find the best candidate among all available paths")(
        "validate-alignments", po::value<bool>()->default_value(false)->implicit_value(true),
        "Use information in the input bam read names to collect statistics about the accuracy of alignments. "
        "Requires bam file produced with simulate-reads.sh")(
        "output-detailed-read-counts", po::value<bool>()->default_value(false),
        "Output detailed read counts not just for paths but also for each node/edge on the paths.")(
        "output-variants,v", po::value<bool>()->default_value(false), "Output variants not present in the graph.")(
        "output-path-coverage", po::value<bool>()->default_value(false), "Output coverage for paths")(
        "output-node-coverage", po::value<bool>()->default_value(false), "Output coverage for nodes")(
        "output-read-haplotypes", po::value<bool>()->default_value(false),
        "Output graph haplotypes supported by reads.")(
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
        "bad-align-nonuniq", po::value<bool>()->default_value(true), "Remove reads that are not mapped uniquely.")(
        "bad-align-frac", po::value<float>()->default_value(0.8f),
        "Fraction of read that needs to be mapped in order for it to be used.")(
        "bad-align-uniq-kmer-len", po::value<int>()->default_value(0),
        "Kmer length for uniqueness check during read filtering.")(
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

        std::vector<string> bam_paths;
        std::vector<string> bam_index_paths;
        if (vm.count("bam") != 0u)
        {
            const auto& tmp_bam_paths = vm["bam"].as<std::vector<string>>();

            for (const auto& bam_path : tmp_bam_paths)
            {
                auto bai_sep_index = bam_path.rfind('#');
                if (bai_sep_index != string::npos)
                {
                    const auto bam_file_path = bam_path.substr(0, bai_sep_index);
                    const auto bai_file_path = bam_path.substr(bai_sep_index + 1);
                    bam_paths.push_back(bam_file_path);
                    bam_index_paths.push_back(bai_file_path);
                }
                else
                {
                    bam_paths.push_back(bam_path);
                    bam_index_paths.push_back("");
                }
            }

            logger->info("Input BAM(s): {}", boost::join(bam_paths, ","));
            assertFilesExist(bam_paths.begin(), bam_paths.end());
        }
        else
        {
            error("ERROR: BAM file is missing.");
        }

        std::vector<string> graph_spec_paths;
        if (vm.count("graph-spec") != 0u)
        {
            graph_spec_paths = vm["graph-spec"].as<std::vector<string>>();
            logger->info("Graph spec: {}", boost::join(graph_spec_paths, ","));
            assertFilesExist(graph_spec_paths.begin(), graph_spec_paths.end());
            if (vm.count("output-folder"))
            {
                // If we're to produce individual output files per input, the input file
                // paths must have unique file names.
                assertFileNamesUnique(graph_spec_paths.begin(), graph_spec_paths.end());
            }
        }
        else
        {
            error("ERROR: File with variant specification is missing.");
        }

        string output_file_path;
        if (vm.count("output-file") != 0u)
        {
            output_file_path = vm["output-file"].as<string>();
        }
        else if (!vm.count("output-folder"))
        {
            output_file_path = "-";
        }

        string output_folder_path;
        if (vm.count("output-folder"))
        {
            output_folder_path = vm["output-folder"].as<string>();
            logger->info("Output folder path: {}", output_folder_path);
            boost::filesystem::create_directory(output_folder_path);
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
        if (vm["output-read-haplotypes"].as<bool>())
        {
            output_options |= Parameters::output_options::HAPLOTYPES;
        }
        if (vm["output-everything"].as<bool>())
        {
            output_options |= Parameters::output_options::ALL;
        }

        Parameters parameters(
            vm["max-reads-per-event"].as<int>(), vm["variant-min-reads"].as<int>(), vm["variant-min-frac"].as<float>(),
            vm["bad-align-frac"].as<float>(), output_options, vm["exact-sequence-matching"].as<bool>(),
            vm["graph-sequence-matching"].as<bool>(), vm["kmer-sequence-matching"].as<bool>(),
            vm["validate-alignments"].as<bool>());

        parameters.set_threads(vm["threads"].as<int>());
        parameters.set_kmer_len(vm["bad-align-uniq-kmer-len"].as<int>());
        parameters.set_remove_nonuniq_reads(vm["bad-align-nonuniq"].as<bool>());

        Workflow workflow(
            1 != bam_paths.size(), bam_paths, bam_index_paths, graph_spec_paths, output_file_path, output_folder_path,
            parameters, reference_path, target_regions);
        workflow.run();
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
