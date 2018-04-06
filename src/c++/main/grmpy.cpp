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

#include "grmpy/Parameters.hh"
#include "grmpy/Workflow.hh"

#include "common/Error.hh"
#include "common/Program.hh"

using std::string;
namespace po = boost::program_options;

using namespace grmpy;

class Options : public common::Options
{
public:
    Options();

    common::Options::Action parse(int argc, const char* argv[]);
    void postProcess(boost::program_options::variables_map& vm);

    string reference_path;
    std::vector<std::string> graph_spec_paths;
    string output_file_path;
    string output_folder_path;
    genotyping::Samples manifest;
    string genotyping_parameter_path;
    int sample_threads = std::thread::hardware_concurrency();
    int max_reads_per_event = 10000;
    float bad_align_frac = 0.8f;
    bool exact_sequence_matching = true;
    bool graph_sequence_matching = true;
    bool kmer_sequence_matching = false;
    int bad_align_uniq_kmer_len = 0;
    bool gzip_output = false;

    std::string usagePrefix() const { return "grmpy -r <reference> -g <graphs> -m <manifest> [optional arguments]"; }
};

Options::Options()
{
    namedOptions_.add_options()("help,h", "Produce help message.")(
        "reference,r", po::value<string>(), "Reference genome fasta file.")(
        "graph-spec,g", po::value<std::vector<string>>()->multitoken(), "JSON file(s) describing the graph(s)")(
        "genotyping-parameters,G", po::value<string>(), "JSON file with genotyping model parameters")(
        "manifest,m", po::value<string>(), "Manifest of samples with path and bam stats.")(
        "output-file,o", po::value<string>(),
        "Output file name. Will output tabular format to stdout if omitted or '-'.")(
        "output-folder,O", po::value<string>(),
        "Output folder path. paragraph will attempt to create "
        "the folder but not the entire path. Will output to stdout if neither of output-file or "
        "output-folder provided. If specified, paragraph will produce one output file for each "
        "input file bearing the same name.")(
        "max-reads-per-event,M", po::value<int>(&max_reads_per_event)->default_value(max_reads_per_event),
        "Maximum number of reads to process for a single event.")(
        "bad-align-frac", po::value<float>(&bad_align_frac)->default_value(bad_align_frac),
        "Fraction of read that needs to be mapped in order for it to be used.")(
        "exact-sequence-matching",
        po::value<bool>(&exact_sequence_matching)->default_value(exact_sequence_matching)->implicit_value(true),
        "Use exact sequence match to find the best candidate among all available paths")(
        "graph-sequence-matching",
        po::value<bool>(&graph_sequence_matching)->default_value(graph_sequence_matching)->implicit_value(true),
        "Enables smith waterman graph alignment")(
        "kmer-sequence-matching",
        po::value<bool>(&kmer_sequence_matching)->default_value(kmer_sequence_matching)->implicit_value(true),
        "Use kmer aligner.")(
        "bad-align-uniq-kmer-len", po::value<int>(&bad_align_uniq_kmer_len)->default_value(bad_align_uniq_kmer_len),
        "Kmer length for uniqueness check during read filtering.")(
        "sample-threads,t", po::value<int>(&sample_threads)->default_value(sample_threads),
        "Number of threads for parallel sample processing.")(
        "gzip-output,z", po::value<bool>(&gzip_output)->default_value(gzip_output)->implicit_value(true),
        "gzip-compress output files."
        "If -O is used, output file names are appended with .gz");
    unnamedOptions_.add_options()("alignment-threads", po::value<int>()->default_value(1), "Deprecated. Don't use.");
}

/**
 * \brief remembers the original argv array and hands over to the base implementation
 */
Options::Action Options::parse(int argc, const char* argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    std::cerr << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;

    common::Options::Action ret = common::Options::parse(argc, argv);
    return ret;
}

void Options::postProcess(boost::program_options::variables_map& vm)
{
    std::shared_ptr<spdlog::logger> logger = LOG();

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

    if (vm.count("output-file") != 0u)
    {
        output_file_path = vm["output-file"].as<string>();
    }
    else if (!vm.count("output-folder"))
    {
        output_file_path = "-";
    }

    if (vm.count("output-folder"))
    {
        output_folder_path = vm["output-folder"].as<string>();
        logger->info("Output folder path: {}", output_folder_path);
        boost::filesystem::create_directory(output_folder_path);
    }

    if (vm.count("manifest"))
    {
        const string manifest_path = vm["manifest"].as<string>();
        logger->info("Manifest path: {}", manifest_path);
        assertFileExists(manifest_path);
        manifest = genotyping::loadManifest(manifest_path);
        if (graph_spec_paths.empty())
        {
            // If no graphs given, all manifest samples must have paragraph column set
            for (const genotyping::SampleInfo& sample : manifest)
            {
                if (sample.get_alignment_data().isNull())
                {
                    error(
                        "Error: No graphs given on the command line and sample '%s' has empty paragraph "
                        "column in the manifest.",
                        sample.sample_name().c_str());
                }
            }
        }
        else if (1 < graph_spec_paths.size())
        {
            for (const genotyping::SampleInfo& sample : manifest)
            {
                if (!sample.get_alignment_data().isNull())
                {
                    error(
                        "ERROR: Pre-aligned samples are allowed only when genotyping for a single variant. %d "
                        "graphs provided.",
                        graph_spec_paths.size());
                }
            }
        }
    }
    else
    {
        error("Error: Manifest file is missing.");
    }

    if (vm.count("genotyping-parameters"))
    {
        genotyping_parameter_path = vm["genotyping-parameters"].as<string>();
    }

    if (vm.count("alignment-threads"))
    {
        LOG()->warn("--alignment-threads is deprecated and ignored");
    }
}

static void runGrmpy(const Options& options)
{
    Parameters parameters(
        options.sample_threads, options.max_reads_per_event, options.bad_align_frac, options.exact_sequence_matching,
        options.graph_sequence_matching, options.kmer_sequence_matching, options.bad_align_uniq_kmer_len);
    grmpy::Workflow workflow(
        options.graph_spec_paths, options.genotyping_parameter_path, options.manifest, options.output_file_path,
        options.output_folder_path, options.gzip_output, parameters, options.reference_path);
    workflow.run();
}

int main(int argc, const char* argv[])
{
    common::run(runGrmpy, "Genotyping", argc, argv);
    return 0;
}
