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
#include "common/Program.hh"
#include "common/StringUtil.hh"

#include "paragraph/Parameters.hh"
#include "paragraph/Workflow.hh"

using std::string;
namespace po = boost::program_options;

using namespace paragraph;

class Options : public common::Options
{
public:
    Options();

    common::Options::Action parse(const char* moduleName, int argc, const char* argv[]);
    void postProcess(boost::program_options::variables_map& vm) override;

    string reference_path;
    std::vector<string> graph_spec_paths;
    string output_file_path;
    string output_folder_path;
    string target_regions;
    int threads = std::thread::hardware_concurrency();
    bool path_sequence_matching = true;
    bool graph_sequence_matching = true;
    bool klib_sequence_matching = false;
    bool kmer_sequence_matching = false;
    bool gzip_output = false;
    int output_options = Parameters::output_options::NODE_READ_COUNTS | Parameters::output_options::EDGE_READ_COUNTS
        | Parameters::output_options::PATH_READ_COUNTS;
    std::vector<string> bam_paths;
    std::vector<string> bam_index_paths;
    int max_reads_per_event = 10000;
    int variant_min_reads = 3;
    float variant_min_frac = 0.01f;
    float bad_align_frac = 0.8f;
    bool validate_alignments = false;
    int bad_align_uniq_kmer_len = 0;
    bool bad_align_nonuniq = true;

    std::string usagePrefix() const override
    {
        return "paragraph -r <reference> -g <graph(s)> -b <input cram(s)/bam(s)> [optional arguments]";
    }
};

Options::Options()
{
    // clang-format off
    namedOptions_.add_options()
        ("bam,b", po::value<std::vector<string>>()->multitoken(),
         "Input BAM file(s) for read extraction. We align all reads to all graphs.")
        ("graph-spec,g", po::value<std::vector<string>>(&graph_spec_paths)->multitoken(), "JSON file(s) describing the graph(s)")
        ("output-file,o", po::value<string>(&output_file_path),
         "Output file name. Will output to stdout if '-' or neither of output-file or output-folder provided.")
        ("output-folder,O", po::value<string>(&output_folder_path),
         "Output folder path. paragraph will attempt to create "
         "the folder but not the entire path. Will output to stdout if neither of output-file or "
         "output-folder provided. If specified, paragraph will produce one output file for each "
         "input file bearing the same name.")
        ("target-regions,T", po::value<string>(&target_regions),
         "Comma-separated list of target regions, e.g. chr1:1-20,chr2:2-40. "
         "This overrides the target regions in the graph spec.")
        ("path-sequence-matching",
         po::value<bool>(&path_sequence_matching)->default_value(path_sequence_matching),
         "Enable path seeding aligner")
        ("graph-sequence-matching",
         po::value<bool>(&graph_sequence_matching)->default_value(graph_sequence_matching),
         "Enables smith waterman graph alignment")
        ("klib-sequence-matching",
         po::value<bool>(&klib_sequence_matching)->default_value(klib_sequence_matching),
         "Use klib smith-waterman aligner.")
        ("kmer-sequence-matching",
         po::value<bool>(&kmer_sequence_matching)->default_value(kmer_sequence_matching),
         "Use kmer aligner.")
        ("validate-alignments", po::value<bool>(&validate_alignments)->default_value(validate_alignments)->implicit_value(true),
         "Use information in the input bam read names to collect statistics about the accuracy of alignments. "
         "Requires bam file produced with simulate-reads.sh")
        ("output-detailed-read-counts", po::value<bool>()->default_value(false),
         "Output detailed read counts not just for paths but also for each node/edge on the paths.")
        ("output-variants,v", po::value<bool>()->default_value(false), "Output variants not present in the graph.")
        ("output-path-coverage", po::value<bool>()->default_value(false), "Output coverage for paths")
        ("output-node-coverage", po::value<bool>()->default_value(false), "Output coverage for nodes")
        ("output-read-haplotypes", po::value<bool>()->default_value(false), "Output graph haplotypes supported by reads.")
        ("output-alignments,a", po::value<bool>()->default_value(false), "Output alignments for every read (large).")
        ("output-filtered-alignments,A", po::value<bool>()->default_value(false),
         "Output alignments for every read even when it was filtered (larger).")
        ("output-everything,E", po::value<bool>()->default_value(false),
         "Write all information we have into JSON. (=enable all --output-* above)")
        ("max-reads-per-event,M", po::value<int>(&max_reads_per_event)->default_value(max_reads_per_event),
         "Maximum number of reads to process for a single event.")
        ("variant-min-reads", po::value<int>(&variant_min_reads)->default_value(variant_min_reads),
         "Minimum number of reads required to report a variant.")
        ("variant-min-frac", po::value<float>(&variant_min_frac)->default_value(variant_min_frac),
         "Minimum fraction of reads required to report a variant.")
        ("bad-align-nonuniq", po::value<bool>(&bad_align_nonuniq)->default_value(bad_align_nonuniq), "Remove reads that are not mapped uniquely.")
        ("bad-align-frac", po::value<float>(&bad_align_frac)->default_value(bad_align_frac),
         "Fraction of read that needs to be mapped in order for it to be used.")
        ("bad-align-uniq-kmer-len", po::value<int>(&bad_align_uniq_kmer_len)->default_value(bad_align_uniq_kmer_len),
         "Kmer length for uniqueness check during read filtering.")
        ("reference,r", po::value<string>(&reference_path), "Reference genome fasta file.")
        ("threads", po::value<int>(&threads)->default_value(threads), "Number of threads to use for parallel alignment.")
        ("gzip-output,z", po::value<bool>(&gzip_output)->default_value(gzip_output)->implicit_value(true),
         "gzip-compress output files. If -O is used, output file names are appended with .gz");
}

/**
 * \brief remembers the original argv array and hands over to the base implementation
 */
Options::Action Options::parse(const char* moduleName, int argc, const char* argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
#ifdef GRMPY_TRACE
    LOG()->info("argc: {} argv: {}", argc, boost::join(allOptions, " "));
#endif
    common::Options::Action ret = common::Options::parse(moduleName, argc, argv);
    return ret;
}

void Options::postProcess(boost::program_options::variables_map& vm)
{
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

        LOG()->info("Input BAM(s): {}", boost::join(bam_paths, ","));
        assertFilesExist(bam_paths.begin(), bam_paths.end());
    }
    else
    {
        error("ERROR: BAM file is missing.");
    }

    if (!graph_spec_paths.empty())
    {
        LOG()->info("Graph spec: {}", boost::join(graph_spec_paths, ","));
        assertFilesExist(graph_spec_paths.begin(), graph_spec_paths.end());
        if (!output_folder_path.empty())
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

    if (output_file_path.empty() && output_folder_path.empty())
    {
        output_file_path = "-";
    }

    if (!output_folder_path.empty())
    {
        LOG()->info("Output folder path: {}", output_folder_path);
        boost::filesystem::create_directory(output_folder_path);
    }

    if (!reference_path.empty())
    {
        LOG()->info("Reference: {}", reference_path);
        assertFileExists(reference_path);
    }
    else
    {
        error("ERROR: Reference genome is missing.");
    }

    if (!target_regions.empty())
    {
        LOG()->info("Overriding target regions: {}", target_regions);
    }

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
}

static void runParagraph(const Options& options)
{
    Parameters parameters(
        options.max_reads_per_event, options.variant_min_reads, options.variant_min_frac,
        options.bad_align_frac, options.output_options, options.path_sequence_matching, options.graph_sequence_matching,
        options.klib_sequence_matching, options.kmer_sequence_matching, options.validate_alignments);

    parameters.set_threads(options.threads);
    parameters.set_kmer_len(options.bad_align_uniq_kmer_len);
    parameters.set_remove_nonuniq_reads(options.bad_align_nonuniq);

    Workflow workflow(
            1 != options.bam_paths.size(), options.bam_paths, options.bam_index_paths, options.graph_spec_paths,
            options.output_file_path, options.output_folder_path,
            options.gzip_output, parameters, options.reference_path, options.target_regions);
    workflow.run();
}

int main(int argc, const char* argv[])
{
    common::run(runParagraph, "paragraph", argc, argv);
    return 0;
}
