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
 * Functions to run graph alignment in grmpy
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

// avoid boost lambda issues with boost::adaptors::transformed on boost 1.53
#define BOOST_RESULT_OF_USE_DECLTYPE

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <chrono>
#include <iostream>
#include <random>
#include <regex>
#include <sstream>
#include <thread>

#include "genotyping/SampleInfo.hh"
#include "grmpy/AlignSamples.hh"

#include "common/JsonHelpers.hh"
#include "paragraph/Disambiguation.hh"

#include "common/Error.hh"

namespace grmpy
{

static void writeAlignments(
    Json::Value& output, const Parameters& parameters, const paragraph::Parameters& paragraph_parameters,
    const std::string& referencePath, genotyping::SampleInfo& sample)
{
    using namespace boost::filesystem;
    using namespace boost::algorithm;
    using namespace boost::adaptors;
    output["sample"] = sample.sample_name();
    output["reference"] = referencePath;

    static const std::regex unsafe_characters{ "[^A-Za-z0-9.-]" };

    const std::string safe_sample_name = std::regex_replace(sample.sample_name(), unsafe_characters, "_");
    const std::string safe_target_regions = std::regex_replace(
        join(
            paragraph_parameters.target_regions()
                | transformed([](common::Region const& region) -> std::string { return std::string(region); }),
            "_"),
        unsafe_characters, "_");

    Json::Value const& graph = paragraph_parameters.description();

    std::string graph_id;
    if (graph.isMember("ID"))
    {
        graph_id = graph["ID"].asString();
    }
    else if (graph.isMember("model_name"))
    {
        graph_id = graph["model_name"].asString();
    }
    else
    {
        graph_id = boost::uuids::to_string(boost::uuids::uuid());
    }
    const std::string safe_graph_id = std::regex_replace(graph_id, unsafe_characters, "_");

    const path output_path = path(parameters.alignment_output_folder())
        / (safe_sample_name + "-" + safe_graph_id + "-" + safe_target_regions + ".json.gz");

    boost::iostreams::basic_file_sink<char> of(output_path.string());
    if (!of.is_open())
    {
        error(
            "ERROR: Failed to open output file '%s'. Error: '%s'", output_path.string().c_str(), std::strerror(errno));
    }

    boost::iostreams::filtering_ostream fos;
    fos.push(boost::iostreams::gzip_compressor());
    fos.push(of);

    fos << common::writeJson(output);
}

/**
 * Run single sample alignment
 * @param sample sample data structure
 */
void alignSingleSample(
    const Parameters& parameters, const std::string& graphPath, const std::string& referencePath,
    common::BamReader& reader, genotyping::SampleInfo& sample)
{
    auto logger = LOG();
    const bool write_alignments = !parameters.alignment_output_folder().empty()
        && boost::filesystem::is_directory(parameters.alignment_output_folder());

    // set up paragraph aligner
    auto output_options = write_alignments
        ? paragraph::Parameters::ALL
        : paragraph::Parameters::NODE_READ_COUNTS | paragraph::Parameters::EDGE_READ_COUNTS
            | paragraph::Parameters::PATH_READ_COUNTS | paragraph::Parameters::DETAILED_READ_COUNTS;
    if (parameters.infer_read_haplotypes())
    {
        output_options |= paragraph::Parameters::HAPLOTYPES;
    }
    else
    {
        output_options = output_options & (~paragraph::Parameters::HAPLOTYPES);
    }

    paragraph::Parameters paragraph_parameters(
        parameters.max_reads(),
        write_alignments ? 3 : parameters.max_reads() + 1, // disable variants unless we actually also write them out
                                                           // by using min reads for a variant > max reads we read
        0.01, parameters.bad_align_frac(), output_options, parameters.path_sequence_matching(),
        parameters.graph_sequence_matching(), parameters.klib_sequence_matching(), parameters.kmer_sequence_matching(),
        false);
    paragraph_parameters.set_threads(static_cast<uint32_t>(parameters.threads()));
    paragraph_parameters.set_kmer_len(parameters.bad_align_uniq_kmer_len());

    logger->info("Loading parameters for sample {} graph {}", sample.sample_name(), graphPath);
    paragraph_parameters.load(graphPath, referencePath);
    logger->info("Done loading parameters");

    common::ReadBuffer all_reads;

    common::extractReads(
        reader, paragraph_parameters.target_regions(), parameters.max_reads(),
        paragraph_parameters.longest_alt_insertion(), all_reads);
    Json::Value output = paragraph::alignAndDisambiguate(paragraph_parameters, all_reads);
    output["bam"] = sample.filename();

    if (write_alignments)
    {
        writeAlignments(output, parameters, paragraph_parameters, referencePath, sample);
    }

    // alignments take a lot of memory and are not required for downstream processing.
    output.removeMember("alignments");
    output.removeMember("node_coverage");
    output.removeMember("path_coverage");
    output.removeMember("phasing");
    output.removeMember("variants");

    sample.set_alignment_data(output);
}
}
