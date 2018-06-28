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
 * \brief Graph kmer statistics tool
 *
 * \file kmerstats.cpp
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/JsonHelpers.hh"
#include "graphalign/KmerIndex.hh"
#include "graphalign/KmerIndexOperations.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/Parameters.hh"

#include "common/Error.hh"

namespace po = boost::program_options;

using namespace paragraph;
using std::cerr;
using std::endl;
using std::list;
using std::string;

int main(int argc, char const* argv[])
{

    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
    ("help,h", "produce help message")
    ("graph-spec,g", po::value<string>()->required(), "JSON file describing the graph")
    ("output,o", po::value<string>()->default_value(""), "Output file name. Will output to stdout if omitted.")
    ("reference,r", po::value<string>()->required(), "FASTA with reference genome")
    ("kmer-length,k", po::value<int>()->default_value(-1), "Kmer length (use negative value to autodetect kmer that covers nodes with minimum number of kmers).")
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
            "kmerstats", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<string>().c_str());
        logger = LOG();

        const string graph_spec_path = vm["graph-spec"].as<string>();
        logger->info("Graph spec: {}", graph_spec_path);
        assertFileExists(graph_spec_path);

        const string reference_path = vm["reference"].as<string>();
        logger->info("Reference: {}", reference_path);
        assertFileExists(reference_path);

        const string output_path = vm["output"].as<string>();

        if (!output_path.empty())
        {
            logger->info("Output path: {}", output_path);
        }

        Json::Value input;
        {
            std::ifstream input_file(graph_spec_path);
            if (!input_file.good())
            {
                error("Cannot open %s", graph_spec_path.c_str());
            }
            input_file >> input;
        }

        int32_t kmer_len = vm["kmer-length"].as<int>();

        graphtools::Graph graph = grm::graphFromJson(input, reference_path);

        if (kmer_len < 0)
        {
            logger->info(
                "Auto-detecting kmer length that covers all nodes + edges with at least {} unique kmers.", -kmer_len);
            kmer_len = graphtools::findMinCoveringKmerLength(
                &graph, static_cast<size_t>(-kmer_len), static_cast<size_t>(-kmer_len));
            if (kmer_len <= 0)
            {
                error("Cannot detect kmer length!");
            }
            logger->info("Auto-detected kmer length is {}", kmer_len);
        }

        graphtools::KmerIndex index(graph, kmer_len);

        Json::Value output;
        output["graph"] = graph_spec_path;
        output["kmer_length"] = kmer_len;
        output["kmers"] = Json::arrayValue;
        output["total_kmers"] = 0;
        output["unique_kmers"] = 0;

        output["kmers_by_node"] = Json::objectValue;

        for (graphtools::NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
        {
            output["kmers_by_node"][graph.nodeName(node_id)]
                = static_cast<Json::UInt64>(index.numUniqueKmersOverlappingNode((uint32_t)node_id));
        }

        if (output_path.empty())
        {
            std::cout << common::writeJson(output);
        }
        else
        {
            std::ofstream output_stream(output_path);
            output_stream << common::writeJson(output, false);
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
