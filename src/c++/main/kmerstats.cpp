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
#include <map>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <google/protobuf/util/json_util.h>

#include "graphs/GraphSearch.hh"
#include "graphs/KmerIndex.hh"

#include "paragraph/Disambiguation.hh"
#include "paragraph/Parameters.hh"

#include "common/Error.hh"

namespace po = boost::program_options;

using namespace paragraph;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "graph-spec,g", po::value<string>(), "JSON file describing the graph")(
        "output,o", po::value<string>(), "Output file name. Will output to stdout if omitted.")(
        "output-alignments,a", po::value<string>()->default_value(""),
        "Output kmer matches as alignments into this JSON file")(
        "reference,r", po::value<string>(),
        "FASTA with reference genome")("kmer-size,k", po::value<int>()->default_value(3), "Kmer length")(
        "positions,p", po::value<bool>()->default_value(false), "Output positions as well as counts")(
        "convert-reference-to-sequence,C", po::value<bool>()->default_value(false),
        "Convert reference nodes to sequence nodes")(
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
            "kmerstats", vm["log-file"].as<string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<string>().c_str());
        logger = LOG();

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

        Json::Value input;
        {
            std::ifstream input_file(graph_spec_path);
            if (!input_file.good())
            {
                error("Cannot open %s", graph_spec_path.c_str());
            }
            Json::Reader reader;
            reader.parse(input_file, input);
        }

        graphs::Graph g;
        graphs::fromJson(input, reference_path, g);
        graphs::WalkableGraph wg(g);
        graphs::KmerIndex index(wg, vm["kmer-size"].as<int>());

        Json::Value output;
        output["kmers"] = Json::arrayValue;

        if (!input.isMember("alignments"))
        {
            input["alignments"] = Json::arrayValue;
        }

        auto add_alignment = [&input, &index, &wg](std::string kmer, uint64_t node, int pos, bool unique) {
            static int id = 0;

            common::Read read;
            read.setCoreInfo(std::string("k") + std::to_string(id++), kmer, std::string(kmer.size(), '#'));
            const auto gcigar = graphs::prefixMatch(wg, node, pos, kmer);
            read.set_graph_cigar(gcigar);
            read.set_graph_pos(pos);
            read.set_graph_mapping_status(reads::MAPPED);
            read.set_graph_mapq(unique ? 60 : 0);
            read.set_graph_alignment_score(static_cast<google::protobuf::int32>(kmer.size()));
            read.set_is_graph_alignment_unique(unique);

            string str;
            google::protobuf::util::MessageToJsonString(*((google::protobuf::Message*)&read), &str);
            Json::Reader reader;
            Json::Value val;
            reader.parse(str, val);
            input["alignments"].append(val);
        };

        for (auto const& k : index.kmers())
        {
            Json::Value kmerinfo;
            kmerinfo["s"] = k;
            kmerinfo["n"] = index.kmerCount(k);

            if (vm["positions"].as<bool>() || !vm["output-alignments"].as<string>().empty())
            {
                if (vm["positions"].as<bool>())
                {
                    kmerinfo["matches"] = Json::arrayValue;
                }
                index.search(k);
                for (const auto& pos : index.matches())
                {
                    if (vm["positions"].as<bool>())
                    {
                        Json::Value json_pos;
                        json_pos["node"] = (Json::UInt64)pos->node();
                        json_pos["pos"] = pos->pos();
                        kmerinfo["matches"].append(json_pos);
                    }
                    if (!vm["output-alignments"].as<string>().empty())
                    {
                        add_alignment(k, pos->node(), pos->pos(), index.count() > 1);
                    }
                }
            }

            output["kmers"].append(kmerinfo);
        }
        output["badkmers"] = Json::arrayValue;
        for (auto const& k : index.badkmers())
        {
            Json::Value kmerinfo;
            kmerinfo["s"] = k;
            kmerinfo["n"] = index.kmerCount(k);
            for (const auto& pos : index.matches())
            {
                if (vm["positions"].as<bool>())
                {
                    Json::Value json_pos;
                    json_pos["node"] = (Json::UInt64)pos->node();
                    json_pos["pos"] = pos->pos();
                    kmerinfo["matches"].append(json_pos);
                }
                if (!vm["output-alignments"].as<string>().empty())
                {
                    add_alignment(k, pos->node(), pos->pos(), false);
                }
            }
            output["badkmers"].append(kmerinfo);
        }
        if (output_path.empty())
        {
            Json::StyledWriter fastWriter;
            std::cout << fastWriter.write(output);
        }
        else
        {
            Json::FastWriter fastWriter;
            std::ofstream output_stream(output_path);
            output_stream << fastWriter.write(output);
        }
        if (!vm["output-alignments"].as<string>().empty())
        {
            if (vm["convert-reference-to-sequence"].as<bool>())
            {
                for (int n = 0; n < (int)input["nodes"].size(); ++n)
                {
                    // convert reference nodes into sequence nodes
                    if (input["nodes"][n].isMember("reference"))
                    {
                        input["nodes"][n]["reference_location"] = input["nodes"][n]["reference"];
                        input["nodes"][n].removeMember("reference");
                        input["nodes"][n]["sequence"] = wg.node(static_cast<uint64_t>(n))->sequence();
                    }
                }
            }
            Json::FastWriter fastWriter;
            std::ofstream output_stream(vm["output-alignments"].as<string>());
            output_stream << fastWriter.write(input);
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
