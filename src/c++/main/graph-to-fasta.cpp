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
 * \brief Produces path sequences in fasta format
 *
 * \file paragraph.cpp
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/Region.hh"
#include "grm/GraphInput.hh"
#include "json/json.h"

#include "common/Error.hh"

namespace po = boost::program_options;

void dumpPath(
    const graphtools::Graph& graph, const std::unordered_map<std::string, const Json::Value&>& nodes,
    const graphtools::Path& path, std::ostream& os, const std::string& file)
{
    std::string path_id = path.encode();
    std::string cigar;
    std::string sequence;
    common::Region lastRegion;
    //    std::string dbg;
    for (const auto& n : path.nodeIds())
    {
        const auto nodeName = graph.nodeName(n);
        //        dbg += nodeName + ' ';
        if ("source" == nodeName || "sink" == nodeName)
        {
            continue;
        }
        const auto& node = nodes.find(nodeName);
        assert(node != nodes.cend());

        const std::string& s = graph.nodeSeq(n);

        sequence += s;
        if (node->second.isMember("reference"))
        {
            const std::string reg = node->second["reference"].asString();
            common::Region region(reg);
            if (region.chrom == lastRegion.chrom)
            {
                if (region.start == lastRegion.end + 1)
                {
                    lastRegion.end = region.end;
                    continue;
                }
                else if (region.start > lastRegion.end + 1)
                {
                    cigar += std::to_string(lastRegion.length()) + "M" + std::to_string(region.start - lastRegion.start)
                        + "D";
                    lastRegion = region;
                    continue;
                }
            }
            else
            {
                if (!lastRegion.chrom.empty())
                {
                    cigar += std::to_string(lastRegion.end + 1 - lastRegion.start) + "M,";
                }
            }
            lastRegion = region;
            cigar += lastRegion.chrom + '+' + std::to_string(lastRegion.start + 1) + '-';
        }
        else
        {
            cigar += std::to_string(lastRegion.length()) + "M" + std::to_string(s.length()) + "I";
        }
    }

    cigar += std::to_string(lastRegion.length()) + "M";
    //    os << '>' << dbg << std::endl;
    os << '>' << path_id << '_' << cigar << "\n" << sequence << std::endl;
    if (!os)
    {
        error("ERROR: Failed to write output to '%s' error: '%s'", file.c_str(), std::strerror(errno));
    }
}

void dumpPathContigs(
    const Json::Value& root, const std::string& referencePath,
    //                     TargetRegion& targetRegion,
    std::ostream& os, const std::string& file)
{
    const Json::Value& pathsIn = root["paths"];
    // Initialize the graph aligner.
    graphtools::Graph graph = grm::graphFromJson(root, referencePath);

    std::unordered_map<std::string, const Json::Value&> nodes;

    for (const auto& n : root["nodes"])
    {
        const std::string nodeName = n["name"].asString();
        assert(!nodes.count(nodeName));
        nodes.emplace(nodeName, n);
    }

    const std::list<graphtools::Path> paths = grm::pathsFromJson(&graph, pathsIn);

    for (const auto& path : paths)
    {
        dumpPath(
            graph, nodes, path,
            // targetRegion,
            os, file);
    }
}

int main(int argc, char const* argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "graph-spec,g", po::value<std::vector<std::string>>()->multitoken(), "JSON file(s) describing the graph")(
        "output-file,o", po::value<std::string>(),
        "Output file name. "
        "Will output to stdout if '-' or neither of output-file or output-folder provided.")(
        "output-folder,O", po::value<std::string>(),
        "Output folder path. paragraph will attempt to create "
        "the folder but not the entire path. Will output to stdout if neither of output-file or "
        "output-folder provided. If specified, paragraph will produce one output file for each "
        "input file bearing the same name.")("reference,r", po::value<std::string>(), "FASTA with reference genome")(
        "log-level", po::value<std::string>()->default_value("info"), "Set log level (error, warning, info).")(
        "log-file", po::value<std::string>()->default_value(""), "Log to a file instead of stderr.")(
        "log-async", po::value<bool>()->default_value(true), "Enable / disable async logging.");

    po::variables_map vm;
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
            "paragraph", vm["log-file"].as<std::string>().c_str(), vm["log-async"].as<bool>(),
            vm["log-level"].as<std::string>().c_str());

        std::vector<std::string> graph_spec_paths;
        if (vm.count("graph-spec") != 0u)
        {
            graph_spec_paths = vm["graph-spec"].as<std::vector<std::string>>();
            LOG()->info("Graph spec: {}", boost::join(graph_spec_paths, ","));
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

        std::ofstream outputFileStream;
        std::string output_file_path;
        if (vm.count("output-file") != 0u)
        {
            output_file_path = vm["output-file"].as<std::string>();
        }
        else if (!vm.count("output-folder"))
        {
            output_file_path = "-";
        }
        if (!output_file_path.empty())
        {
            if ("-" != output_file_path)
            {
                LOG()->info("Output file path: {}", output_file_path);
                outputFileStream.open(output_file_path);
                if (!outputFileStream)
                {
                    error(
                        "ERROR: Failed to open output file '%s'. Error: '%s'", output_file_path.c_str(),
                        std::strerror(errno));
                }
            }
            else
            {
                LOG()->info("Output to stdout");
            }
        }

        std::string output_folder_path;
        if (vm.count("output-folder"))
        {
            output_folder_path = vm["output-folder"].as<std::string>();
            LOG()->info("Output folder path: {}", output_folder_path);
            boost::filesystem::create_directory(output_folder_path);
        }

        std::string reference_path;
        if (vm.count("reference") != 0u)
        {
            reference_path = vm["reference"].as<std::string>();
            LOG()->info("Reference: {}", reference_path);
            assertFileExists(reference_path);
        }
        else
        {
            error("ERROR: Reference genome is missing.");
        }

        for (auto graph_spec_path : graph_spec_paths)
        {
            LOG()->info("Loading parameters {}", graph_spec_path);
            Json::Value root;
            std::ifstream graph_desc(graph_spec_path);
            graph_desc >> root;

            if ("-" == output_file_path)
            {
                dumpPathContigs(
                    root, reference_path,
                    // targetRegion,
                    std::cout, "standard output");
            }
            else if (!output_file_path.empty())
            {
                dumpPathContigs(
                    root, reference_path,
                    // targetRegion,
                    std::cout, output_file_path);
            }

            if (!output_folder_path.empty())
            {
                boost::filesystem::path inputPath(graph_spec_path);
                boost::filesystem::path outputPath = boost::filesystem::path(output_folder_path) / inputPath.filename();
                std::ofstream ofs(outputPath.string());
                if (!ofs)
                {
                    error(
                        "ERROR: Failed to open output file '%s'. Error: '%s'", outputPath.string().c_str(),
                        std::strerror(errno));
                }
                dumpPathContigs(
                    root, reference_path,
                    // targetRegion,
                    std::cout, outputPath.string());
            }
        }
    }
    catch (const std::exception& e)
    {
        if (LOG())
        {
            LOG()->critical(e.what());
        }
        else
        {
            std::cerr << e.what() << std::endl;
        }
        return 1;
    }

    return 0;
}
