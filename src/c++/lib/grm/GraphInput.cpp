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
 *  \brief Graph Helper
 *
 * \file Graph.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <algorithm>
#include <string>

#include "grm/GraphInput.hh"

#include "common/Error.hh"
#include "common/Fasta.hh"
#include "common/StringUtil.hh"

using graphtools::Graph;
using graphtools::NodeId;
using std::string;

namespace grm
{

/**
 * Initialize graph from JSON
 * @param in Input JSON node
 * @param reference reference fasta name
 * @param out Output graph to initialize
 */
Graph graphFromJson(Json::Value const& in, string const& reference, bool store_ref_sequence)
{
    common::FastaFile ref(reference);
    Json::Value const* in_graph = &in;
    if (in.isMember("graph"))
    {
        in_graph = &in["graph"];
    }

    assert((*in_graph)["nodes"].type() == Json::ValueType::arrayValue);

    Graph result{ (*in_graph)["nodes"].size(), false };

    if ((*in_graph)["edges"].type() != Json::ValueType::nullValue)
    {
        assert((*in_graph)["edges"].type() == Json::ValueType::arrayValue);
    }

    std::map<string, NodeId> node_map;
    for (NodeId i = 0; i < (*in_graph)["nodes"].size(); ++i)
    {
        auto const& in_n = (*in_graph)["nodes"][(int)i];

        const string name = in_n.isMember("name") ? in_n["name"].asString() : string("node-") + std::to_string(i + 1);
        assert(!node_map.count(name));
        node_map[name] = i;
        result.setNodeName(i, name);

        string uc_name = name;
        common::stringutil::toUpper(uc_name);
        const bool is_source_or_sink
            = (i == 0 || i == (*in_graph)["nodes"].size() - 1) && (uc_name == "SOURCE" || uc_name == "SINK");

        assert(in_n.isMember("sequence") || in_n.isMember("reference") || is_source_or_sink);

        if (is_source_or_sink)
        {
            result.setNodeSeq(i, "X");
        }
        else if (in_n.isMember("sequence"))
        {
            result.setNodeSeq(i, in_n["sequence"].asString());
        }
        else
        {
            string reference_sequence;
            if (in_n["reference"].type() == Json::ValueType::stringValue)
            {
                const string reference_location = in_n["reference"].asString();
                reference_sequence = ref.query(reference_location);
            }
            else
            {
                assert(in_n["reference"].type() == Json::ValueType::arrayValue);

                for (const auto& ref_inst : in_n["reference"])
                {
                    assert(ref_inst.type() == Json::ValueType::stringValue);

                    const string reference_location = ref_inst.asString();
                    const string curr_seq(ref.query(reference_location));

                    if (!reference_sequence.empty())
                    {
                        assert(reference_sequence == curr_seq);
                    }
                    reference_sequence = curr_seq;
                }
            }
            if (store_ref_sequence)
            {
                assert(!reference_sequence.empty());
                result.setNodeSeq(i, reference_sequence);
            }
        }
    }

    for (size_t i = 0; i < (*in_graph)["edges"].size(); ++i)
    {
        auto const& in_e = (*in_graph)["edges"][(int)i];
        assert(node_map.count(in_e["from"].asString()));
        assert(node_map.count(in_e["to"].asString()));
        const auto from_node = node_map[in_e["from"].asString()];
        const auto to_node = node_map[in_e["to"].asString()];
        result.addEdge(from_node, to_node);

        for (auto& sequence : in_e["sequences"])
        {
            result.addLabelToEdge(from_node, to_node, sequence.asString());
        }
    }

    // Node label shortcut: add label to all in and out edges
    for (NodeId i = 0; i < (*in_graph)["nodes"].size(); ++i)
    {
        auto const& in_n = (*in_graph)["nodes"][(int)i];
        for (auto& sequence : in_n["sequences"])
        {
            for (auto h : result.predecessors(i))
            {
                result.addLabelToEdge(h, i, sequence.asString());
            }
            for (auto j : result.successors(i))
            {
                result.addLabelToEdge(i, j, sequence.asString());
            }
        }
    }
    return result;
}

/**
 * Read paths from JSON
 * @param graph graph to use for the paths
 * @param in_paths Input JSON node with paths
 */
std::list<graphtools::Path> pathsFromJson(graphtools::Graph const* graph, Json::Value const& in_paths)
{
    std::list<graphtools::Path> paths;

    assert(in_paths.type() == Json::ValueType::arrayValue);

    std::unordered_map<std::string, NodeId> node_id_map;
    for (NodeId node_id = 0; node_id != graph->numNodes(); ++node_id)
    {
        node_id_map[graph->nodeName(node_id)] = node_id;
    }

    for (const auto& path : in_paths)
    {
        auto const& path_nodes = path["nodes"];
        assert(path_nodes.type() == Json::ValueType::arrayValue);

        std::vector<NodeId> nodes;
        nodes.reserve(path_nodes.size());
        for (const auto& node_json : path_nodes)
        {
            const auto& node_name = node_json.asString();
            assert(node_id_map.count(node_name) == 1);
            nodes.push_back(node_id_map[node_name]);
        }

        paths.emplace_back(graph, 0, nodes, graph->nodeSeq(nodes.back()).size() - 1);
    }
    return paths;
}
}
