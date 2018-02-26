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

#include "common/Fasta.hh"
#include "common/Region.hh"
#include "common/StringUtil.hh"
#include "graphs/Graph.hh"

#include "common/Error.hh"

using std::string;

namespace graphs
{

namespace graphs_internal
{

    /**
     * Clean up protobuf so we don't cause memory leaks
     */
    struct PB_Cleanup
    {
        ~PB_Cleanup() { ::google::protobuf::ShutdownProtobufLibrary(); }
    };

    volatile PB_Cleanup cleaner;
}

static void setReferenceNode(
    p_Node& node_inst, string const& reference_location, bool store_ref_sequence, common::FastaFile const& ref)
{
    if (store_ref_sequence)
    {
        node_inst->set_sequence(ref.query(reference_location));
    }

    auto ref_pos = node_inst->add_reference_pos();
    string chr;
    int64_t start = -1;
    int64_t end = -1;
    common::stringutil::parsePos(reference_location, chr, start, end);
    ref_pos->set_chrom(chr);
    ref_pos->set_pos(static_cast<google::protobuf::uint64>(start));
    ref_pos->set_end(static_cast<google::protobuf::uint64>(end));
}

/**
 * Initialize graph from JSON
 * @param in Input JSON node
 * @param reference reference fasta name
 * @param out Output graph to initialize
 */
void fromJson(Json::Value const& in, string const& reference, Graph& out, bool store_ref_sequence)
{
    common::FastaFile ref(reference);

    p_Header h{ new graphs::GraphHeader };
    out.header = std::move(h);
    out.nodes.clear();
    out.edges.clear();

    Json::Value const* in_graph = &in;
    if (in.isMember("graph"))
    {
        in_graph = &in["graph"];
    }

    std::map<string, size_t> sequencenames;
    if ((*in_graph).isMember("sequencenames"))
    {
        assert((*in_graph)["sequencenames"].type() == Json::ValueType::arrayValue);
        for (size_t i = 0; i < (*in_graph)["sequencenames"].size(); ++i)
        {
            auto const& sequencename = (*in_graph)["sequencenames"][(int)i].asString();
            assert(sequencenames.count(sequencename) == 0);
            out.header->add_sequencenames(sequencename);
            sequencenames[sequencename] = i;
        }
    }

    assert((*in_graph)["nodes"].type() == Json::ValueType::arrayValue);

    if ((*in_graph)["edges"].type() != Json::ValueType::nullValue)
    {
        assert((*in_graph)["edges"].type() == Json::ValueType::arrayValue);
    }

    std::map<string, uint64_t> node_map;
    for (size_t i = 0; i < (*in_graph)["nodes"].size(); ++i)
    {
        auto const& in_n = (*in_graph)["nodes"][(int)i];
        p_Node out_n{ new Node };
        for (auto& sequence : in_n["sequences"])
        {
            assert(sequencenames.find(sequence.asString()) != sequencenames.end());
            out_n->add_sequence_ids(sequencenames[sequence.asString()]);
        }

        const string name = in_n.isMember("name") ? in_n["name"].asString() : string("node-") + std::to_string(i + 1);

        assert(!node_map.count(name));

        out_n->set_name(name);

        assert(in_n.isMember("sequence") || in_n.isMember("reference"));

        if (in_n.isMember("sequence"))
        {
            out_n->set_sequence(in_n["sequence"].asString());
        }
        else
        {
            if (in_n["reference"].type() == Json::ValueType::stringValue)
            {
                const string reference_location = in_n["reference"].asString();
                setReferenceNode(out_n, reference_location, store_ref_sequence, ref);
            }
            else
            {
                assert(in_n["reference"].type() == Json::ValueType::arrayValue);

                string prev_seq;
                for (const auto& ref_inst : in_n["reference"])
                {
                    assert(ref_inst.type() == Json::ValueType::stringValue);

                    const string reference_location = ref_inst.asString();
                    const string curr_seq(ref.query(reference_location));

                    if (!prev_seq.empty())
                    {
                        assert(prev_seq == curr_seq);
                    }

                    setReferenceNode(out_n, reference_location, store_ref_sequence, ref);
                    prev_seq = curr_seq;
                }
            }
        }

        out.nodes.insert(std::pair<uint64_t, p_Node>(i, out_n));
        node_map[name] = i;
    }

    for (size_t i = 0; i < (*in_graph)["edges"].size(); ++i)
    {
        auto const& in_e = (*in_graph)["edges"][(int)i];
        p_Edge out_e{ new Edge };

        assert(node_map.count(in_e["from"].asString()));
        assert(node_map.count(in_e["to"].asString()));
        out_e->set_from(node_map[in_e["from"].asString()]);
        out_e->set_to(node_map[in_e["to"].asString()]);
        for (auto& sequence : in_e["sequences"])
        {
            assert(sequencenames.find(sequence.asString()) != sequencenames.end());
            out_e->add_sequence_ids(sequencenames[sequence.asString()]);
        }
        out.edges.push_back(out_e);
    }
}

void reverse(const Graph& g, Graph& g_rev)
{
    // Make header.
    p_Header h{ new graphs::GraphHeader };
    g_rev.header = std::move(h);
    g_rev.nodes.clear();
    g_rev.edges.clear();

    const uint64_t num_nodes = g.nodes.size();

    // Make nodes.
    for (const auto& id_node : g.nodes)
    {
        p_Node node_rev{ new graphs::Node };
        string node_rev_seq = id_node.second->sequence();
        std::reverse(node_rev_seq.begin(), node_rev_seq.end());
        node_rev->set_sequence(node_rev_seq);
        const uint64_t node_rev_id = num_nodes - id_node.first - 1;
        node_rev->set_id(node_rev_id);
        g_rev.nodes[node_rev_id] = node_rev;
    }

    // Make edges.
    for (const auto& e : g.edges)
    {
        p_Edge e_rev{ new graphs::Edge };
        e_rev->set_from(num_nodes - e->to() - 1);
        e_rev->set_to(num_nodes - e->from() - 1);
        g_rev.edges.push_back(e_rev);
    }
}
}
