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

#include "graphs/GraphBuilders.hh"

#include <cmath>

using std::ceil;
using std::list;
using std::map;
using std::string;
using std::vector;

namespace graphs
{
void makeDeletionGraph(const string& left_flank, const string& deletion, const string& right_flank, Graph& g)
{
    // Make header
    p_Header h{ new graphs::GraphHeader };
    g.header = std::move(h);
    g.nodes.clear();
    g.edges.clear();

    // Make nodes.
    p_Node n1{ new graphs::Node };
    n1->set_id(0);
    n1->set_sequence(left_flank.c_str());
    n1->set_name("LF");
    g.nodes[0] = n1;

    p_Node n2{ new graphs::Node };
    n2->set_id(1);
    n2->set_sequence(deletion.c_str());
    n2->set_name("DEL");
    g.nodes[1] = n2;

    p_Node n3{ new graphs::Node };
    n3->set_id(2);
    n3->set_sequence(right_flank.c_str());
    n3->set_name("RF");
    g.nodes[2] = n3;

    // Make edges.
    p_Edge e1{ new graphs::Edge };
    e1->set_from(0);
    e1->set_to(1);
    g.edges.push_back(e1);

    p_Edge e2{ new graphs::Edge };
    e2->set_from(0);
    e2->set_to(2);
    g.edges.push_back(e2);

    p_Edge e3{ new graphs::Edge };
    e3->set_from(1);
    e3->set_to(2);
    g.edges.push_back(e3);
}

void makeSimpleSwapGraph(
    const std::string& left_flank, const std::string& deletion, const std::string& insertion,
    const std::string& right_flank, Graph& g)
{
    // Make header
    p_Header h{ new graphs::GraphHeader };
    g.header = std::move(h);
    g.nodes.clear();
    g.edges.clear();

    // Make nodes.
    p_Node lf_node{ new graphs::Node };
    lf_node->set_id(0);
    lf_node->set_sequence(left_flank.c_str());
    lf_node->set_name("LF");
    g.nodes[0] = lf_node;

    p_Node del_node{ new graphs::Node };
    del_node->set_id(1);
    del_node->set_sequence(deletion.c_str());
    del_node->set_name("DEL");
    g.nodes[1] = del_node;

    p_Node ins_node{ new graphs::Node };
    ins_node->set_id(2);
    ins_node->set_sequence(insertion.c_str());
    ins_node->set_name("INS");
    g.nodes[2] = ins_node;

    p_Node rf_node{ new graphs::Node };
    rf_node->set_id(3);
    rf_node->set_sequence(right_flank.c_str());
    rf_node->set_name("RF");
    g.nodes[3] = rf_node;

    // Make edges.
    p_Edge edge_0{ new graphs::Edge };
    edge_0->set_from(0);
    edge_0->set_to(1);
    g.edges.push_back(edge_0);

    p_Edge edge_1{ new graphs::Edge };
    edge_1->set_from(0);
    edge_1->set_to(2);
    g.edges.push_back(edge_1);

    p_Edge edge_2{ new graphs::Edge };
    edge_2->set_from(1);
    edge_2->set_to(3);
    g.edges.push_back(edge_2);

    p_Edge edge_3{ new graphs::Edge };
    edge_3->set_from(2);
    edge_3->set_to(3);
    g.edges.push_back(edge_3);
}

void makeSwapGraph(
    const string& left_flank, const string& deletion, const list<string>& alternatives, const string& right_flank,
    Graph& g)
{
    // Make header
    p_Header h{ new graphs::GraphHeader };
    g.header = std::move(h);
    g.nodes.clear();
    g.edges.clear();

    // Make nodes.
    p_Node n1{ new graphs::Node };
    n1->set_id(0);
    n1->set_sequence(left_flank.c_str());
    n1->set_name("LF");
    g.nodes[0] = n1;

    const uint64_t min_alt_id = (deletion.empty() ? 1 : 2);
    const uint64_t rf_id = (uint64_t)(min_alt_id + alternatives.size());
    if (!deletion.empty())
    {
        p_Node n2{ new graphs::Node };
        n2->set_id(1);
        n2->set_sequence(deletion.c_str());
        n2->set_name("DEL");
        g.nodes[1] = n2;

        // LF -> DEL
        p_Edge e1{ new graphs::Edge };
        e1->set_from(0);
        e1->set_to(1);
        g.edges.push_back(e1);

        // DEL -> RF
        p_Edge e3{ new graphs::Edge };
        e3->set_from(1);
        e3->set_to(rf_id);
        g.edges.push_back(e3);
    }

    p_Node nx{ new graphs::Node };
    nx->set_id(rf_id);
    nx->set_sequence(right_flank.c_str());
    nx->set_name("RF");
    g.nodes[rf_id] = nx;

    // LF -> RF
    p_Edge e2{ new graphs::Edge };
    e2->set_from(0);
    e2->set_to(rf_id);
    g.edges.push_back(e2);

    uint64_t alt_id = min_alt_id;
    for (const auto& a : alternatives)
    {
        p_Node na{ new graphs::Node };
        na->set_id(alt_id);
        na->set_sequence(a.c_str());
        if (alt_id == min_alt_id && alternatives.size() == 1)
        {
            na->set_name("INS");
        }
        else
        {
            na->set_name(string("ALT") + std::to_string(alt_id - 1));
        }
        g.nodes[alt_id] = na;

        // LF -> ALTX/INS
        p_Edge e_in{ new graphs::Edge };
        e_in->set_from(0);
        e_in->set_to(alt_id);
        g.edges.push_back(e_in);

        // ALTX/INS -> RF
        p_Edge e_out{ new graphs::Edge };
        e_out->set_from(alt_id);
        e_out->set_to(rf_id);
        g.edges.push_back(e_out);
        alt_id++;
    }
}

void makeDoubleSwapGraph(
    const std::string& left_flank, const std::string& deletion1, const std::string& insertion1,
    const std::string& middle, const std::string& deletion2, const std::string& insertion2,
    const std::string& right_flank, Graph& g)
{
    // Make header
    p_Header h{ new graphs::GraphHeader };
    g.header = std::move(h);
    g.nodes.clear();
    g.edges.clear();

    // Make nodes.
    p_Node lf_node{ new graphs::Node };
    lf_node->set_id(0);
    lf_node->set_sequence(left_flank.c_str());
    lf_node->set_name("LF");
    g.nodes[0] = lf_node;

    p_Node del1_node{ new graphs::Node };
    del1_node->set_id(1);
    del1_node->set_sequence(deletion1.c_str());
    del1_node->set_name("DEL1");
    g.nodes[1] = del1_node;

    p_Node ins1_node{ new graphs::Node };
    ins1_node->set_id(2);
    ins1_node->set_sequence(insertion1.c_str());
    ins1_node->set_name("INS1");
    g.nodes[2] = ins1_node;

    p_Node mid_node{ new graphs::Node };
    mid_node->set_id(3);
    mid_node->set_sequence(middle.c_str());
    mid_node->set_name("MID");
    g.nodes[3] = mid_node;

    p_Node del2_node{ new graphs::Node };
    del2_node->set_id(4);
    del2_node->set_sequence(deletion2.c_str());
    del2_node->set_name("DEL2");
    g.nodes[4] = del2_node;

    p_Node ins2_node{ new graphs::Node };
    ins2_node->set_id(5);
    ins2_node->set_sequence(insertion2.c_str());
    ins2_node->set_name("INS2");
    g.nodes[5] = ins2_node;

    p_Node rf_node{ new graphs::Node };
    rf_node->set_id(6);
    rf_node->set_sequence(right_flank.c_str());
    rf_node->set_name("RF");
    g.nodes[6] = rf_node;

    // Make edges.
    p_Edge edge_0_1{ new graphs::Edge };
    edge_0_1->set_from(0);
    edge_0_1->set_to(1);
    g.edges.push_back(edge_0_1);

    p_Edge edge_0_2{ new graphs::Edge };
    edge_0_2->set_from(0);
    edge_0_2->set_to(2);
    g.edges.push_back(edge_0_2);

    p_Edge edge_1_3{ new graphs::Edge };
    edge_1_3->set_from(1);
    edge_1_3->set_to(3);
    g.edges.push_back(edge_1_3);

    p_Edge edge_2_3{ new graphs::Edge };
    edge_2_3->set_from(2);
    edge_2_3->set_to(3);
    g.edges.push_back(edge_2_3);

    p_Edge edge_3_4{ new graphs::Edge };
    edge_3_4->set_from(3);
    edge_3_4->set_to(4);
    g.edges.push_back(edge_3_4);

    p_Edge edge_3_5{ new graphs::Edge };
    edge_3_5->set_from(3);
    edge_3_5->set_to(5);
    g.edges.push_back(edge_3_5);

    p_Edge edge_4_6{ new graphs::Edge };
    edge_4_6->set_from(4);
    edge_4_6->set_to(6);
    g.edges.push_back(edge_4_6);

    p_Edge edge_5_6{ new graphs::Edge };
    edge_5_6->set_from(5);
    edge_5_6->set_to(6);
    g.edges.push_back(edge_5_6);
}

std::map<StrSpecUnit::UnitType, const char*> StrSpecUnit::unit_type_to_string
    = { { StrSpecUnit::UnitType::STR, "STR" }, { StrSpecUnit::UnitType::SEQ, "SEQ" } };

/**
 * Create a multi-unit STR graph.
 *
 * @param left_flank left flanking reference sequence
 * @param right_flank right flanking reference sequence
 * @param g output graph structure
 */
void makeStrGraph(
    const string& left_flank, const string& right_flank, size_t read_len, std::vector<StrSpecUnit> spec, Graph& g)
{
    StrSpecUnit unit_lf, unit_rf;
    unit_lf.seq = left_flank;
    unit_lf.type = StrSpecUnit::UnitType::SEQ;
    unit_rf.seq = right_flank;
    unit_rf.type = StrSpecUnit::UnitType::SEQ;

    vector<StrSpecUnit> str_spec = { unit_lf };
    str_spec.insert(str_spec.end(), spec.begin(), spec.end());
    str_spec.push_back(unit_rf);

    size_t num_nodes = 0; // Flanks are always present.
    map<string, size_t> num_nodes_per_str_loop;

    for (const auto& node_spec : str_spec)
    {
        if (node_spec.type == StrSpecUnit::UnitType::STR)
        {
            num_nodes_per_str_loop[node_spec.seq] = (size_t)ceil(read_len / (double)node_spec.seq.size());
            num_nodes += num_nodes_per_str_loop[node_spec.seq];
        }
        else
        {
            num_nodes += 1;
        }
    }

    p_Header h{ new graphs::GraphHeader };
    g.header = std::move(h);
    g.nodes.clear();
    g.edges.clear();

    size_t node_index = 0;
    for (size_t unit_index = 0; unit_index != str_spec.size() - 1; ++unit_index)
    {
        const string& cur_unit_seq = str_spec[unit_index].seq;
        const bool is_str_unit = str_spec[unit_index].type == StrSpecUnit::UnitType::STR;
        const size_t num_loop_nodes = is_str_unit ? num_nodes_per_str_loop[cur_unit_seq] : 0;
        const size_t next_unit_node_index = is_str_unit ? node_index + num_loop_nodes : node_index + 1;

        string node_name = "LF";
        if (unit_index != 0)
        {
            node_name = is_str_unit ? "U" : "S";
            node_name += std::to_string(unit_index);
        }

        p_Node cur_node{ new graphs::Node };
        cur_node->set_id(node_index);
        cur_node->set_sequence(cur_unit_seq.c_str());
        cur_node->set_name(node_name.c_str());
        g.nodes[node_index] = cur_node;

        p_Edge cur_edge{ new graphs::Edge };
        cur_edge->set_from(node_index);
        cur_edge->set_to(next_unit_node_index);
        g.edges.push_back(cur_edge);

        ++node_index;

        for (size_t loop_index = 1; loop_index < num_loop_nodes; ++loop_index)
        {
            p_Node loop_node{ new graphs::Node };
            loop_node->set_id(node_index);
            loop_node->set_sequence(cur_unit_seq.c_str());
            loop_node->set_name(node_name.c_str());
            g.nodes[node_index] = loop_node;

            p_Edge loop_edge_to{ new graphs::Edge };
            loop_edge_to->set_from(node_index - 1);
            loop_edge_to->set_to(node_index);
            g.edges.push_back(loop_edge_to);

            p_Edge loop_edge_from{ new graphs::Edge };
            loop_edge_from->set_from(node_index);
            loop_edge_from->set_to(next_unit_node_index);
            g.edges.push_back(loop_edge_from);

            ++node_index;
        }
    }

    p_Node rf_node{ new graphs::Node };
    rf_node->set_id(node_index);
    rf_node->set_sequence(right_flank.c_str());
    rf_node->set_name("RF");
    g.nodes[node_index] = rf_node;
}
}