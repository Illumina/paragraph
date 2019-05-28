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

#include "grm/GraphAligner.hh"

#include <cstdlib>
#include <gssw.h>
#include <iostream>
#include <list>
#include <sstream>
#include <string>

#include "common/Error.hh"
#include "common/StringUtil.hh"
#include "graphcore/GraphOperations.hh"
#include "graphutils/SequenceOperations.hh"

using namespace grm;
using namespace common;
using std::string;

using graphtools::Graph;
using graphtools::NodeId;
using graphtools::reverseComplement;
using graphtools::reverseGraph;
using graphtools::reverseString;

extern "C" {
#include "gssw.h"
}

static inline void safe_gssw_graph_destroy(gssw_graph* g)
{
    if (g)
    {
        gssw_graph_destroy(g);
    }
}

static inline void safe_gssw_graph_mapping_destroy(gssw_graph_mapping* g)
{
    if (g)
    {
        gssw_graph_mapping_destroy(g);
    }
}

struct GraphAligner::GraphAlignerImpl
{
    typedef std::unique_ptr<int8_t, decltype(&free)> p_int8_t;
    typedef std::unique_ptr<gssw_graph, decltype(&safe_gssw_graph_destroy)> p_gssw_graph;
    typedef std::unique_ptr<gssw_graph_mapping, decltype(&safe_gssw_graph_mapping_destroy)> p_gssw_graph_mapping;

    GraphAlignerImpl()
        : nt_table_(p_int8_t(gssw_create_nt_table(), free))
        , mat_(p_int8_t(gssw_create_score_matrix(match_, mismatch_), free))
        ,
        /* graph_ owns all nodes_[] elements. There probably is a better way to keep track of this */
        graph_(nullptr, safe_gssw_graph_destroy)
        , graph_reversed_(nullptr, safe_gssw_graph_destroy)
        , original_graph_(nullptr)
    {
    }

    /**
     * Helper function that extracts cigar strings from
     * gssw_graph_mapping structs; largely copied from GSSW.c and
     * modified to print cigar strings to std::string instead of FILE*.
     * @param gm  graph mapping to extract from
     * @return cigar string
     */
    static std::string extractCigar(gssw_graph_mapping* gm, std::vector<NodeId> const& node_map)
    {
        std::stringstream cigar;
        // cigar << gm->score << "@" << gm->position;
        gssw_graph_cigar* g = &gm->cigar;
        gssw_node_cigar* nc = g->elements;

        for (uint32_t i = 0; i < g->length; ++i, ++nc)
        {
            cigar << node_map[nc->node->id] << "[";
            gssw_cigar* c = nc->cigar;
            gssw_cigar_element* e = c->elements;
            for (int32_t j = 0; j < c->length; ++j, ++e)
            {
                cigar << e->length << e->type;
            }
            cigar << "]";
        }

        return cigar.str();
    }

    void initializeGraph(
        Graph const& graph, p_gssw_graph& gssw_graph, std::vector<gssw_node*>& nodes, std::vector<NodeId>& node_map,
        std::vector<uint32_t>& first_gssw_node)
    {
        nodes.clear();
        node_map.clear();
        first_gssw_node.clear();
        first_gssw_node.resize(graph.numNodes());

        uint32_t gssw_node_id = 0;
        for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
        {
            first_gssw_node[node_id] = gssw_node_id;
            if (node_id != 0 && node_id != graph.numNodes() - 1)
            {
                for (auto sequence : graph.nodeSeqExpansion(node_id))
                {
                    stringutil::toUpper(sequence);

                    gssw_node* nd
                        = gssw_node_create(nullptr, gssw_node_id++, sequence.c_str(), nt_table_.get(), mat_.get());
                    nodes.push_back(nd);
                    node_map.push_back(node_id);
                }
            }
            else
            {
                std::string sequence = graph.nodeSeq(node_id);
                stringutil::toUpper(sequence);
                gssw_node* nd
                    = gssw_node_create(nullptr, gssw_node_id++, sequence.c_str(), nt_table_.get(), mat_.get());
                nodes.push_back(nd);
                node_map.push_back(node_id);
            }

            for (auto to_gssw_node = first_gssw_node[node_id]; to_gssw_node < gssw_node_id; ++to_gssw_node)
            {
                for (auto pred : graph.predecessors(node_id))
                {
                    // topological ordering
                    assert(pred < node_id);
                    auto from_gssw_node = first_gssw_node[pred];
                    while (from_gssw_node < node_map.size() && node_map[from_gssw_node] == pred)
                    {
                        gssw_nodes_add_edge(nodes[from_gssw_node], nodes[to_gssw_node]);
                        ++from_gssw_node;
                    }
                }
            }
        }

        // Make the graph itself.
        gssw_graph.reset(gssw_graph_create((uint32_t)nodes.size()));
        for (gssw_node* node : nodes)
        {
            gssw_graph_add_node(gssw_graph.get(), node);
        }
    }

    /** Assumes that DP matrix has been filled out. */
    static bool alignsEndAtMultNodes(
        gssw_graph* graph, std::vector<gssw_node*> const& nodes, std::vector<NodeId> const& node_map, int32_t read_len)
    {
        /** Top local-alignment score. */
        uint16_t top_score = graph->max_node->alignment->score1;

        std::set<NodeId> nodes_where_top_aligns_end;

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            auto* mH = (uint8_t*)nodes[i]->alignment->mH;
            bool top_align_ends_at_node = false;
            /** Rows correspond to reference and columns to read bases. */
            for (int32_t row = 0; row < nodes[i]->len; ++row)
            {
                for (int32_t col = 0; col < read_len; ++col)
                {
                    uint8_t score = mH[read_len * row + col];
                    if (score == top_score)
                    {
                        top_align_ends_at_node = true;
                        break;
                    }
                }
                if (top_align_ends_at_node)
                {
                    break;
                }
            }

            if (top_align_ends_at_node)
            {
                nodes_where_top_aligns_end.insert(node_map[i]);
            }

            if (nodes_where_top_aligns_end.size() > 1)
            {
                return true;
            }
        }

        return false;
    }

    bool alignString(
        gssw_graph* g, std::vector<gssw_node*> const& nodes, std::vector<NodeId> const& node_map, std::string str,
        p_gssw_graph_mapping& output)
    {
        stringutil::toUpper(str);
        // Compute dynamic programming matrix.
        gssw_graph_fill(g, str.c_str(), nt_table_.get(), mat_.get(), gap_open_, gap_extension_, 15, 2);

        // Compute optimal local alignment.
        gssw_graph_mapping* gm = gssw_graph_trace_back(
            g, str.c_str(), (int32_t)str.length(), nt_table_.get(), mat_.get(), gap_open_, gap_extension_);
        output = p_gssw_graph_mapping(gm, gssw_graph_mapping_destroy);
        return alignsEndAtMultNodes(g, nodes, node_map, (int32_t)str.size());
    }

    /** Default alignment parameters. */
    int8_t match_ = 1;
    int8_t mismatch_ = 4;
    uint8_t gap_open_ = 6;
    uint8_t gap_extension_ = 1;

    p_int8_t nt_table_;
    p_int8_t mat_;

    p_gssw_graph graph_;
    /** nodes are owned by graph */
    std::vector<gssw_node*> nodes_;

    /** each node in the input graph may create more than one gssw node */
    std::vector<NodeId> node_map_;
    /** stores the first gssw node for each graph node */
    std::vector<uint32_t> first_gssw_node_;

    p_gssw_graph graph_reversed_;
    /** nodes are owned by graph */
    std::vector<gssw_node*> nodes_reversed_;
    /** each node in the input graph may create more than one gssw node */
    std::vector<NodeId> node_map_reversed_;
    /** stores the first gssw node for each graph node */
    std::vector<uint32_t> first_gssw_node_reversed_;

    /** keep original graphtools graph */
    Graph const* original_graph_;
};

GraphAligner::GraphAligner()
    : _impl(new GraphAlignerImpl())
{
}

GraphAligner::~GraphAligner() = default;

GraphAligner::GraphAligner(GraphAligner&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

GraphAligner& GraphAligner::operator=(GraphAligner&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

void GraphAligner::setGraph(Graph const* graph)
{
    _impl->original_graph_ = graph;
    _impl->initializeGraph(*graph, _impl->graph_, _impl->nodes_, _impl->node_map_, _impl->first_gssw_node_);
    Graph graph_reverse = reverseGraph(*graph);
    _impl->initializeGraph(
        graph_reverse, _impl->graph_reversed_, _impl->nodes_reversed_, _impl->node_map_reversed_,
        _impl->first_gssw_node_reversed_);
}

string GraphAligner::align(const string& read, int& mapq, int& position, int& score) const
{
    Read temp_read;
    temp_read.set_bases(read);
    alignRead(temp_read, AF_CIGAR);
    mapq = temp_read.graph_mapq();
    position = temp_read.graph_pos();
    score = temp_read.graph_alignment_score();
    return temp_read.graph_cigar();
}

/**
 * Align a read to the graph and update the graph_* fields.
 *
 * We will align both the forward and the reverse strand and return
 * the alignment which is unique, or with the better score, or default
 * to the forward strand.
 *
 * @param read read structure
 * @param alignment_flags flags for alignment (see above)
 */
void GraphAligner::alignRead(Read& read, unsigned int alignment_flags) const
{
    GraphAlignerImpl::p_gssw_graph_mapping gm_fwd_strand{ nullptr, safe_gssw_graph_mapping_destroy };
    GraphAlignerImpl::p_gssw_graph_mapping gm_reverse_strand{ nullptr, safe_gssw_graph_mapping_destroy };
    GraphAlignerImpl::p_gssw_graph_mapping rgm_fwd_strand{ nullptr, safe_gssw_graph_mapping_destroy };
    GraphAlignerImpl::p_gssw_graph_mapping rgm_reverse_strand{ nullptr, safe_gssw_graph_mapping_destroy };

    const std::string rev_cmp_bases = reverseComplement(read.bases());

    const bool fwd_strand_multi_align
        = _impl->alignString(_impl->graph_.get(), _impl->nodes_, _impl->node_map_, read.bases(), gm_fwd_strand);
    const bool reverse_strand_multi_align = (alignment_flags & AF_BOTH_STRANDS) != 0u
        ? _impl->alignString(_impl->graph_.get(), _impl->nodes_, _impl->node_map_, rev_cmp_bases, gm_reverse_strand)
        : false;

    bool rfwd_strand_multi_align = false;
    bool rreverse_strand_multi_align = false;

    if ((alignment_flags & AF_REVERSE_GRAPH) != 0u)
    {
        string bases_rev = read.bases();
        std::reverse(bases_rev.begin(), bases_rev.end());

        rfwd_strand_multi_align = _impl->alignString(
            _impl->graph_reversed_.get(), _impl->nodes_reversed_, _impl->node_map_reversed_, bases_rev, rgm_fwd_strand);
        rreverse_strand_multi_align = (alignment_flags & AF_BOTH_STRANDS) != 0u
            ? _impl->alignString(
                  _impl->graph_reversed_.get(), _impl->nodes_reversed_, _impl->node_map_reversed_,
                  reverseComplement(bases_rev), rgm_reverse_strand)
            : false;
    }

    bool fwd_strand_is_unique = (!fwd_strand_multi_align) && (!rfwd_strand_multi_align);
    bool reverse_strand_is_unique = (!reverse_strand_multi_align) && (!rreverse_strand_multi_align);

    // prefer unique alignment. if not unique, or if both unique, prefer alignment with better score
    bool return_reverse = false;
    if ((!fwd_strand_is_unique) && reverse_strand_is_unique && gm_reverse_strand)
    {
        return_reverse = true;
    }
    else if (fwd_strand_is_unique && (!reverse_strand_is_unique))
    {
        return_reverse = false;
    }
    else if (gm_reverse_strand)
    {
        return_reverse = gm_fwd_strand->score < gm_reverse_strand->score;
    }

    const bool resulting_strand_is_reverse = read.is_reverse_strand() != return_reverse;
    read.set_is_graph_reverse_strand(resulting_strand_is_reverse);

    auto make_node_list = [this](gssw_graph_mapping* gm) -> std::list<int> {
        std::list<int> result;
        auto const& g = gm->cigar;
        auto nc = g.elements;
        for (uint32_t i = 0; i < g.length; ++i, ++nc)
        {
            result.push_back(_impl->node_map_[nc->node->id]);
        }
        return result;
    };

    std::list<int> nodes_passed;
    if (return_reverse)
    {
        read.set_bases(rev_cmp_bases);
        std::string rev_quals = read.quals();
        std::reverse(rev_quals.begin(), rev_quals.end());
        read.set_quals(rev_quals);

        read.set_graph_pos(gm_reverse_strand->position);
        read.set_graph_alignment_score(gm_reverse_strand->score);
        read.set_is_graph_alignment_unique(reverse_strand_is_unique);
        read.set_graph_mapq(reverse_strand_is_unique ? 60 : 0);
        if (alignment_flags & AF_CIGAR)
        {
            read.set_graph_cigar(_impl->extractCigar(gm_reverse_strand.get(), _impl->node_map_));
        }
        nodes_passed = make_node_list(gm_reverse_strand.get());
    }
    else
    {
        read.set_graph_pos(gm_fwd_strand->position);
        read.set_graph_alignment_score(gm_fwd_strand->score);
        read.set_is_graph_alignment_unique(fwd_strand_is_unique);
        read.set_graph_mapq(fwd_strand_is_unique ? 60 : 0);
        if (alignment_flags & AF_CIGAR)
        {
            read.set_graph_cigar(_impl->extractCigar(gm_fwd_strand.get(), _impl->node_map_));
        }
        nodes_passed = make_node_list(gm_fwd_strand.get());
    }

    LOG()->trace("GraphAligner::alignRead read.cigar={}", read.graph_cigar());
}
