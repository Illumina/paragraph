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

#include "graphs/Graph.hh"

#include "common.hh"
#include "gtest/gtest.h"
#include "json/json.h"
#include <fstream>

#include "common/Error.hh"

using namespace graphs;
using std::string;

Json::Value getJsonRoot(const string& gspath)
{
    Json::Value root;
    Json::Reader reader;
    std::ifstream graph_spec(gspath);

    assert(graph_spec.good());

    bool success = reader.parse(graph_spec, root);

    assert(success);
    return root;
}

TEST(Graph, LoadsGraphWithEdgesAndNodes)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-nodes.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);
    graphs::fromJson(root, reference_path, test_graph);

    EXPECT_EQ(test_graph.edges.size(), (size_t)5);
    EXPECT_EQ(test_graph.nodes.size(), (size_t)5);
    EXPECT_EQ(test_graph.header->sequencenames().size(), 2);

    for (const auto& nit : test_graph.nodes)
    {
        EXPECT_GT(nit.second->sequence().size(), (size_t)0);
    }
}

TEST(Graph, LoadsGraphWithEdgesAndNodesNoRefSeq)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-nodes.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);
    bool no_store_ref_sequence = false;
    graphs::fromJson(root, reference_path, test_graph, no_store_ref_sequence);
    for (const auto& nit : test_graph.nodes)
    {
        if (nit.second->name() != "source" && nit.second->name() != "sink")
        {
            EXPECT_EQ(nit.second->sequence().size(), (size_t)0);
        }
    }
}

TEST(Graph, LoadsGraphWithNodesOnly)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-nodes-only.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);
    graphs::fromJson(root, reference_path, test_graph);

    EXPECT_EQ(test_graph.edges.size(), (size_t)0);
    EXPECT_EQ(test_graph.nodes.size(), (size_t)3);
    EXPECT_EQ(test_graph.header->sequencenames().size(), 1);
}

TEST(Graph, LoadsGraphWithMultiNode)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-ref-node-array.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);
    graphs::fromJson(root, reference_path, test_graph);

    EXPECT_EQ(test_graph.edges.size(), (size_t)0);
    EXPECT_EQ(test_graph.nodes.size(), (size_t)4);
    EXPECT_EQ(test_graph.header->sequencenames().size(), 2);

    for (const auto& nit : test_graph.nodes)
    {
        if (nit.second->reference_pos_size() > 1)
        {
            EXPECT_EQ(nit.second->reference_pos_size(), 3);
        }
    }
}

TEST(Graph, NoReferenceOrSequenceNodeIdDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-no-ref-or-seq-node-key.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}

TEST(Graph, DuplicatedSeqsDeathTest)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-duplicate-seqs.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}

TEST(Graph, EdgesOnlyDeathTest)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-only.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}

TEST(Graph, BadEdgesValueDeathTest)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-bad-edges-value.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}

TEST(Graph, BadNodeSequenceIdsDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-bad-node-seq-ids.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}

TEST(Graph, DuplicateNodeNameDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-duplicate-node-names.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Graph test_graph;
    Json::Value root = getJsonRoot(graph_spec_path);

    EXPECT_THROW(graphs::fromJson(root, reference_path, test_graph), std::exception);
}