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

#include "grm/GraphInput.hh"

#include "common.hh"
#include "gtest/gtest.h"
#include "json/json.h"
#include <fstream>

#include "common/Error.hh"

using namespace grm;
using namespace graphtools;
using std::string;

Json::Value getJsonRoot(const string& gspath)
{
    Json::Value root;
    std::ifstream graph_spec(gspath);
    assert(graph_spec.good());

    Json::CharReaderBuilder rBuilder;
    std::string errors;
    bool success = Json::parseFromStream(rBuilder, graph_spec, &root, &errors);
    assert(success);
    return root;
}

TEST(Graph, LoadsGraphWithEdgesAndNodes)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-nodes.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    Graph test_graph = graphFromJson(root, reference_path);

    EXPECT_EQ(test_graph.numEdges(), (size_t)5);
    EXPECT_EQ(test_graph.numNodes(), (size_t)5);

    for (NodeId node_id = 0; node_id < test_graph.numNodes(); ++node_id)
    {
        EXPECT_GT(test_graph.nodeSeq(node_id).size(), (size_t)0);
    }
}

TEST(Graph, LoadsGraphWithEdgesAndNodesNoRefSeq)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-nodes.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    ;
    Json::Value root = getJsonRoot(graph_spec_path);
    Graph test_graph = graphFromJson(root, reference_path, false);
    for (NodeId node_id = 0; node_id < test_graph.numNodes(); ++node_id)
    {
        if (test_graph.nodeName(node_id) != "source" && test_graph.nodeName(node_id) != "sink")
        {
            EXPECT_EQ(test_graph.nodeSeq(node_id).size(), (size_t)0);
        }
    }
}

TEST(Graph, LoadsGraphWithNodesOnly)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-nodes-only.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    Graph test_graph = graphFromJson(root, reference_path);

    EXPECT_EQ(test_graph.numEdges(), (size_t)0);
    EXPECT_EQ(test_graph.numNodes(), (size_t)3);
}

TEST(Graph, LoadsGraphWithMultiNode)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-ref-node-array.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    Graph test_graph = graphFromJson(root, reference_path);

    EXPECT_EQ(test_graph.numEdges(), (size_t)0);
    EXPECT_EQ(test_graph.numNodes(), (size_t)4);

    for (NodeId node_id = 0; node_id < test_graph.numNodes(); ++node_id)
    {
        EXPECT_GT(test_graph.nodeSeq(node_id).size(), (size_t)0);
    }
}

TEST(Graph, NoReferenceOrSequenceNodeIdDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-no-ref-or-seq-node-key.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    EXPECT_THROW(graphFromJson(root, reference_path), std::exception);
}

TEST(Graph, EdgesOnlyDeathTest)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-edges-only.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    EXPECT_THROW(graphFromJson(root, reference_path), std::exception);
}

TEST(Graph, BadEdgesValueDeathTest)
{
    const string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-bad-edges-value.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    EXPECT_THROW(graphFromJson(root, reference_path), std::exception);
}

TEST(Graph, BadNodeSequenceIdsDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-bad-node-seq-ids.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    EXPECT_THROW(graphFromJson(root, reference_path), std::exception);
}

TEST(Graph, DuplicateNodeNameDeathTest)
{
    const string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/basic/del-with-duplicate-node-names.json";
    const string reference_path = g_testenv->getBasePath() + "/../share/test-data/basic/dummy.fa";

    Json::Value root = getJsonRoot(graph_spec_path);
    EXPECT_THROW(graphFromJson(root, reference_path), std::exception);
}