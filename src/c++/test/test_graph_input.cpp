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