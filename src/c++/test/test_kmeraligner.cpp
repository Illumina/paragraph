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

#include "common/ReadReader.hh"
#include "graphs/Graph.hh"
#include "grm/KmerAligner.hh"
#include "paragraph/Disambiguation.hh"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;

using namespace testing;
using namespace common;
using namespace graphs;

class KmerAlignerTest : public Test
{
public:
    vector<Read> reads{ 6 };

    WalkableGraph graph;

    void SetUp() override
    {
        reads[0].setCoreInfo("f1", "AAAAAAAATTTTTTTTAAAAAAAA", "########################");
        reads[1].setCoreInfo("f2", "TTTTTTAAAAAAAATTTTTTT", "#####################");
        reads[2].setCoreInfo("f3", "AAAAAGGGGGGGGAAAAAA", "###################");
        reads[3].setCoreInfo("f4", "AAAAGGGGGGGGAAAAAA", "##################");
        reads[4].setCoreInfo("f5", "TTTTTTCCCCCCCCTTTTT", "###################");
        // this read won't align uniquely here since we can move the A matches around
        reads[5].setCoreInfo("f6", "AAAAAAAAAAAAAAAAAAA", "###################");

        Graph init_graph;
        init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());
        init_graph.header->add_sequencenames("P");
        init_graph.header->add_sequencenames("Q");
        init_graph.header->add_sequencenames("D");

        init_graph.nodes[0] = std::make_shared<Node>();
        init_graph.nodes[0]->set_id(0);
        init_graph.nodes[0]->set_name("LF");
        init_graph.nodes[0]->set_sequence("AAAAAAAAAAA");

        init_graph.nodes[1] = std::make_shared<Node>();
        init_graph.nodes[1]->set_id(1);
        init_graph.nodes[1]->set_name("P1");
        init_graph.nodes[1]->set_sequence("TTTTTTTT");
        init_graph.nodes[1]->add_sequence_ids(0);

        init_graph.nodes[2] = std::make_shared<Node>();
        init_graph.nodes[2]->set_id(2);
        init_graph.nodes[2]->set_name("Q1");
        init_graph.nodes[2]->set_sequence("GGGGGGGG");
        init_graph.nodes[2]->add_sequence_ids(1);

        init_graph.nodes[3] = std::make_shared<Node>();
        init_graph.nodes[3]->set_id(3);
        init_graph.nodes[3]->set_name("RF");
        init_graph.nodes[3]->set_sequence("AAAAAAAAAAA");

        /**
         *
         * LF --------> RF
         *  |           ^
         *  |           |
         *  *-> P1 -----*
         *  |           |
         *  *-> Q1 -----*
         */

        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[0]->set_from(0);
        init_graph.edges[0]->set_to(1);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[1]->set_from(0);
        init_graph.edges[1]->set_to(2);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[2]->set_from(1);
        init_graph.edges[2]->set_to(3);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[3]->set_from(2);
        init_graph.edges[3]->set_to(3);
        init_graph.edges.emplace_back(new Edge());
        init_graph.edges[4]->set_from(0);
        init_graph.edges[4]->set_to(3);
        init_graph.edges[4]->add_sequence_ids(2);

        auto paths = Json::Value(Json::arrayValue);

        auto p_path = Json::Value(Json::objectValue);
        p_path["path_id"] = "P|1";
        p_path["sequence"] = "P";
        p_path["nodes"] = Json::Value(Json::arrayValue);
        p_path["nodes"].append("LF");
        p_path["nodes"].append("P1");
        p_path["nodes"].append("RF");
        paths.append(p_path);

        auto q_path = Json::Value(Json::objectValue);
        q_path["path_id"] = "Q|1";
        q_path["sequence"] = "Q";
        q_path["nodes"] = Json::Value(Json::arrayValue);
        q_path["nodes"].append("LF");
        q_path["nodes"].append("Q1");
        q_path["nodes"].append("RF");
        paths.append(q_path);

        auto d_path = Json::Value(Json::objectValue);
        d_path["path_id"] = "D|1";
        d_path["sequence"] = "D";
        d_path["nodes"] = Json::Value(Json::arrayValue);
        d_path["nodes"].append("LF");
        d_path["nodes"].append("RF");
        paths.append(d_path);

        LOG()->set_level(spdlog::level::err);

        grm::KmerAligner<16> aligner;
        aligner.setGraph(init_graph, paths);
        for (size_t i = 0; i < reads.size(); ++i)
        {
            aligner.alignRead(reads[i]);
        }
        auto rb_reads = toReadBuffer(reads);
        graph = init_graph;
        paragraph::disambiguateReads(graph, rb_reads, nullptr, nullptr, paths);
        for (size_t i = 0; i < rb_reads.size(); ++i)
        {
            reads[i] = *(rb_reads[i]);
        }
    }
};

TEST_F(KmerAlignerTest, Aligns)
{
    const std::string expected[]
        = { "{\"fragmentId\":\"f1\",\"bases\":\"AAAAAAAATTTTTTTTAAAAAAAA\",\"quals\":\"########################\","
            "\"chromId\":-1,\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":3,"
            "\"graphCigar\":\"0[8M]1[8M]3[8M]\",\"graphMapq\":60,\"graphAlignmentScore\":24,\"isGraphAlignmentUnique\":"
            "true,\"graphNodesSupported\":[\"LF\",\"P1\",\"RF\"],\"graphEdgesSupported\":[\"LF_P1\",\"P1_RF\"],"
            "\"graphSequencesSupported\":[\"P\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f2\",\"bases\":\"AAAAAAATTTTTTTTAAAAAA\",\"quals\":\"#####################\","
            "\"chromId\":-1,\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":4,"
            "\"graphCigar\":\"0[7M]1[8M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":21,\"isGraphAlignmentUnique\":"
            "true,\"isGraphReverseStrand\":true,\"graphNodesSupported\":[\"LF\",\"P1\",\"RF\"],\"graphEdgesSupported\":"
            "[\"LF_P1\",\"P1_RF\"],\"graphSequencesSupported\":[\"P\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f3\",\"bases\":\"AAAAAGGGGGGGGAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":6,\"graphCigar\":\"0[5M]2["
            "8M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":19,\"isGraphAlignmentUnique\":true,"
            "\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],"
            "\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f4\",\"bases\":\"AAAAGGGGGGGGAAAAAA\",\"quals\":\"##################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":7,\"graphCigar\":\"0[4M]2["
            "8M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":18,\"isGraphAlignmentUnique\":true,"
            "\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],"
            "\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f5\",\"bases\":\"AAAAAGGGGGGGGAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":6,\"graphCigar\":\"0[5M]2["
            "8M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":19,\"isGraphAlignmentUnique\":true,"
            "\"isGraphReverseStrand\":true,\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":["
            "\"LF_Q1\",\"Q1_RF\"],\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            // this is the repeat one
            "{\"fragmentId\":\"f6\",\"bases\":\"AAAAAAAAAAAAAAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphCigar\":\"0[11M]3[8M]\","
            "\"graphAlignmentScore\":19,\"graphMappingStatus\":\"BAD_ALIGN\"}" };
    ASSERT_EQ(6ull, reads.size());
    int i = 0;
    for (auto const& read : reads)
    {
        std::string str;
        google::protobuf::util::MessageToJsonString(*((google::protobuf::Message*)&read), &str);
        // std::cerr << str << std::endl;
        ASSERT_EQ(expected[i++], str);
    }
}
