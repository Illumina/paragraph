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

#include "common/JsonHelpers.hh"
#include "common/ReadReader.hh"
#include "grm/KlibAligner.hh"
#include "paragraph/Disambiguation.hh"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "common/Error.hh"
#include "gtest/gtest.h"

using graphtools::Graph;

using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;

using namespace testing;
using namespace common;

class KlibAlignerTest : public Test
{
public:
    vector<Read> reads{ 7 };

    Graph graph{ 0 };

    void SetUp() override
    {
        reads[0].setCoreInfo("f1", "AAAAAAAATTTTTTTTAAAAAAAA", "########################");
        reads[1].setCoreInfo("f2", "TTTTTTAAAAAAAATTTTTTT", "#####################");
        reads[2].setCoreInfo("f3", "AAAAAGGGGGGGGAAAAAA", "###################");
        reads[3].setCoreInfo("f4", "AAAAGGGGGGGGAAAAAA", "##################");
        reads[4].setCoreInfo("f5", "TTTTTTCCCCCCCCTTTTT", "###################");
        // test for full clipping of flank nodes
        reads[5].setCoreInfo("f7", "TTTTTTCCCCCCCCGGGGG", "###################");
        reads[6].setCoreInfo("f8", "GGGGGGCCCCCCCCTTTTT", "###################");
        // this read won't align uniquely here since we can move the A matches around
        // reads[7].setCoreInfo("f6", "AAAAAAAAAAAAAAAAAAA", "###################");
        for (auto& read : reads)
        {
            read.set_mate_chrom_id(0);
            read.set_mate_pos(0);
        }

        graph = Graph{ 4 };
        graph.setNodeName(0, "LF");
        graph.setNodeSeq(0, "AAAAAAAAAAA");

        graph.setNodeName(1, "P1");
        graph.setNodeSeq(1, "TTTTTTTT");

        graph.setNodeName(2, "Q1");
        graph.setNodeSeq(2, "GGGGGGGG");

        graph.setNodeName(3, "RF");
        graph.setNodeSeq(3, "AAAAAAAAAAA");

        /**
         *
         * LF --------> RF
         *  |           ^
         *  |           |
         *  *-> P1 -----*
         *  |           |
         *  *-> Q1 -----*
         */

        graph.addEdge(0, 1);
        graph.addEdge(0, 2);
        graph.addEdge(0, 3);
        graph.addEdge(1, 3);
        graph.addEdge(2, 3);
        graph.addLabelToEdge(0, 1, "P");
        graph.addLabelToEdge(1, 3, "P");
        graph.addLabelToEdge(0, 2, "Q");
        graph.addLabelToEdge(2, 3, "Q");
        graph.addLabelToEdge(0, 3, "D");

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

        grm::KlibAligner aligner;
        const std::list<graphtools::Path> grmPaths = grm::pathsFromJson(&graph, paths);
        aligner.setGraph(&graph, grmPaths);
        for (auto& read : reads)
        {
            aligner.alignRead(read);
        }
        auto rb_reads = toReadBuffer(reads);
        paragraph::disambiguateReads(&graph, rb_reads);
        for (size_t i = 0; i < rb_reads.size(); ++i)
        {
            reads[i] = *(rb_reads[i]);
        }
    }
};

TEST_F(KlibAlignerTest, Aligns)
{
    const std::string expected[]
        // clang-format off
        = { "{\"bases\":\"AAAAAAAATTTTTTTTAAAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f1\","
                    "\"graphAlignmentScore\":24,\"graphCigar\":\"0[8M]1[8M]3[8M]\","
                    "\"graphEdgesSupported\":[\"LF_P1\",\"P1_RF\"],\"graphMappingStatus\":\"MAPPED\","
                    "\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"P1\",\"RF\"],"
                    "\"graphPos\":3,\"graphSequencesSupported\":[\"P\"],"
                    "\"isFirstMate\":true,"
                    "\"isGraphAlignmentUnique\":true,"
                    "\"pos\":-1,\"quals\":\"########################\"}",
            "{\"bases\":\"AAAAAAATTTTTTTTAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f2\","
                    "\"graphAlignmentScore\":21,\"graphCigar\":\"0[7M]1[8M]3[6M]\","
                    "\"graphEdgesSupported\":[\"LF_P1\",\"P1_RF\"],\"graphMappingStatus\":\"MAPPED\","
                    "\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"P1\",\"RF\"],"
                    "\"graphPos\":4,\"graphSequencesSupported\":[\"P\"],"
                    "\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,"
                    "\"isGraphReverseStrand\":true,\"pos\":-1,\"quals\":\"#####################\"}",
            "{\"bases\":\"AAAAAGGGGGGGGAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f3\","
                    "\"graphAlignmentScore\":19,\"graphCigar\":\"0[5M]2[8M]3[6M]\","
                    "\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],\"graphMappingStatus\":\"MAPPED\","
                    "\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],"
                    "\"graphPos\":6,\"graphSequencesSupported\":[\"Q\"],"
                    "\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,"
                    "\"pos\":-1,\"quals\":\"###################\"}",
            "{\"bases\":\"AAAAGGGGGGGGAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f4\",\"graphAlignmentScore\":18,\"graphCigar\":\"0[4M]2[8M]3[6M]\",\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],\"graphMappingStatus\":\"MAPPED\",\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphPos\":7,\"graphSequencesSupported\":[\"Q\"],\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,\"pos\":-1,\"quals\":\"##################\"}",
            "{\"bases\":\"AAAAAGGGGGGGGAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f5\",\"graphAlignmentScore\":19,\"graphCigar\":\"0[5M]2[8M]3[6M]\",\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],\"graphMappingStatus\":\"MAPPED\",\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphPos\":6,\"graphSequencesSupported\":[\"Q\"],\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,\"isGraphReverseStrand\":true,\"pos\":-1,\"quals\":\"###################\"}",
            // test for full clipping of flank nodes
            "{\"bases\":\"CCCCCGGGGGGGGAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f7\",\"graphAlignmentScore\":14,\"graphCigar\":\"2[5S8M]3[6M]\",\"graphEdgesSupported\":[\"Q1_RF\"],\"graphMappingStatus\":\"MAPPED\",\"graphMapq\":60,\"graphNodesSupported\":[\"Q1\",\"RF\"],\"graphSequencesSupported\":[\"Q\"],\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,\"isGraphReverseStrand\":true,\"pos\":-1,\"quals\":\"###################\"}",
            "{\"bases\":\"AAAAAGGGGGGGGCCCCCC\",\"chromId\":-1,\"fragmentId\":\"f8\",\"graphAlignmentScore\":13,\"graphCigar\":\"0[5M]2[8M6S]\",\"graphEdgesSupported\":[\"LF_Q1\"],\"graphMappingStatus\":\"MAPPED\",\"graphMapq\":60,\"graphNodesSupported\":[\"LF\",\"Q1\"],\"graphPos\":6,\"graphSequencesSupported\":[\"Q\"],\"isFirstMate\":true,\"isGraphAlignmentUnique\":true,\"isGraphReverseStrand\":true,\"pos\":-1,\"quals\":\"###################\"}"
            // this is the repeat one
//            "{\"bases\":\"AAAAAAAAAAAAAAAAAAA\",\"chromId\":-1,\"fragmentId\":\"f6\",\"graphAlignmentScore\":19,\"graphCigar\":\"0[11M]3[8M]\",\"graphMappingStatus\":\"BAD_ALIGN\",\"isFirstMate\":true,\"pos\":-1,\"quals\":\"###################\"}"//,
            };
    // clang-format on
    ASSERT_EQ(7ull, reads.size());
    int i = 0;

    for (auto const& read : reads)
    {
        const std::string str = common::writeJson(read.toJson(), false);
        // std::cerr << str << std::endl;
        ASSERT_EQ(expected[i++], str);
    }
}
