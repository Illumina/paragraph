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

#include "common/Threads.hh"
#include "grm/Align.hh"
#include "grm/GraphAligner.hh"
#include "grm/GraphInput.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/GraphVariants.hh"

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;

using namespace testing;
using namespace common;
using namespace grm;
using namespace graphtools;

class ParagraphTest : public Test
{
public:
    vector<Read> reads{ 6 };

    Graph graph;

    virtual void SetUp()
    {
        reads[0].setCoreInfo("f1", "AAAAAAAATTTTCTTTAAAAAAAA", "########################");
        reads[1].setCoreInfo("f2", "TTTTTTAAAGAAAATTTTTTT", "#####################");
        reads[2].setCoreInfo("f3", "AAAAAGCGGGGGGAAAAAA", "###################");
        reads[3].setCoreInfo("f4", "AAAAGCGGGGGGAAAAAA", "##################");
        reads[4].setCoreInfo("f5", "TTTTTTCCCCCCGCTTTTT", "###################");
        reads[5].setCoreInfo("f6", "AAAAAAAAAAAAAAAAAAA", "###################");

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

        LOG()->set_level(spdlog::level::err);

        auto rb_reads = toReadBuffer(reads);
        // this ensures results will be in predictable order
        common::CPU_THREADS().reset(1);
        std::list<Path> paths;
        grm::alignReads(&graph, paths, rb_reads, nullptr, false, true, false, false, false);
        paragraph::disambiguateReads(&graph, rb_reads);
        for (size_t i = 0; i < rb_reads.size(); ++i)
        {
            reads[i] = *(rb_reads[i]);
        }
    }
};

TEST_F(ParagraphTest, Aligns)
{
    const std::string expected[]
        = { "{\"fragmentId\":\"f1\",\"bases\":\"AAAAAAAATTTTCTTTAAAAAAAA\",\"quals\":\"########################\","
            "\"chromId\":-1,\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":3,"
            "\"graphCigar\":\"0[8M]1[4M1X3M]3[8M]\",\"graphMapq\":60,\"graphAlignmentScore\":19,"
            "\"isGraphAlignmentUnique\":true,\"graphNodesSupported\":[\"LF\",\"P1\",\"RF\"],\"graphEdgesSupported\":["
            "\"LF_P1\",\"P1_RF\"],\"graphSequencesSupported\":[\"P\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f2\",\"bases\":\"AAAAAAATTTTCTTTAAAAAA\",\"quals\":\"#####################\","
            "\"chromId\":-1,\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":4,"
            "\"graphCigar\":\"0[7M]1[4M1X3M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":16,"
            "\"isGraphAlignmentUnique\":true,\"isGraphReverseStrand\":true,\"graphNodesSupported\":[\"LF\",\"P1\","
            "\"RF\"],\"graphEdgesSupported\":[\"LF_P1\",\"P1_RF\"],\"graphSequencesSupported\":[\"P\"],"
            "\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f3\",\"bases\":\"AAAAAGCGGGGGGAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":6,\"graphCigar\":\"0[5M]2["
            "1M1X6M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":14,\"isGraphAlignmentUnique\":true,"
            "\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],"
            "\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f4\",\"bases\":\"AAAAGCGGGGGGAAAAAA\",\"quals\":\"##################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":7,\"graphCigar\":\"0[4M]2["
            "1M1X6M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":13,\"isGraphAlignmentUnique\":true,"
            "\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":[\"LF_Q1\",\"Q1_RF\"],"
            "\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f5\",\"bases\":\"AAAAAGCGGGGGGAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphPos\":6,\"graphCigar\":\"0[5M]2["
            "1M1X6M]3[6M]\",\"graphMapq\":60,\"graphAlignmentScore\":14,\"isGraphAlignmentUnique\":true,"
            "\"isGraphReverseStrand\":true,\"graphNodesSupported\":[\"LF\",\"Q1\",\"RF\"],\"graphEdgesSupported\":["
            "\"LF_Q1\",\"Q1_RF\"],\"graphSequencesSupported\":[\"Q\"],\"graphMappingStatus\":\"MAPPED\"}",
            "{\"fragmentId\":\"f6\",\"bases\":\"AAAAAAAAAAAAAAAAAAA\",\"quals\":\"###################\",\"chromId\":-1,"
            "\"pos\":-1,\"isFirstMate\":true,\"mateChromId\":-1,\"matePos\":-1,\"graphCigar\":\"0[11M]3[8M]\","
            "\"graphMapq\":60,\"graphAlignmentScore\":19,\"isGraphAlignmentUnique\":true,\"graphNodesSupported\":["
            "\"LF\",\"RF\"],\"graphEdgesSupported\":[\"LF_RF\"],\"graphSequencesSupported\":[\"D\"],"
            "\"graphMappingStatus\":\"MAPPED\"}" };
    ASSERT_EQ(6ull, reads.size());
    int i = 0;

    for (auto const& read : reads)
    {
        Json::Value in_val;
        std::stringstream ss(expected[i++]);
        ss >> in_val;

        const std::string expected_str = in_val.toStyledString();
        const std::string str = read.toJson().toStyledString();
        // std::cerr << str << std::endl;
        ASSERT_EQ(expected_str, str);
    }
}

TEST_F(ParagraphTest, FindsVariants)
{
    paragraph::NodeCandidates variant_candidates;

    for (auto const& read : reads)
    {
        paragraph::updateVariantCandidateLists(&graph, read, variant_candidates);
    }

    ASSERT_EQ(4ull, variant_candidates.size());

    const size_t number_of_variants[] = { 0, 1, 1, 0 };

    const int ref_fwd_expected[4][11] = { {
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              3,
                                              4,
                                              4,
                                              4,
                                              4,
                                          },
                                          {
                                              1,
                                              1,
                                              1,
                                              1,
                                              0,
                                              1,
                                              1,
                                              1,
                                              -1,
                                              -1,
                                              -1,
                                          },
                                          {
                                              2,
                                              0,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              -1,
                                              -1,
                                              -1,
                                          },
                                          {
                                              4,
                                              4,
                                              4,
                                              4,
                                              4,
                                              4,
                                              2,
                                              2,
                                              0,
                                              0,
                                              0,
                                          } };
    const int ref_reverse_expected[4][11] = { {
                                                  0,
                                                  0,
                                                  0,
                                                  0,
                                                  1,
                                                  1,
                                                  2,
                                                  2,
                                                  2,
                                                  2,
                                                  2,
                                              },
                                              {
                                                  1,
                                                  1,
                                                  1,
                                                  1,
                                                  0,
                                                  1,
                                                  1,
                                                  1,
                                                  -1,
                                                  -1,
                                                  -1,
                                              },
                                              {
                                                  1,
                                                  0,
                                                  1,
                                                  1,
                                                  1,
                                                  1,
                                                  1,
                                                  1,
                                                  -1,
                                                  -1,
                                                  -1,
                                              },
                                              {
                                                  2,
                                                  2,
                                                  2,
                                                  2,
                                                  2,
                                                  2,
                                                  0,
                                                  0,
                                                  0,
                                                  0,
                                                  0,
                                              } };

    // clang-format off
    const int nonref_fwd_expected[4][11] = {
        {
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        },
        {
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
        },
        {
            0,
            2,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        },
        {
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        },
    };
    const int nonref_reverse_expected[4][11] = { {
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                 },
                                                 {
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     1,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                 },
                                                 {
                                                     0,
                                                     1,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                 },
                                                 {
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                 } };
    // clang-format on

    for (uint64_t i = 0; i <= 3; ++i)
    {
        auto cand_id = variant_candidates.find(i);
        ASSERT_NE(cand_id, variant_candidates.end());
        const auto& candidates_for_node = cand_id->second;
        const auto& variant_list = candidates_for_node.getVariants();
        ASSERT_EQ(number_of_variants[i], variant_list.size());

        for (size_t pos = 0; pos < candidates_for_node.getReference().size(); ++pos)
        {
            ASSERT_EQ(ref_fwd_expected[i][pos], candidates_for_node.getRefPileup((int)pos).stranded_DP[0]);
            ASSERT_EQ(ref_reverse_expected[i][pos], candidates_for_node.getRefPileup((int)pos).stranded_DP[1]);
            ASSERT_EQ(nonref_fwd_expected[i][pos], candidates_for_node.getNonrefPileup((int)pos).stranded_DP[0]);
            ASSERT_EQ(nonref_reverse_expected[i][pos], candidates_for_node.getNonrefPileup((int)pos).stranded_DP[1]);
        }

        // std::cerr << "Node " << i << std::endl;
        // std::cerr << "ref fwd: ";
        // for(size_t pos = 0; pos < candidates_for_node.getReference().size(); ++pos)
        // {
        //     std::cerr << candidates_for_node.getRefPileup((int) pos).stranded_DP[0] << ", ";
        // }
        // std::cerr << std::endl;
        // std::cerr << "ref reverse: ";
        // for(size_t pos = 0; pos < candidates_for_node.getReference().size(); ++pos)
        // {
        //     std::cerr << candidates_for_node.getRefPileup(pos).stranded_DP[1] << ", ";
        // }
        // std::cerr << std::endl;
    }
}
