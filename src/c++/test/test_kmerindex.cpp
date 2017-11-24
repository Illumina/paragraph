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
#include "graphs/KmerIndex.hh"

#include <grm/Align.hh>
#include <iostream>
#include <map>
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
using namespace graphs;

TEST(KmerIndex, NucleotideCount)
{
    Graph init_graph;
    init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());
    init_graph.header->add_sequencenames("P");
    init_graph.header->add_sequencenames("Q");
    init_graph.header->add_sequencenames("D");

    init_graph.nodes[0] = std::make_shared<Node>();
    init_graph.nodes[0]->set_id(0);
    init_graph.nodes[0]->set_name("LF");
    init_graph.nodes[0]->set_sequence("AAAA");

    init_graph.nodes[1] = std::make_shared<Node>();
    init_graph.nodes[1]->set_id(1);
    init_graph.nodes[1]->set_name("P1");
    init_graph.nodes[1]->set_sequence("TTTT");
    init_graph.nodes[1]->add_sequence_ids(0);

    init_graph.nodes[2] = std::make_shared<Node>();
    init_graph.nodes[2]->set_id(2);
    init_graph.nodes[2]->set_name("Q1");
    init_graph.nodes[2]->set_sequence("GGGG");
    init_graph.nodes[2]->add_sequence_ids(1);

    init_graph.nodes[3] = std::make_shared<Node>();
    init_graph.nodes[3]->set_id(3);
    init_graph.nodes[3]->set_name("RF");
    init_graph.nodes[3]->set_sequence("AAAAA");

    /**
     *
     * LF           RF
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

    // the graph above should only have unique 5-mers
    WalkableGraph graph(init_graph);
    KmerIndex index(graph, 1);
    const std::set<std::string> expected_kmers{ "A", "G", "T" };
    for (auto const& k : index.kmers())
    {
        ASSERT_EQ(1ull, expected_kmers.count(k));
        if (k == "A")
        {
            ASSERT_EQ(4 + 5, index.kmerCount(k));
        }
        else if (k == "G")
        {
            ASSERT_EQ(4, index.kmerCount(k));
        }
        else if (k == "T")
        {
            ASSERT_EQ(4, index.kmerCount(k));
        }
    }
    ASSERT_EQ(expected_kmers.size(), index.kmers().size());
    ASSERT_EQ(0ull, index.badkmers().size());
}

TEST(KmerIndex, AllUniqueKmers)
{
    Graph init_graph;
    init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());
    init_graph.header->add_sequencenames("P");
    init_graph.header->add_sequencenames("Q");
    init_graph.header->add_sequencenames("D");

    init_graph.nodes[0] = std::make_shared<Node>();
    init_graph.nodes[0]->set_id(0);
    init_graph.nodes[0]->set_name("LF");
    init_graph.nodes[0]->set_sequence("AAAA");

    init_graph.nodes[1] = std::make_shared<Node>();
    init_graph.nodes[1]->set_id(1);
    init_graph.nodes[1]->set_name("P1");
    init_graph.nodes[1]->set_sequence("TTTT");
    init_graph.nodes[1]->add_sequence_ids(0);

    init_graph.nodes[2] = std::make_shared<Node>();
    init_graph.nodes[2]->set_id(2);
    init_graph.nodes[2]->set_name("Q1");
    init_graph.nodes[2]->set_sequence("GGGG");
    init_graph.nodes[2]->add_sequence_ids(1);

    init_graph.nodes[3] = std::make_shared<Node>();
    init_graph.nodes[3]->set_id(3);
    init_graph.nodes[3]->set_name("RF");
    init_graph.nodes[3]->set_sequence("AAAAA");

    /**
     *
     * LF           RF
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

    // the graph above should only have unique 5-mers
    WalkableGraph graph(init_graph);
    KmerIndex index(graph, 5);
    const std::set<std::string> expected_kmers{ "AAAAT", "AAAAG", "AAATT", "AAAGG", "AATTT", "AAGGG",
                                                "ATTTT", "AGGGG", "TTTTA", "GGGGA", "TTTAA", "GGGAA",
                                                "TTAAA", "GGAAA", "TAAAA", "GAAAA", "AAAAA" };
    for (auto const& k : index.kmers())
    {
        ASSERT_EQ(1ull, expected_kmers.count(k));
        ASSERT_EQ(1, index.kmerCount(k));
    }
    ASSERT_EQ(expected_kmers.size(), index.kmers().size());
    ASSERT_EQ(0ull, index.badkmers().size());
}

TEST(KmerIndex, BadKmers)
{
    Graph init_graph;
    init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());
    init_graph.header->add_sequencenames("P");
    init_graph.header->add_sequencenames("Q");
    init_graph.header->add_sequencenames("D");

    init_graph.nodes[0] = std::make_shared<Node>();
    init_graph.nodes[0]->set_id(0);
    init_graph.nodes[0]->set_name("LF");
    init_graph.nodes[0]->set_sequence("AAAA");

    init_graph.nodes[1] = std::make_shared<Node>();
    init_graph.nodes[1]->set_id(1);
    init_graph.nodes[1]->set_name("P1");
    init_graph.nodes[1]->set_sequence("GTTTT");
    init_graph.nodes[1]->add_sequence_ids(0);

    init_graph.nodes[2] = std::make_shared<Node>();
    init_graph.nodes[2]->set_id(2);
    init_graph.nodes[2]->set_name("Q1");
    init_graph.nodes[2]->set_sequence("GGGGT");
    init_graph.nodes[2]->add_sequence_ids(1);

    init_graph.nodes[3] = std::make_shared<Node>();
    init_graph.nodes[3]->set_id(3);
    init_graph.nodes[3]->set_name("RF");
    init_graph.nodes[3]->set_sequence("AAAAA");

    /**
     *
     * LF           RF
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

    // the graph above should only have unique 5-mers
    WalkableGraph graph(init_graph);
    KmerIndex index(graph, 5);
    const std::map<std::string, int> expected_kmers{
        { "AAAGT", 1 }, { "AAAGG", 1 }, { "AAGTT", 1 }, { "AAGGG", 1 }, { "AGTTT", 1 },
        { "GTTTT", 1 }, { "AGGGG", 1 }, { "TTTTA", 1 }, { "GGGGT", 1 }, { "TTTAA", 1 },
        { "GGGTA", 1 }, { "TTAAA", 1 }, { "GGTAA", 1 }, { "TAAAA", 2 }, // two paths / locations for this one
        { "GTAAA", 1 }, { "AAAAA", 1 },
    };
    for (auto const& k : index.kmers())
    {
        // std::cerr << k << std::endl;
        ASSERT_EQ(1ull, expected_kmers.count(k));
        ASSERT_EQ(expected_kmers.find(k)->second, index.kmerCount(k));
    }
    ASSERT_EQ(expected_kmers.size(), index.kmers().size());
    const std::map<std::string, int> expected_badkmers{
        { "AAAAG", 1 } // matches both subsequent nodes, but only at one location
    };
    ASSERT_EQ(expected_badkmers.size(), index.badkmers().size());
    for (auto const& k : index.badkmers())
    {
        // std::cerr << k << std::endl;
        ASSERT_EQ(1ull, expected_badkmers.count(k));
        ASSERT_EQ(expected_badkmers.find(k)->second, index.kmerCount(k));
    }
}

TEST(KmerIndex, KmerLocations)
{
    Graph init_graph;
    init_graph.header = std::unique_ptr<GraphHeader>(new GraphHeader());
    init_graph.header->add_sequencenames("P");
    init_graph.header->add_sequencenames("Q");
    init_graph.header->add_sequencenames("D");

    init_graph.nodes[0] = std::make_shared<Node>();
    init_graph.nodes[0]->set_id(0);
    init_graph.nodes[0]->set_name("LF");
    init_graph.nodes[0]->set_sequence("CAAA");

    init_graph.nodes[1] = std::make_shared<Node>();
    init_graph.nodes[1]->set_id(1);
    init_graph.nodes[1]->set_name("P1");
    init_graph.nodes[1]->set_sequence("GTTTT");
    init_graph.nodes[1]->add_sequence_ids(0);

    init_graph.nodes[2] = std::make_shared<Node>();
    init_graph.nodes[2]->set_id(2);
    init_graph.nodes[2]->set_name("Q1");
    init_graph.nodes[2]->set_sequence("GGGGT");
    init_graph.nodes[2]->add_sequence_ids(1);

    init_graph.nodes[3] = std::make_shared<Node>();
    init_graph.nodes[3]->set_id(3);
    init_graph.nodes[3]->set_name("RF");
    init_graph.nodes[3]->set_sequence("AAAAAAAAAA");

    /**
     *
     * LF --------- RF
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

    // the graph above should only have unique 5-mers
    WalkableGraph graph(init_graph);
    KmerIndex index(graph, 9);
    const std::map<std::string, std::pair<uint64_t, int>> expected_kmers{
        { "GGTAAAAAA", { 2, 2 } }, { "GGGTAAAAA", { 2, 1 } }, { "TTAAAAAAA", { 1, 3 } }, { "TTTAAAAAA", { 1, 2 } },
        { "TTTTAAAAA", { 1, 1 } }, { "GTTTTAAAA", { 1, 0 } }, { "AAGTTTTAA", { 0, 2 } }, { "AAAGGGGTA", { 0, 1 } },
        { "AGTTTTAAA", { 0, 3 } }, { "TAAAAAAAA", { 1, 4 } }, { "AAGGGGTAA", { 0, 2 } }, { "AAAGTTTTA", { 0, 1 } },
        { "GGGGTAAAA", { 2, 0 } }, { "CAAAGTTTT", { 0, 0 } }, { "GTAAAAAAA", { 2, 3 } }, { "AGGGGTAAA", { 0, 3 } },
        { "CAAAGGGGT", { 0, 0 } }, { "CAAAAAAAA", { 0, 0 } }, { "AAAAAAAAA", { 0, 1 } },
    };
    for (auto const& k : index.kmers())
    {
        // std::cerr << k << std::endl;
        ASSERT_EQ(1ull, expected_kmers.count(k));
        auto expected_pos = expected_kmers.find(k)->second;
        index.search(k);
        if (k != "AAAAAAAAA")
        {
            // TAAAAAAAA occurs twice, once from P1 and once from Q1
            if (k != "TAAAAAAAA")
            {
                ASSERT_EQ(1ull, index.count());
                auto match = std::move(index.matches().front());
                ASSERT_EQ(match->node(), expected_pos.first);
                ASSERT_EQ(match->pos(), expected_pos.second);
                ASSERT_EQ(match->length(), 9ull);
            }
            else
            {
                // this kmer may start at two different nodes
                const std::set<std::pair<uint64_t, int>> expected_two_positions{ { 1, 4 }, { 2, 4 } };
                ASSERT_EQ(expected_two_positions.size(), index.count());
                for (const auto& match : index.matches())
                {
                    ASSERT_EQ(1ull, expected_two_positions.count(std::make_pair(match->node(), match->pos())));
                }
            }
        }
        else
        {
            const std::set<std::pair<uint64_t, int>> expected_five_positions{
                { 0, 1 }, { 0, 2 }, { 0, 3 }, { 3, 0 }, { 3, 1 }
            };
            ASSERT_EQ(expected_five_positions.size(), index.count());
            for (const auto& match : index.matches())
            {
                ASSERT_EQ(1ull, expected_five_positions.count(std::make_pair(match->node(), match->pos())));
            }
        }
    }
    ASSERT_EQ(expected_kmers.size(), index.kmers().size());
    ASSERT_EQ(0ull, index.badkmers().size());
}
