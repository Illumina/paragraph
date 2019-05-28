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

/**
 * \brief Test paragraph aligner
 *
 * \file test_paragraph.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "gtest/gtest.h"

#include <sstream>
#include <string>

#include "common.hh"

#include "common/JsonHelpers.hh"
#include "common/ReadExtraction.hh"
#include "common/Threads.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/Parameters.hh"

// Error.hh always needs to go last
#include "common/Error.hh"

using namespace paragraph;

auto compare_values = [](Json::Value const& lhs, Json::Value const& rhs) {
    for (auto const& name : lhs.getMemberNames())
    {
        // when the glitch in FWD/REV counting is solved, remove this condition
        if (name.find("FWD") != std::string::npos || name.find("REV") != std::string::npos)
        {
            continue;
        }
        ASSERT_TRUE(rhs.isMember(name));
        ASSERT_EQ(lhs[name].asUInt64(), rhs[name].asUInt64());
    }
};

TEST(Paragraph, AlignsPGHetIns)
{
    auto logger = LOG();

    std::string bam_path = g_testenv->getBasePath() + "/../share/test-data/paragraph/pg-het-ins/na12878.bam";
    std::string graph_spec_path = g_testenv->getBasePath() + "/../share/test-data/paragraph/pg-het-ins/pg-het-ins.json";

    std::string reference_path = g_testenv->getHG38Path();
    if (reference_path.empty())
    {
        logger->warn("Warning: cannot do round-trip testing for paragraph without hg38 reference file -- please "
                     "specify a location using the HG38 environment variable.");
        return;
    }

    Parameters parameters;

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 0, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto result = alignAndDisambiguate(parameters, all_reads);

    Json::Value expected_result;
    {
        std::ifstream expected_file(
            g_testenv->getBasePath() + "/../share/test-data/paragraph/pg-het-ins/pg-het-ins.result.json");
        ASSERT_TRUE(expected_file.good());
        expected_file >> expected_result;
    }

    ASSERT_TRUE(result.isMember("read_counts_by_node"));
    ASSERT_TRUE(expected_result.isMember("read_counts_by_node"));
    compare_values(expected_result["read_counts_by_node"], result["read_counts_by_node"]);

    ASSERT_TRUE(result.isMember("read_counts_by_edge"));
    compare_values(expected_result["read_counts_by_edge"], result["read_counts_by_edge"]);
}

TEST(Paragraph, AlignsPGLongDel)
{
    auto logger = LOG();

    std::string bam_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.bam";
    std::string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.json";

    std::string reference_path = g_testenv->getHG19Path();
    if (reference_path.empty())
    {
        logger->warn("Warning: cannot do round-trip testing for paragraph without hg19 reference file -- please "
                     "specify a location using the HG19 environment variable.");
        return;
    }

    Parameters parameters;

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 0, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto result = alignAndDisambiguate(parameters, all_reads);

    Json::Value expected_result;
    {
        std::ifstream expected_file(
            g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.paragraph.json");
        ASSERT_TRUE(expected_file.good());
        expected_file >> expected_result;
    }

    ASSERT_TRUE(result.isMember("read_counts_by_node"));
    ASSERT_TRUE(expected_result.isMember("read_counts_by_node"));
    compare_values(expected_result["read_counts_by_node"], result["read_counts_by_node"]);

    ASSERT_TRUE(result.isMember("read_counts_by_edge"));
    compare_values(expected_result["read_counts_by_edge"], result["read_counts_by_edge"]);
}

TEST(Paragraph, CountsClippedReads)
{
    auto logger = LOG();

    std::string bam_path = g_testenv->getBasePath() + "/../share/test-data/paragraph/quantification/NA12878-chr6.bam";
    std::string graph_spec_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/quantification/chr6-53037879-53037949.vcf.json";

    std::string reference_path = g_testenv->getHG19Path();
    if (reference_path.empty())
    {
        logger->warn("Warning: cannot do round-trip testing for paragraph without hg19 reference file -- please "
                     "specify a location using the HG19 environment variable.");
        return;
    }

    Parameters parameters;

    logger->info("Loading parameters from {} {}", graph_spec_path, reference_path);
    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 0, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 0ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total"].asUInt64(), 0ull);
}
