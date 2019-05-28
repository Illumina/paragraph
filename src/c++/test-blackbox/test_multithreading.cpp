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

#include "gtest/gtest.h"

#include "common.hh"
#include <fstream>
#include <sstream>
#include <string>

#include "common/Error.hh"
#include "common/ReadExtraction.hh"
#include "common/Threads.hh"
#include "grm/GraphInput.hh"
#include "paragraph/Disambiguation.hh"

using common::Read;
using common::ReadBuffer;
using paragraph::Parameters;

using paragraph::alignAndDisambiguate;

auto compare_values = [](Json::Value const& lhs, Json::Value const& rhs) {
    for (auto const& name : lhs.getMemberNames())
    {
        if (name.find("FWD") != std::string::npos || name.find("REV") != std::string::npos)
        {
            continue;
        }
        ASSERT_TRUE(rhs.isMember(name));
        ASSERT_EQ(lhs[name].asUInt64(), rhs[name].asUInt64());
    }
};

TEST(Paragraph, AlignsSequentially)
{
    auto logger = LOG();
    // test corner cases for complex graph alignment
    const std::string bam_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.bam";
    const std::string spec_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.json";

    const std::string reference_path = g_testenv->getHG19Path();
    if (reference_path.empty())
    {
        logger->warn("Warning: cannot do round-trip testing for paragraph without hg19 reference file -- please "
                     "specify a location using the HG19 environment variable.");
        return;
    }

    Parameters parameters;
    parameters.set_threads(1);
    common::CPU_THREADS().reset(1);
    parameters.load(spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 0, all_reads);

    auto result = alignAndDisambiguate(parameters, all_reads);

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

TEST(Paragraph, AlignsMultithreaded)
{
    auto logger = LOG();
    // test corner cases for complex graph alignment
    const std::string bam_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.bam";
    const std::string spec_path
        = g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.json";

    const std::string reference_path = g_testenv->getHG19Path();
    if (reference_path.empty())
    {
        logger->warn("Warning: cannot do round-trip testing for paragraph without hg19 reference file -- please "
                     "specify a location using the HG19 environment variable.");
        return;
    }

    Parameters parameters;
    parameters.set_threads(4);
    common::CPU_THREADS().reset(4);
    parameters.load(spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 0, all_reads);

    auto result = alignAndDisambiguate(parameters, all_reads);

    Json::Value expected_result;
    {
        std::ifstream expected_file(
            g_testenv->getBasePath() + "/../share/test-data/paragraph/long-del/chr4-21369091-21376907.paragraph.json");
        ASSERT_TRUE(expected_file.good());
        expected_file >> expected_result;
    }

    ASSERT_TRUE(result.isMember("read_counts_by_node"));
    compare_values(expected_result["read_counts_by_node"], result["read_counts_by_node"]);

    ASSERT_TRUE(result.isMember("read_counts_by_edge"));
    compare_values(expected_result["read_counts_by_edge"], result["read_counts_by_edge"]);

    ASSERT_TRUE(result.isMember("read_counts_by_sequence"));
    for (auto const& expected_name : expected_result["read_counts_by_sequence"].getMemberNames())
    {
        ASSERT_TRUE(result["read_counts_by_sequence"].isMember(expected_name));
        compare_values(
            expected_result["read_counts_by_sequence"][expected_name],
            result["read_counts_by_sequence"][expected_name]);
    }
}
