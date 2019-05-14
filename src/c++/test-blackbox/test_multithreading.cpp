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
