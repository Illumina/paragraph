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

#include "common/ReadExtraction.hh"
#include "common/Threads.hh"
#include "paragraph/Disambiguation.hh"
#include "paragraph/Parameters.hh"

// Error.hh always needs to go last
#include "common/Error.hh"

using namespace paragraph;

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

    Parameters parameters(
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, true);

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    ASSERT_EQ(output["read_counts_by_edge"]["total"].asUInt64(), 30ull);
    ASSERT_EQ(output["read_counts_by_edge"]["total:FWD"].asUInt64(), 15ull);
    ASSERT_EQ(output["read_counts_by_edge"]["total:REV"].asUInt64(), 15ull);

    ASSERT_EQ(output["read_counts_by_node"]["total"].asUInt64(), 30ull);
    ASSERT_EQ(output["read_counts_by_node"]["total:FWD"].asUInt64(), 15ull);
    ASSERT_EQ(output["read_counts_by_node"]["total:REV"].asUInt64(), 15ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 15ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:FWD"].asUInt64(), 8ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:REV"].asUInt64(), 7ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total"].asUInt64(), 8ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total:FWD"].asUInt64(), 4ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total:REV"].asUInt64(), 4ull);
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

    Parameters parameters(
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, true);

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    ASSERT_EQ(output["read_counts_by_edge"]["total"].asUInt64(), 203ull);
    ASSERT_EQ(output["read_counts_by_edge"]["total:READS"].asUInt64(), 212ull);
    ASSERT_EQ(output["read_counts_by_edge"]["total:FWD"].asUInt64(), 105ull);
    ASSERT_EQ(output["read_counts_by_edge"]["total:REV"].asUInt64(), 107ull);

    ASSERT_EQ(output["read_counts_by_node"]["total"].asUInt64(), 203ull);
    ASSERT_EQ(output["read_counts_by_node"]["total:READS"].asUInt64(), 212ull);
    ASSERT_EQ(output["read_counts_by_node"]["total:FWD"].asUInt64(), 105ull);
    ASSERT_EQ(output["read_counts_by_node"]["total:REV"].asUInt64(), 107ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 133ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:READS"].asUInt64(), 141ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:FWD"].asUInt64(), 71ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:REV"].asUInt64(), 70ull);

    ASSERT_FALSE(output["read_counts_by_sequence"].isMember("ALT"));
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

    Parameters parameters(
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, true);

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 12ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:FWD"].asUInt64(), 6ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:REV"].asUInt64(), 6ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total"].asUInt64(), 0ull);
}
