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

#include "common/JsonHelpers.hh"
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
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, false, true);

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 1000, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    auto& edges = output["read_counts_by_edge"];
    ASSERT_EQ(edges["chr1:939571-939570:CCCTGGAGGACC_ref-chr1:939571-939720"].asUInt64(), 10ull);
    ASSERT_EQ(edges["chr1:939571-939570:CCCTGGAGGACC_ref-chr1:939571-939720:FWD"].asUInt64(), 4ull);
    ASSERT_EQ(edges["chr1:939571-939570:CCCTGGAGGACC_ref-chr1:939571-939720:REV"].asUInt64(), 6ull);
    ASSERT_EQ(edges["ref-chr1:939420-939570_chr1:939571-939570:CCCTGGAGGACC"].asUInt64(), 14ull);
    ASSERT_EQ(edges["ref-chr1:939420-939570_ref-chr1:939571-939720"].asUInt64(), 12ull);

    auto& nodes = output["read_counts_by_node"];
    ASSERT_EQ(nodes["chr1:939571-939570:CCCTGGAGGACC"].asUInt64(), 15ull);
    ASSERT_EQ(nodes["chr1:939571-939570:CCCTGGAGGACC:FWD"].asUInt64(), 7ull);
    ASSERT_EQ(nodes["chr1:939571-939570:CCCTGGAGGACC:REV"].asUInt64(), 8ull);
    ASSERT_EQ(nodes["ref-chr1:939420-939570"].asUInt64(), 28ull);
    ASSERT_EQ(nodes["ref-chr1:939571-939720"].asUInt64(), 29ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 12ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:FWD"].asUInt64(), 7ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total:REV"].asUInt64(), 5ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["ref-chr1:939420-939570"].asUInt64(), 12ull);
    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total"].asUInt64(), 16ull);
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
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, false, true);

    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 1000, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    auto& edges = output["read_counts_by_edge"];
    ASSERT_EQ(edges["ref-chr4:21368941-21369090_ref-chr4:21369091-21369240"].asUInt64(), 44ull);
    ASSERT_EQ(edges["ref-chr4:21376758-21376907_ref-chr4:21376908-21377057"].asUInt64(), 42ull);

    auto& nodes = output["read_counts_by_node"];
    ASSERT_EQ(nodes["ref-chr4:21368941-21369090"].asUInt64(), 70ull);
    ASSERT_EQ(nodes["ref-chr4:21369091-21369240"].asUInt64(), 68ull);
    ASSERT_EQ(nodes["ref-chr4:21376758-21376907"].asUInt64(), 93ull);
    ASSERT_EQ(nodes["ref-chr4:21376908-21377057"].asUInt64(), 82ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 86ull);
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
        10000, 3, 0.01f, 0.8f, Parameters::output_options::ALL & ~Parameters::output_options::HAPLOTYPES, false, true,
        true, true);

    logger->info("Loading parameters from {} {}", graph_spec_path, reference_path);
    parameters.load(graph_spec_path, reference_path);

    common::ReadBuffer all_reads;
    common::extractReads(
        bam_path, "", reference_path, parameters.target_regions(), (int)parameters.max_reads(), 1000, all_reads);
    // this ensures results will be in predictable order
    common::CPU_THREADS().reset(1);
    const auto output = alignAndDisambiguate(parameters, all_reads);

    ASSERT_EQ(output["read_counts_by_sequence"]["REF"]["total"].asUInt64(), 0ull);

    ASSERT_EQ(output["read_counts_by_sequence"]["ALT"]["total"].asUInt64(), 0ull);
}
