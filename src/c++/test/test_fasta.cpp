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
 *  \brief Test Fasta file reading
 *
 * \file test_fasta.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Error.hh"
#include "common/Fasta.hh"
#include "gtest/gtest.h"

#include "common.hh"

#include <iostream>

using namespace common;

TEST(Fasta, ReadsFasta)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";
    FastaFile f(tp.c_str());
    ASSERT_EQ("CCAAA", f.query("chrQ:5-9"));
}

TEST(Fasta, DoesntReadFastaPastEnd)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";
    FastaFile f(tp.c_str());
    ASSERT_EQ("", f.query("chrS:151"));
}

TEST(Fasta, ReadsLineWrapped)
{
    // test if we retrieve correct results if the line wraps in the
    // input fasta
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";
    FastaFile f(tp.c_str());
    ASSERT_EQ("TTCAGTGTTCTTTTTACTTAAGCCTTCTTTCTGGTACGTATGAGGTGTGCTGTCATACGTATGTCGTTATT", f.query("chrT:50-120"));
}

TEST(Fasta, ReadsLineWrappedEOC)
{
    // test if we retrieve correct results if the line wraps in the
    // input fasta and we encounter the contig end
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";
    FastaFile f(tp.c_str());
    ASSERT_EQ(
        "TTCAGTGTTCTTTTTACTTAAGCCTTCTTTCTGGTACGTATGAGGTGTGCTGTCATACGTATGTCGTTATTTCTCTTTTCAGATTAGTCATGTCCCTAATT",
        f.query("chrT:50-200"));
}
