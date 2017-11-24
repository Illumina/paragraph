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
