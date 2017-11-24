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

#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "genotyping/EMbreakpointGenotyper.hh"
#include "genotyping/Genotype.hh"
#include "genotyping/Idxdepth.hh"

using namespace genotyping;

TEST(emGenotyper, FixParameter)
{
    EMbreakpointGenotyper genotyper(0.01, 2, 2, 10, 16);
    std::vector<EdgeCounts> counts;
    counts.push_back({ 20, 20 });
    counts.push_back({ 1, 40 });
    Idxdepth idxdepth(100, 40, 2);
    std::vector<Genotype> result = genotyper.genotype(counts, idxdepth, false);
    EXPECT_EQ("0/1", (std::string)result[0]);
    EXPECT_EQ("1/1", (std::string)result[1]);
}

TEST(emGenotyper, emParameter)
{
    EMbreakpointGenotyper genotyper(0.01, 2, 2, 10, 16);
    std::vector<EdgeCounts> counts;
    counts.push_back({ 20, 20 });
    counts.push_back({ 1, 40 });
    Idxdepth idxdepth(100, 40, 2);
    std::vector<Genotype> result = genotyper.genotype(counts, idxdepth, true);
    EXPECT_EQ("0/1", (std::string)result[0]);
    EXPECT_EQ("1/1", (std::string)result[1]);
}

TEST(emGenotyper, MissingIfEmpty)
{
    EMbreakpointGenotyper genotyper(0.01, 2, 2, 10, 16);
    std::vector<EdgeCounts> counts;
    counts.push_back({});
    counts.push_back({ 1, 40 });
    Idxdepth idxdepth(100, 40, 2);
    std::vector<Genotype> result = genotyper.genotype(counts, idxdepth, true);
    EXPECT_EQ("./.", (std::string)result[0]);
    EXPECT_EQ("1/1", (std::string)result[1]);
}
