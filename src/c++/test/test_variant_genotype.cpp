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

#include "genotyping/Genotype.hh"
#include "genotyping/VariantGenotype.hh"
#include "gtest/gtest.h"
#include <string>

using namespace genotyping;

TEST(VariantGenotype, SimplePass)
{
    Genotype gt1;
    gt1.gt = { 1, 1 };
    gt1.gl_name = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    gt1.gl = { -10, -10, -0.1 };

    VariantGenotype vg;
    for (size_t i = 0; i < 2; i++)
    {
        vg.addInfoFromSingleBreakpoint(gt1);
    }
    vg.genotype();

    EXPECT_EQ("1/1", std::string(vg));
}

TEST(VariantGenotype, EdgeWithTwoTags)
{
    Genotype gt1;
    gt1.gt = { 0, 1 };
    gt1.gl_name = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    gt1.gl = { -10, -0.1, -10 };
    Genotype gt2;
    gt2.gt = { 1, 1 };
    gt2.equivalent_gt.push_back({ 0, 0 });
    gt2.equivalent_gt.push_back({ 0, 1 });
    gt2.gl_name = { { 1, 1 }, { 1, 2 }, { 2, 2 } };
    gt2.gl = { 0.1, -10, -10 };

    VariantGenotype vg;
    vg.addInfoFromSingleBreakpoint(gt1);
    vg.addInfoFromSingleBreakpoint(gt2);
    vg.genotype();

    EXPECT_EQ("0/1", std::string(vg));
    EXPECT_EQ("PASS", vg.filter());
}