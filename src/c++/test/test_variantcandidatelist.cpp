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
 *  \brief Test cases for haplotype comparison
 *
 *
 * \file test_refvar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common.hh"
#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "variant/Variant.hh"

using namespace variant;

TEST(VariantCandidateList, BasicCandidateListTest)
{
    VariantCandidateList vl("CCACATATATATATATATATA");

    // add a SNP
    RefVar rv;
    rv.start = 3;
    rv.end = 3;
    rv.alt = "T";
    vl.addRefVarObservation(rv, false);
    vl.addRefVarObservation(rv, true);
    vl.addRefVarObservation(rv, true, 0, 10);

    rv.start = 3;
    rv.end = 5;
    rv.alt = ".";

    vl.addRefVarObservation(rv, false);
    vl.addRefVarObservation(rv, false);
    vl.addRefVarObservation(rv, false);
    vl.addRefVarObservation(rv, false, 0, 10);
    vl.addRefVarObservation(rv, true);
    vl.addRefVarObservation(rv, true);
    vl.addRefVarObservation(rv, true);
    vl.addRefVarObservation(rv, true);
    vl.addRefVarObservation(rv, true, 0, 10);

    ASSERT_EQ(0, vl.getRefPileup(2).stranded_DP[0]);
    ASSERT_EQ(0, vl.getRefPileup(2).stranded_DP[1]);
    ASSERT_EQ(4, vl.getRefPileup(3).stranded_DP[0]);
    ASSERT_EQ(5, vl.getRefPileup(3).stranded_DP[1]);
    ASSERT_EQ(4, vl.getRefPileup(4).stranded_DP[0]);
    ASSERT_EQ(5, vl.getRefPileup(4).stranded_DP[1]);
    ASSERT_EQ(4, vl.getRefPileup(5).stranded_DP[0]);
    ASSERT_EQ(5, vl.getRefPileup(5).stranded_DP[1]);
    ASSERT_EQ(0, vl.getRefPileup(6).stranded_DP[0]);
    ASSERT_EQ(0, vl.getRefPileup(6).stranded_DP[1]);
    ASSERT_EQ(0, vl.getNonrefPileup(2).stranded_DP[0]);
    ASSERT_EQ(0, vl.getNonrefPileup(2).stranded_DP[1]);
    ASSERT_EQ(1, vl.getNonrefPileup(3).stranded_DP[0]);
    ASSERT_EQ(2, vl.getNonrefPileup(3).stranded_DP[1]);
    ASSERT_EQ(0, vl.getNonrefPileup(4).stranded_DP[0]);
    ASSERT_EQ(0, vl.getNonrefPileup(4).stranded_DP[1]);
    ASSERT_EQ(0, vl.getNonrefPileup(5).stranded_DP[0]);
    ASSERT_EQ(0, vl.getNonrefPileup(5).stranded_DP[1]);

    ASSERT_FLOAT_EQ(1 - 1e-6, vl.getNonrefPileup(3).qual_weighted_DP[0]);
    ASSERT_FLOAT_EQ(1 * (1 - 1e-6) + 0.9, vl.getNonrefPileup(3).qual_weighted_DP[1]);
    ASSERT_FLOAT_EQ(3 * (1 - 1e-6) + 0.9, vl.getRefPileup(3).qual_weighted_DP[0]);
    ASSERT_FLOAT_EQ(4 * (1 - 1e-6) + 0.9, vl.getRefPileup(3).qual_weighted_DP[1]);

    ASSERT_EQ(1ull, vl.getVariants().size());
    Variant* v0 = vl.getVariants().front();
    ASSERT_EQ(3, v0->start());
    ASSERT_EQ(3, v0->end());
    ASSERT_EQ("T", v0->alt());

    ASSERT_EQ(4, v0->adr_forward());
    ASSERT_FLOAT_EQ(3 * (1 - 1e-6) + 0.9, v0->wadr_forward());
    ASSERT_EQ(5, v0->adr_backward());
    ASSERT_FLOAT_EQ(4 * (1 - 1e-6) + 0.9, v0->wadr_backward());

    ASSERT_FLOAT_EQ(1 - 1e-6 + 0.9, v0->wada_backward());
    ASSERT_EQ(1, v0->ada_forward());
    ASSERT_FLOAT_EQ(1 - 1e-6, v0->wada_forward());
    ASSERT_EQ(2, v0->ada_backward());
    ASSERT_FLOAT_EQ(1 - 1e-6 + 0.9, v0->wada_backward());

    ASSERT_EQ(0, v0->ado_forward());
    ASSERT_FLOAT_EQ(0, v0->wado_forward());
    ASSERT_EQ(0, v0->ado_backward());
    ASSERT_FLOAT_EQ(0, v0->wado_backward());
}

TEST(VariantCandidateList, CandidateListTestIndel)
{
    VariantCandidateList vl("CCACATATATATATATATATA");

    // add indel at leftmost pos
    RefVar rv;
    rv.start = 3;
    rv.end = 5;
    rv.alt = "C";
    vl.addRefVarObservation(rv, false);

    rv.start = 14;
    rv.end = 17;
    rv.alt = "AT";

    vl.addRefVarObservation(rv, true);

    for (int i = 0; i < 4; ++i)
    {
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[1]);
        ASSERT_EQ(0, vl.getNonrefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getNonrefPileup(i).stranded_DP[1]);
    }

    for (int i = 4; i < (int)vl.getReference().size(); ++i)
    {
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[1]);
        ASSERT_EQ(1, vl.getNonrefPileup(i).stranded_DP[0]);
        ASSERT_EQ(1, vl.getNonrefPileup(i).stranded_DP[1]);
    }

    ASSERT_EQ(1ull, vl.getVariants().size());
    Variant* v0 = vl.getVariants().front();
    ASSERT_EQ(4, v0->start());
    ASSERT_EQ(5, v0->end());
    ASSERT_EQ("", v0->alt());

    ASSERT_FLOAT_EQ(1 - 1e-6, v0->wada_backward());
    ASSERT_EQ(1, v0->ada_forward());
    ASSERT_FLOAT_EQ(1 - 1e-6, v0->wada_backward());
    ASSERT_EQ(1, v0->ada_backward());

    ASSERT_EQ(0, v0->adr_forward());
    ASSERT_FLOAT_EQ(0, v0->wadr_forward());
    ASSERT_EQ(0, v0->adr_backward());
    ASSERT_FLOAT_EQ(0, v0->wadr_backward());

    ASSERT_EQ(0, v0->ado_forward());
    ASSERT_FLOAT_EQ(0, v0->wado_forward());
    ASSERT_EQ(0, v0->ado_backward());
    ASSERT_FLOAT_EQ(0, v0->wado_backward());
}

TEST(VariantCandidateList, CandidateListTestInsertion)
{
    VariantCandidateList vl("CCACATATATATATATATATA");

    // add indel at leftmost pos
    RefVar rv;
    rv.start = 4;
    rv.end = 3;
    rv.alt = "AT";
    vl.addRefVarObservation(rv, false);

    rv.start = 14;
    rv.end = 13;
    rv.alt = "AT";

    vl.addRefVarObservation(rv, true);

    for (int i = 0; i < 4; ++i)
    {
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[1]);
        ASSERT_EQ(0, vl.getNonrefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getNonrefPileup(i).stranded_DP[1]);
    }

    for (int i = 4; i < (int)vl.getReference().size(); ++i)
    {
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[0]);
        ASSERT_EQ(0, vl.getRefPileup(i).stranded_DP[1]);
        ASSERT_EQ(1, vl.getNonrefPileup(i).stranded_DP[0]);
        ASSERT_EQ(1, vl.getNonrefPileup(i).stranded_DP[1]);
    }

    ASSERT_EQ(1ull, vl.getVariants().size());
    Variant* v0 = vl.getVariants().front();
    ASSERT_EQ(4, v0->start());
    ASSERT_EQ(3, v0->end());
    ASSERT_EQ("AT", v0->alt());

    ASSERT_FLOAT_EQ(1 - 1e-6, v0->wada_backward());
    ASSERT_EQ(1, v0->ada_forward());
    ASSERT_FLOAT_EQ(1 - 1e-6, v0->wada_backward());
    ASSERT_EQ(1, v0->ada_backward());

    ASSERT_EQ(0, v0->adr_forward());
    ASSERT_FLOAT_EQ(0, v0->wadr_forward());
    ASSERT_EQ(0, v0->adr_backward());
    ASSERT_FLOAT_EQ(0, v0->wadr_backward());

    ASSERT_EQ(0, v0->ado_forward());
    ASSERT_FLOAT_EQ(0, v0->wado_forward());
    ASSERT_EQ(0, v0->ado_backward());
    ASSERT_FLOAT_EQ(0, v0->wado_backward());
}
