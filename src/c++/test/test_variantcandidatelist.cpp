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
