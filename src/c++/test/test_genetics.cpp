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
 *
 * \file test_genetics.cpp
 * \author Richard Shaw
 * \email rshaw@illumina.com
 *
 */

#include "gtest/gtest.h"

#include "common/Genetics.hh"

using namespace common;

TEST(Genetics, testIsRealBase)
{
    ASSERT_TRUE(isRealBase('A') == true);
    ASSERT_TRUE(isRealBase('C') == true);
    ASSERT_TRUE(isRealBase('G') == true);
    ASSERT_TRUE(isRealBase('T') == true);
    ASSERT_TRUE(isRealBase('N') == false);
    ASSERT_TRUE(isRealBase('D') == false);
    ASSERT_TRUE(isRealBase('E') == false);
}

TEST(Genetics, TransitionBase)
{
    ASSERT_TRUE(transitionBase('A') == 'G');
    ASSERT_TRUE(transitionBase('C') == 'T');
    ASSERT_TRUE(transitionBase('G') == 'A');
    ASSERT_TRUE(transitionBase('T') == 'C');
}

TEST(Genetics, SnvIsTransversion)
{
    bool isValidSnv(false);

    ASSERT_TRUE(snvIsTransversion('A', 'C', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('A', 'T', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('A', 'G', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('A', 'N', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == false);

    ASSERT_TRUE(snvIsTransversion('C', 'A', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('C', 'G', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('C', 'T', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('C', 'N', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == false);

    ASSERT_TRUE(snvIsTransversion('G', 'C', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('G', 'T', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('G', 'A', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('G', 'N', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == false);

    ASSERT_TRUE(snvIsTransversion('T', 'A', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('T', 'G', isValidSnv) == true);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('T', 'C', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == true);

    ASSERT_TRUE(snvIsTransversion('T', 'N', isValidSnv) == false);
    ASSERT_TRUE(isValidSnv == false);
}

TEST(Genetics, ReverseComplement)
{
    ASSERT_EQ(reverseComplement("GGGGaaaaaaaatttatatat"), "atatataaattttttttCCCC");
    ASSERT_EQ(reverseComplement("NNNNnnnn"), "nnnnNNNN");
}
