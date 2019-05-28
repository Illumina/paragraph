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
 *  \brief Phred conversion
 *
 * \file Phred.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Error.hh"
#include "gtest/gtest.h"

#include "common/Phred.hh"

using namespace common::phred;

TEST(Phred, ScaleConversion)
{
    ASSERT_DOUBLE_EQ(60, errorProbToPhred(1e-6));
    ASSERT_NEAR(1e-6, phredToErrorProb(60), 1e-10);
    ASSERT_DOUBLE_EQ(-2, phredToLogErrorProb(20));
    ASSERT_DOUBLE_EQ(10, logErrorProbToPhred(-1));
}
