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
 *  \brief Test main file
 *
 * \file test_main.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <libgen.h>
#include <limits.h>
#include <stdlib.h>
#include <string>

#include "common.hh"
#include "gtest/gtest.h"

GTestEnvironment* g_testenv = nullptr;

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // determine base path
    char* bn = dirname(argv[0]);
    char actualpath[PATH_MAX + 1];
    if (!realpath(bn, actualpath))
    {
        exit(1);
    }

    ::testing::AddGlobalTestEnvironment(g_testenv = new GTestEnvironment(actualpath));
    int res = RUN_ALL_TESTS();

    return res;
}
