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
 * Test reading VG files
 *
 * \file test_serialisation.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <fstream>
#include <iostream>
#include <limits.h>

#include "gtest/gtest.h"

#include "common.hh"
//
// TEST(Graphs, RW)
//{
//    //    using namespace grmpy;
//    //    DeletionGraphFactory factory((g_testenv->getBasePath() + "/share/data/chrQ.fa").c_str());
//    //    grmpy::Graph g;
//    //    factory.make("chrQ", 3, 5, g);
//    //
//    //    std::ofstream tmpf("test.gsq");
//    //    serializeGraph(tmpf, g);
//}