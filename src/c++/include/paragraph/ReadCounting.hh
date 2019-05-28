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

#pragma once

#include "Parameters.hh"
#include "common/Read.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphCoordinates.hh"
#include "json/json.h"

namespace paragraph
{

/**
 * Output disambiguated reads counts for graph elements
 * @param coordinates graph and coordinate information
 * @param reads list of reads
 * @param output output JSON node
 * @param by_node output per-node counts
 * @param by_edge output per-edge counts
 * @param by_pathFam output per-pathFamily counts
 * @param pathFam_detailed output node and edge counts for each path-family
 */
void countReads(
    graphtools::GraphCoordinates const& coordinates, common::ReadBuffer const& reads, Json::Value& output,
    bool by_node = true, bool by_edge = true, bool by_pathFam = true, bool pathFam_detailed = false);
}
