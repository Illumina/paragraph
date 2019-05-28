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
 * \brief Graph class to collect all information for a single graph, path helpers, ...
 *
 * \file Graph.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

#include "json/json.h"

namespace grm
{

/**
 * Initialize graph from JSON
 * @param in Input JSON node
 * @param reference reference fasta name
 * @param store_ref_sequence store sequence for reference nodes
 */
graphtools::Graph graphFromJson(Json::Value const& in, std::string const& reference, bool store_ref_sequence = true);

/**
 * Read paths from JSON
 * @param graph graph to use for the paths
 * @param in_paths Input JSON node with paths
 */
std::list<graphtools::Path> pathsFromJson(graphtools::Graph const* graph, Json::Value const& in_paths);
};
