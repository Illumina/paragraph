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
 * \Summary of graph alignment statistics for validation & filtering purpose
 *
 * \file GraphSummaryStatistics.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Read.hh"
#include "graphcore/Graph.hh"
#include "paragraph/AlignmentStatistics.hh"

#include "json/json.h"

#include <map>
#include <string>
#include <vector>

namespace paragraph
{
/**
 * add summary statistics on graph alignments into JSON output
 */
void summarizeAlignments(graphtools::Graph const& wgraph, common::ReadBuffer const& reads, Json::Value& output);
}
