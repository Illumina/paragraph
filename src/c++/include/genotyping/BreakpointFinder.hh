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
 * Functions to find breakpoints in a graph
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include "genotyping/BreakpointStatistics.hh"
#include "graphcore/Graph.hh"

#include <map>

namespace genotyping
{

typedef std::map<std::string, BreakpointStatistics> BreakpointMap;

/**
 * Create a map of all breakpoints in a graph
 * @param wgraph the graph
 * @return a breakpoint map for the graph
 */
BreakpointMap createBreakpointMap(graphtools::Graph const& wgraph);
}
