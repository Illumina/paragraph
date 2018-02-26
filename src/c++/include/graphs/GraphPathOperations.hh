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
 * \brief Various operations on graph paths
 *
 * \file GraphPath.hh
 * \author Egor Dolzhenko
 * \email edolzhenko@illumina.com
 *
 */

#pragma once

#include <string>
#include <vector>

#include "graphs/GraphPath.hh"

namespace graphs
{

/**
 * @brief Splits sequence into segments corresponding to the path
 *
 * @param path Any valid path
 * @param sequence A string having the same length as the path
 * @return Segments of the sequence corresponding to nodes spanned by the path
 */
std::vector<std::string> splitByPath(const GraphPath& path, const std::string& sequence);

/**
 * Return true if two paths overlap either prefix - suffix or suffix-prefix
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths overlap
 */
bool checkPathPrefixSuffixOverlap(GraphPath const& p1, GraphPath const& p2);

/**
 * Paths can be merged if they overlap prefix-suffix / suffix-prefix.
 *
 * @param p1 first path
 * @param p2 second path
 * @return merged path
 */
GraphPath mergePaths(GraphPath const& p1, GraphPath const& p2);

/**
 * Merge a set of paths
 *
 * This will merge paths iteratively until none of the resulting paths overlap
 *
 * @param paths a list of paths
 */
void greedyMerge(std::list<GraphPath>& paths);

/**
 * Merge a set of paths
 *
 * This will merge paths exhaustively, each path is merged with all
 * paths it overlaps until we cannot merge anymore
 *
 * @param paths a list of paths
 */
void exhaustiveMerge(std::list<GraphPath>& paths);
}