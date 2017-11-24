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
 * \brief Representation of a path in a sequence graph
 *
 * \file GraphPath.hh
 * \author Egor Dolzhenko
 * \email edolzhenko@illumina.com
 *
 */

#pragma once

#include <cstdint>
#include <iostream>
#include <list>
#include <string>

#include "graphs/WalkableGraph.hh"

namespace graphs
{
// A path in a sequence graph is given by (1) a sequence of nodes and (2) start/end position on the first/last node. The
// start/end positions are 0-based and form a closed interval.
class GraphPath
{
public:
    // The constructor does not check if the inputs define a well-formed path; isValid() method can be used to do this.
    GraphPath(
        std::shared_ptr<WalkableGraph> wgraph_ptr, int32_t start_position, const std::list<int32_t>& nodes,
        int32_t end_position);
    ~GraphPath();
    GraphPath(const GraphPath& other);
    GraphPath(GraphPath&& other) noexcept;
    GraphPath& operator=(const GraphPath& other);
    GraphPath& operator=(GraphPath&& other) noexcept;
    std::string seq() const;
    bool isValid() const;
    bool operator==(const GraphPath& other) const;
    std::string encode() const;
    int32_t start_position() const;
    int32_t end_position() const;
    GraphPath extendStartPosition(int32_t extension_len) const;
    GraphPath extendEndPosition(int32_t extension_len) const;
    GraphPath extendStartNodeTo(int32_t node_id) const;
    GraphPath extendEndNodeTo(int32_t node_id) const;
    std::list<GraphPath> extendStartBy(int32_t extension_len) const;
    std::list<GraphPath> extendEndBy(int32_t extension_len) const;
    // Computes all possible extensions of the path by the specified length in both directions.
    std::list<GraphPath> extendBy(int32_t start_extension_len, int32_t end_extension_len) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

std::ostream& operator<<(std::ostream& os, const GraphPath& path);
}