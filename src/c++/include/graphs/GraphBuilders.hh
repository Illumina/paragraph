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

#pragma once

#include <list>
#include <string>
#include <vector>

#include "graphs/Graph.hh"

namespace graphs
{

struct StrSpecUnit
{
    enum class UnitType
    {
        STR,
        SEQ
    };
    static std::map<UnitType, const char*> unit_type_to_string;
    UnitType type;
    std::string seq;
};

/**
 * Create a deletion graph consisting of three nodes
 *
 * @param left_flank left flanking reference sequence
 * @param deletion  deleted reference sequence
 * @param right_flank right flanking reference sequence
 * @param g output graph structure
 */
void makeDeletionGraph(
    const std::string& left_flank, const std::string& deletion, const std::string& right_flank, Graph& g);

/**
 * Create a swap graph consisting of four nodes.
 *
 * @param left_flank: left flanking reference sequence
 * @param deletion: sequence deleted from the reference
 * @param insertion: sequences inserted into the reference
 * @param right_flank: right flanking reference sequence
 * @param g output graph structure
 */
void makeSimpleSwapGraph(
    const std::string& left_flank, const std::string& deletion, const std::string& insertion,
    const std::string& right_flank, Graph& g);

/**
 * Create a double swap graph with 7 nodes.
 *
 * @param left_flank left flanking sequence (node id = 0)
 * @param deletion1 first deletion sequence (node id = 1)
 * @param insertion1 first insertion sequence (node id = 2)
 * @param middle middle sequence (node id = 3)
 * @param deletion2 second deletion sequence (node id = 4)
 * @param insertion2 second insertion sequence (node id = 5)
 * @param right_flank right flanking reference sequence (node id = 6)
 * @param g output graph structure
 */
void makeDoubleSwapGraph(
    const std::string& left_flank, const std::string& deletion1, const std::string& insertion1,
    const std::string& middle, const std::string& deletion2, const std::string& insertion2,
    const std::string& right_flank, Graph& g);

/**
 * Create a swap graph consisting of four nodes.
 *
 * @param left_flank left flanking reference sequence
 * @param deletion deleted reference sequence
 * @param alternatives alternate sequences
 * @param right_flank right flanking reference sequence
 * @param g output graph structure
 */
void makeSwapGraph(
    const std::string& left_flank, const std::string& deletion, const std::list<std::string>& alternatives,
    const std::string& right_flank, Graph& g);

/**
 * Create a multi-unit STR graph.
 *
 * @param left_flank left flanking reference sequence
 * @param right_flank right flanking reference sequence
 * @param g output graph structure
 */
void makeStrGraph(
    const std::string& left_flank, const std::string& right_flank, size_t read_len, std::vector<StrSpecUnit> spec,
    Graph& g);
}
