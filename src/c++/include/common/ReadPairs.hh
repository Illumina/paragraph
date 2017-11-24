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

#include <array>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/BamReader.hh"
#include "common/Read.hh"
#include "common/ReadPair.hh"

namespace common
{

/**
 * Read pair container class
 */
class ReadPairs
{
public:
    typedef std::map<std::string, ReadPair>::const_iterator const_iterator;

    const_iterator begin() const { return read_pairs_.begin(); }

    const_iterator end() const { return read_pairs_.end(); }

    ReadPairs() = default;

    void clear();

    void add(const Read& read);

    const ReadPair& operator[](const std::string& fragment_id) const;

    int num_reads() const { return num_reads_; }

    void getReads(std::vector<Read>& reads);
    void getReads(std::vector<p_Read>& reads);

private:
    std::map<std::string, ReadPair> read_pairs_;
    int num_reads_ = 0;
};
}
