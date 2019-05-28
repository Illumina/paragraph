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
