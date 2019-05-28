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

#include "common/ReadPairs.hh"
#include "common/Error.hh"
#include "common/Read.hh"
#include "common/ReadPair.hh"

using std::string;
using std::vector;

namespace common
{

void ReadPairs::add(const Read& read)
{
    ReadPair& mates = read_pairs_[read.fragment_id()];
    const int num_initialized_mates_original
        = (int)mates.first_mate().is_initialized() + (int)mates.second_mate().is_initialized();
    mates.add(read);
    const int num_initialized_mates_after_add
        = (int)mates.first_mate().is_initialized() + (int)mates.second_mate().is_initialized();

    num_reads_ += num_initialized_mates_after_add - num_initialized_mates_original;
}

const ReadPair& ReadPairs::operator[](const string& fragment_id) const
{
    if (read_pairs_.find(fragment_id) == read_pairs_.end())
    {
        error("Fragment %s does not exist", fragment_id.c_str());
    }
    return read_pairs_.at(fragment_id);
}

void ReadPairs::getReads(vector<Read>& reads)
{
    for (const auto& kv : read_pairs_)
    {
        ReadPair mates = kv.second;
        if (mates.first_mate().is_initialized())
        {
            reads.emplace_back(mates.first_mate());
        }
        if (mates.second_mate().is_initialized())
        {
            reads.emplace_back(mates.second_mate());
        }
    }
}

void ReadPairs::getReads(vector<p_Read>& reads)
{
    for (const auto& kv : read_pairs_)
    {
        ReadPair mates = kv.second;
        if (mates.first_mate().is_initialized())
        {
            reads.emplace_back(new Read(mates.first_mate()));
        }
        if (mates.second_mate().is_initialized())
        {
            reads.emplace_back(new Read(mates.second_mate()));
        }
    }
}

void ReadPairs::clear()
{
    read_pairs_.clear();
    num_reads_ = 0;
}
}