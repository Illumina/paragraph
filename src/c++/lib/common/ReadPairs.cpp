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