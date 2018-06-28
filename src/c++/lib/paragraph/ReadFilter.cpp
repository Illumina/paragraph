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
 * \summary Read filter chain implementation
 *
 * \file ReadFilter.cpp
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#include "paragraph/ReadFilter.hh"

#include "readfilters/BadAlign.hh"
#include "readfilters/KmerFilter.hh"
#include "readfilters/NonUniq.hh"

namespace paragraph
{

using graphtools::Graph;

/**
 * Read filter interface
 */
class ReadFilterChain : public ReadFilter
{
public:
    void addFilter(std::unique_ptr<ReadFilter> filter) { filters.emplace_back(std::move(filter)); }
    std::pair<bool, std::string> filterRead(common::Read const& read) override
    {
        for (auto& f : filters)
        {
            auto f_result = f->filterRead(read);
            if (f_result.first)
            {
                LOG()->trace("filtered:{}", f_result.second);
                return f_result;
            }
        }
        return { false, "" };
    }

private:
    std::list<std::unique_ptr<ReadFilter>> filters;
};

/**
 * Create parameterized read filter
 * @param graph the graph the read was aligned to
 * @param remove_nonuniq remove reads that don't have a unique best alignment
 * @param bad_align_frac fraction of read that must be aligned in order to pass
 * @param kmer_len length of kmers for uniqueness
 * @return a read filter
 */
std::unique_ptr<ReadFilter>
createReadFilter(Graph const* graph, bool remove_nonuniq, double bad_align_frac, int32_t kmer_len)
{
    std::unique_ptr<ReadFilter> filters(new ReadFilterChain());

    auto chain_ptr = dynamic_cast<ReadFilterChain*>(filters.get());
    if (remove_nonuniq)
    {
        chain_ptr->addFilter(std::unique_ptr<ReadFilter>(new readfilters::NonUniq));
    }
    chain_ptr->addFilter(std::unique_ptr<ReadFilter>(new readfilters::BadAlign(graph, bad_align_frac)));
    if (kmer_len != 0)
    {
        chain_ptr->addFilter(std::unique_ptr<ReadFilter>(new readfilters::KmerFilter(graph, kmer_len)));
    }

    return filters;
}
}
