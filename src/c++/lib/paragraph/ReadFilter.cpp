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
