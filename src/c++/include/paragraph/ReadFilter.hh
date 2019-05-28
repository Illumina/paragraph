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
 * \summary Read filtering interface
 *
 * \file ReadFilter.hh
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & edolzhenko@illumina.com & pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Read.hh"
#include "graphcore/Graph.hh"

#include <map>
#include <memory>
#include <string>

namespace paragraph
{

class ReadFilter
{
public:
    virtual ~ReadFilter() = default;
    /**
     * Create a read filter
     * @param read the read
     * @return { true if filtered, reason for filtering or "" }
     */
    virtual std::pair<bool, std::string> filterRead(common::Read const& read) = 0;
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
createReadFilter(graphtools::Graph const* graph, bool remove_nonuniq, double bad_align_frac, int32_t kmer_len = 0);
}
