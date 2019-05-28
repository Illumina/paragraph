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

#include "common/StringUtil.hh"

namespace common
{

/**
 * Stores a genomic interval as [start_, end_]. The positions are 0-based
 */
struct Region
{
    Region() {}

    explicit Region(std::string const& _chrom, int64_t _start, int64_t _end)
        : chrom(_chrom)
        , start(_start)
        , end(_end)
    {
    }

    explicit Region(std::string const& str) { stringutil::parsePos(str, chrom, start, end); }

    std::string chrom;
    int64_t start = -1;
    int64_t end = -1;

    /**
     * allow to cast to string
     */
    operator std::string() const { return stringutil::formatPos(chrom, start, end); }

    Region getLeftFlank(int64_t flank_len) const
    {
        Region left_flank;
        left_flank.chrom = chrom;
        left_flank.start = (start - 1) > flank_len ? (start - flank_len - 1) : 0;
        left_flank.end = start - 1;
        return left_flank;
    }

    Region getRightFlank(int64_t flank_len) const
    {
        Region right_flank;
        right_flank.chrom = chrom;
        right_flank.start = end + 1;
        right_flank.end = end + 1 + flank_len;
        return right_flank;
    }

    Region getExtendedRegion(int64_t flank_len) const
    {
        Region extended_region;
        extended_region.chrom = chrom;
        extended_region.start = start > flank_len ? start - flank_len : 0;
        extended_region.end = end + flank_len;
        return extended_region;
    }

    int64_t length() const { return end + 1 - start; }
};
}
