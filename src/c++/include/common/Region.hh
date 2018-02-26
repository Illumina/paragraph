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
