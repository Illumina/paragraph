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
 * \brief common input parameters for pam and paragraph
 *
 * \file Parameters.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#pragma once

#include <list>
#include <string>

#include "common/Region.hh"
#include <json/json.h>

namespace pam
{
class Parameters
{
public:
    explicit Parameters(int max_reads_in = 10000, int outputs = output_options::ALL)
        : max_reads_(static_cast<size_t>(max_reads_in))
        , output_options_(outputs){};

    enum output_options
    {
        ALIGNMENTS = 0x01,
        FILTERED_ALIGNMENTS = 0x02,
        NODE_READ_COUNTS = 0x04,
        NODE_COVERAGE = 0x08,
        ALL = 0xffffffff
    };

    void load(
        const std::string& bam_path, const std::string& event_path, const std::string& reference_path,
        const std::string& override_target_regions = "");

    const std::string& reference_path() const { return reference_path_; }

    const std::string& bam_path() const { return bam_path_; }

    size_t max_reads() const { return max_reads_; }

    std::list<common::Region> target_regions() const { return target_regions_; }

    const Json::Value& description() const { return description_; }

    bool output_enabled(output_options const o) const { return static_cast<const bool>((output_options_ & o) != 0); }

private:
    std::string reference_path_;
    std::string bam_path_;

    size_t max_reads_; ///< maximum number of reads to process per locus

    int output_options_; ///< output options

    Json::Value description_; ///< target description

    std::list<common::Region> target_regions_; ///< target regions for read retrieval
};
}