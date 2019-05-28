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