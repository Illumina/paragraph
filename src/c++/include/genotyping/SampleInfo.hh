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
 * @brief data structure
 *
 * @author Sai Chen & Peter Krusche & Egor Dolzhenko
 * @email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <list>
#include <map>
#include <string>
#include <vector>

#include "json/json.h"

namespace genotyping
{

/**
 * Stats and metadata for a single BAM/CRAM
 */
class SampleInfo
{
public:
    SampleInfo() = default;
    /**
     * Getters / setters for input file names
     */
    std::string const& sample_name() const { return sample_name_; }
    void set_sample_name(std::string const& sample_name) { sample_name_ = sample_name; }
    std::string const& filename() const { return filename_; }
    void set_filename(std::string const& filename) { filename_ = filename; }
    std::string const& index_filename() const { return index_filename_; }
    void set_index_filename(std::string const& index_filename) { index_filename_ = index_filename; }

    /**
     * Getters / setters for BAM statistics
     */
    unsigned int read_length() const { return read_length_; }
    void set_read_length(unsigned int read_length) { read_length_ = read_length; }
    double autosome_depth() const { return autosome_depth_; }
    void set_autosome_depth(double autosome_depth);
    double depth_sd() const { return depth_sd_; }
    void set_depth_sd(double depth_sd) { depth_sd_ = depth_sd; }

    /**
     * Getters / setters for sex information
     */
    enum Sex
    {
        UNKNOWN = 0,
        MALE = 1,
        FEMALE = 2,
    };
    void set_sex(std::string sex_string);
    Sex sex() const { return sex_; }

    /**
     * Getter / setter for the alignment data
     */
    void set_alignment_data(Json::Value const& alignment_data) { alignment_data_ = alignment_data; }
    Json::Value const& get_alignment_data() const { return alignment_data_; }

private:
    std::string sample_name_;
    std::string filename_;
    std::string index_filename_;
    unsigned int read_length_ = 0;
    double autosome_depth_ = 0.0;
    double depth_sd_ = 0.0;
    Sex sex_ = UNKNOWN;
    Json::Value alignment_data_ = Json::nullValue;
};

typedef std::vector<SampleInfo> Samples;

/**
 * Load manifest file. A manifest contains
 *
 *  * Sample names
 *  * Sample BAM/CRAM file locations
 *  * Depth estimates, or location of idxdepth output
 *  * (optionally) location of alignment JSON
 *
 * @param filename file name of manifest
 * @return list of sample info records
 */
Samples loadManifest(const std::string& filename);
};
