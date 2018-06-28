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
    void set_autosome_depth(double autosome_depth) { autosome_depth_ = autosome_depth; }

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
