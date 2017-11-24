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
 * Brief data structure to extract and store bam stats from manifest
 *
 *
 * \author Sai Chen
 * \email schen6@illumina.com
 *
 */

#pragma once
#include <map>
#include <string>
#include <vector>

namespace genotyping
{
/**
 * stats for a single bam/cram
 */
struct SingleIdxdepth
{
    int32_t read_length;
    double autosome_depth;
};

class Idxdepth
{
public:
    Idxdepth(int read_length = 0, int autosome_depth = 0, int sample_size = 0);
    /**
     * load idxdepth information for all the samples in paragraph Json
     * all info stored in idxdepth_info
     */
    void load(const std::string& manifest_path);

    /**
     * Simple getters
     */
    size_t sampleSize() { return sample_names_.size(); }
    std::string getSampleName(int index) { return sample_names_[index]; }
    SingleIdxdepth getIdxStats(int index) { return idx_stats_[index]; }

private:
    /**
     * constants for reading manifest
     */
    const size_t num_tokens_ = 4;
    const int32_t depth_token_index_ = 2;
    const int32_t read_len_token_index_ = 3;
    /**
     * Store idxdepth information for all samples
     */
    std::vector<std::string> sample_names_;
    std::vector<SingleIdxdepth> idx_stats_;

    // support functions
    void is_good_manifest_header(std::string& line);
};
};
