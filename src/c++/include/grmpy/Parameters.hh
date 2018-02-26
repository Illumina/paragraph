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

#include <json/json.h>
#include <list>
#include <string>

#include "genotyping/SampleInfo.hh"

namespace grmpy
{
class Parameters
{
public:
    explicit Parameters(
        int sample_threads = 1, int alignment_threads = 1, int max_reads = 10000, float bad_align_frac = 0.8)
        : sample_threads_(sample_threads)
        , alignment_threads_(alignment_threads)
        , max_reads_(max_reads)
        , bad_align_frac_(bad_align_frac)
    {
    }

    void load(
        std::string const& graph_path, std::string const& reference_path, std::string const& manifest_path,
        std::string const& output_path, std::string const& genotyping_parameter_path);

    /**
     * Input/output paths
     */
    int sample_threads() const { return sample_threads_; }
    int alignment_threads() const { return alignment_threads_; }
    int max_reads() const { return max_reads_; }
    float bad_align_frac() const { return bad_align_frac_; }
    const std::string& graph_path() const { return graph_path_; }
    const std::string& reference_path() const { return reference_path_; }
    const std::string& manifest_path() const { return manifest_path_; }
    const std::string& output_path() const { return output_path_; }

    const std::string& genotypingParameterPath() const { return genotyping_parameter_path_; }

    /**
     * Manifest data
     */
    std::list<genotyping::SampleInfo> const& getSamples() const { return samples_; }
    std::list<genotyping::SampleInfo>& getSamples() { return samples_; }

private:
    int sample_threads_ = 1;
    int alignment_threads_ = 1;
    int max_reads_ = 10000;
    float bad_align_frac_ = 0.8;
    std::string graph_path_;
    std::string reference_path_;
    std::string manifest_path_;
    std::string output_path_;

    std::list<genotyping::SampleInfo> samples_;

    std::string genotyping_parameter_path_;
};
}
