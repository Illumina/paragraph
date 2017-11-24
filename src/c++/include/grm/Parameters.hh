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

#include <iostream>
#include <json/json.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace grm
{
class Parameters
{
public:
    explicit Parameters(
        bool print_csv, double genotype_error_rate, int min_overlap_bases, int max_read_times, bool use_em)
        : print_csv_(print_csv)
        , genotype_error_rate_(genotype_error_rate)
        , min_overlap_bases_(min_overlap_bases)
        , max_read_times_(max_read_times)
        , use_em_(use_em){};

    void load(
        std::string& paragraph_path, std::string& reference_path, std::string& manifest_path, std::string& output_path);
    const std::string input_path() const { return paragraph_path_; }
    const std::string reference_path() const { return reference_path_; }
    const std::string manifest_path() const { return manifest_path_; }
    const std::string output_path() const { return output_path_; }
    double use_em() const { return use_em_; }
    bool output_as_csv() const { return print_csv_; }
    double genotype_error_rate() const { return genotype_error_rate_; }
    int min_overlap_bases() const { return min_overlap_bases_; }
    int max_read_times() const { return max_read_times_; }

private:
    bool print_csv_;
    double genotype_error_rate_;
    int min_overlap_bases_;
    int max_read_times_;
    double use_em_;
    std::string paragraph_path_;
    std::string reference_path_;
    std::string manifest_path_;
    std::string output_path_;
};
}
