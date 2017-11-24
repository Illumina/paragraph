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

#include <string>

namespace idxdepth
{
class Parameters
{
public:
    explicit Parameters(std::string const& bam_path, std::string const& reference_path)
        : bam_path_(bam_path)
        , reference_path_(reference_path)
        , threads_(1)
    {
    }

    const std::string& bam_path() const { return bam_path_; }
    const std::string& reference_path() const { return reference_path_; }

    const std::string& include_regex() const { return include_regex_; }
    void set_include_regex(std::string const& ar) { include_regex_ = ar; }
    const std::string& autosome_regex() const { return autosome_regex_; }
    void set_autosome_regex(std::string const& ar) { autosome_regex_ = ar; }
    const std::string& sex_chromosome_regex() const { return sex_chromosome_regex_; }
    void set_sex_chromosome_regex(std::string const& ar) { sex_chromosome_regex_ = ar; }

    int threads() const { return threads_; }
    void set_threads(int threads) { threads_ = threads; }

private:
    std::string bam_path_;
    std::string reference_path_;
    std::string include_regex_;
    std::string autosome_regex_;
    std::string sex_chromosome_regex_;
    int threads_;
};
};
