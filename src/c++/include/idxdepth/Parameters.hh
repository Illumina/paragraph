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

#include <string>
#include <utility>

namespace idxdepth
{
class Parameters
{
public:
    explicit Parameters(std::string bam_path, std::string bam_index_path, std::string reference_path, bool altcontig)
        : bam_path_(std::move(bam_path))
        , bam_index_path_(std::move(bam_index_path))
        , reference_path_(std::move(reference_path))
        , alt_contig_(altcontig)
        , threads_(1)
    {
    }

    const std::string& bam_path() const { return bam_path_; }
    const std::string& bam_index_path() const { return bam_index_path_; }
    const std::string& reference_path() const { return reference_path_; }
    bool include_alt_contig() const { return alt_contig_; }

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
    std::string bam_index_path_;
    std::string reference_path_;
    bool alt_contig_;
    std::string include_regex_;
    std::string autosome_regex_;
    std::string sex_chromosome_regex_;
    int threads_;
};
};
