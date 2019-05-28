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
        int threads = 1, int max_reads = 10000, float bad_align_frac = 0.8, bool path_sequence_matching = false,
        bool graph_sequence_matching = true, bool klib_sequence_matching = false, bool kmer_sequence_matching = false,
        int bad_align_uniq_kmer_len = 0, std::string const& alignment_output_folder = "",
        bool infer_read_haplotypes = false)
        : threads_(threads)
        , max_reads_(max_reads)
        , bad_align_frac_(bad_align_frac)
        , path_sequence_matching_(path_sequence_matching)
        , graph_sequence_matching_(graph_sequence_matching)
        , klib_sequence_matching_(klib_sequence_matching)
        , kmer_sequence_matching_(kmer_sequence_matching)
        , bad_align_uniq_kmer_len_(bad_align_uniq_kmer_len)
        , alignment_output_folder_(alignment_output_folder)
        , infer_read_haplotypes_(infer_read_haplotypes)
    {
    }

    int threads() const { return threads_; }
    int max_reads() const { return max_reads_; }
    float bad_align_frac() const { return bad_align_frac_; }
    bool path_sequence_matching() const { return path_sequence_matching_; }
    bool graph_sequence_matching() const { return graph_sequence_matching_; }
    bool klib_sequence_matching() const { return klib_sequence_matching_; }
    bool kmer_sequence_matching() const { return kmer_sequence_matching_; }
    int bad_align_uniq_kmer_len() const { return bad_align_uniq_kmer_len_; }
    std::string const& alignment_output_folder() const { return alignment_output_folder_; }
    bool infer_read_haplotypes() const { return infer_read_haplotypes_; }

private:
    int threads_ = 1;
    int max_reads_ = 10000;
    float bad_align_frac_ = 0.8;
    bool path_sequence_matching_;
    bool graph_sequence_matching_ = true;
    bool klib_sequence_matching_ = false;
    bool kmer_sequence_matching_ = false;
    int bad_align_uniq_kmer_len_ = 0;
    std::string alignment_output_folder_;
    bool infer_read_haplotypes_ = false;
};
}
