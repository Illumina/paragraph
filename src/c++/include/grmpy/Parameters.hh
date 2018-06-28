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
