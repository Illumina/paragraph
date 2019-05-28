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
 * \brief Graph caller parameters
 *
 * \file Parameters.hh
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <list>
#include <string>
#include <unordered_map>
#include <vector>

#include <json/json.h>

#include "common/Region.hh"
#include "grm/GraphInput.hh"

namespace paragraph
{

class Parameters
{
public:
    explicit Parameters(
        int max_reads_in = 10000, int min_reads_for_variant = 1, float min_frac_for_variant = 0.0f,
        float bad_align_frac = 0.8, int outputs = output_options::ALL, bool path_sequence_matching = false,
        bool graph_sequence_matching = true, bool klib_sequence_matching = false, bool kmer_sequence_matching = false,
        bool validate_alignments = false)
        : max_reads_(static_cast<size_t>(max_reads_in))
        , min_reads_for_variant_(min_reads_for_variant)
        , min_frac_for_variant_(min_frac_for_variant)
        , bad_align_frac_(bad_align_frac)
        , output_options_(outputs)
        , path_sequence_matching_(path_sequence_matching)
        , graph_sequence_matching_(graph_sequence_matching)
        , klib_sequence_matching_(klib_sequence_matching)
        , kmer_sequence_matching_(kmer_sequence_matching)
        , validate_alignments_(validate_alignments){};

    enum output_options
    {
        ALIGNMENTS = 0x01,
        FILTERED_ALIGNMENTS = 0x02,
        VARIANTS = 0x04,
        NODE_READ_COUNTS = 0x08,
        EDGE_READ_COUNTS = 0x10,
        PATH_READ_COUNTS = 0x20,
        DETAILED_READ_COUNTS = 0x40,
        PATH_COVERAGE = 0x80,
        NODE_COVERAGE = 0x100,
        HAPLOTYPES = 0x200,
        ALL = 0xffffffff
    };

    void load(
        const std::string& graph_path, const std::string& reference_path,
        const std::string& override_target_regions = "");

    const std::string& reference_path() const { return reference_path_; }

    size_t max_reads() const { return max_reads_; }

    int min_reads_for_variant() const { return min_reads_for_variant_; }

    float min_frac_for_variant() const { return min_frac_for_variant_; }

    float bad_align_frac() const { return bad_align_frac_; }

    std::list<common::Region> target_regions() const { return target_regions_; }

    const Json::Value& description() const { return description_; }

    bool output_enabled(output_options const o) const { return static_cast<const bool>((output_options_ & o) != 0); }

    bool path_sequence_matching() const { return path_sequence_matching_; }
    bool graph_sequence_matching() const { return graph_sequence_matching_; }
    bool klib_sequence_matching() const { return klib_sequence_matching_; }
    bool kmer_sequence_matching() const { return kmer_sequence_matching_; }
    bool validate_alignments() const { return validate_alignments_; }
    unsigned longest_alt_insertion() const { return longest_alt_insertion_; }

    uint32_t threads() const { return threads_; }
    void set_threads(uint32_t threads) { threads_ = threads; }

    int kmer_len() const { return kmer_len_; }
    void set_kmer_len(int kmer_len) { kmer_len_ = kmer_len; }

    bool remove_nonuniq_reads() const { return remove_nonuniq_reads_; }
    void set_remove_nonuniq_reads(bool remove_nonuniq_reads) { remove_nonuniq_reads_ = remove_nonuniq_reads; }

private:
    std::string reference_path_;

    size_t max_reads_; ///< maximum number of reads to process per locus

    int min_reads_for_variant_; ///< minimum number of reads required to report a variant
    float min_frac_for_variant_; ///< minimum fraction of reads required to report a variant
    float bad_align_frac_; ///< fraction of read length the alignment score must exceed in order for read to be used

    int output_options_; ///< output options

    bool path_sequence_matching_; ///< enable use of PathAligner
    bool graph_sequence_matching_; ///< enable use of GraphAligner
    bool klib_sequence_matching_; ///< enable use of KlibAligner
    bool kmer_sequence_matching_; ///< enable use of KmerAligner
    bool validate_alignments_;

    Json::Value description_; ///< graph description
    /// if graph contains long insertions, we might want to read mates that are not in the target region.
    unsigned longest_alt_insertion_ = 0;

    std::list<common::Region> target_regions_; ///< target regions for read retrieval

    uint32_t threads_{ 1 }; ///< number of threads for parallel read processing

    int kmer_len_{ 0 }; ///< kmer length for validation

    bool remove_nonuniq_reads_{ true }; // remove reads with no unique alignment
};
}
