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

#include "common/ReadPairs.hh"
#include "common/ReadReader.hh"
#include "common/Region.hh"

#include <list>
#include <utility>
#include <vector>

namespace common
{

/**
 * High-level read extraction interface
 * @param reader An open reader
 * @param target_regions list of target regions
 * @param max_reads maximum number of reads per target region to retrieve
 * @param all_reads output vector to store retrieved reads
 * @param avr_fragment_length decides how long to extend beyond target region
 */
void extractReads(
    BamReader& reader, std::list<Region> const& target_regions, int max_num_reads, std::vector<p_Read>& all_reads,
    int avr_fragment_length = 333);

/**
 * High-level read extraction interface
 * @param bam_path path to BAM file
 * @param bam_index_path path to BAM index file (or empty to use default location)
 * @param reference_path path to FASTA reference file
 * @param target_regions list of target regions
 * @param max_reads maximum number of reads per target region to retrieve
 * @param all_reads output vector to store retrieved reads
 * @param avr_fragment_length decides how long to extend beyond target region
 */
void extractReads(
    const std::string& bam_path, const std::string& bam_index_path, const std::string& reference_path,
    std::list<Region> const& target_regions, int max_num_reads, std::vector<p_Read>& all_reads,
    int avr_fragment_length = 333);

/**
 * Lower-level read extraction interface for specified target region
 * @return <num_original_extracted, num_recovered_mates> when finish
 * @param output Vector to store retrieved uniq ptr of reads
 * @param max_reads Max number of reads to load
 * @param reader Reader that will provide the reads
 * @param region Target region
 * @param avr_fragment_length Decides how long to extend beyond target region
 */
std::pair<int, int> extractReadsFromRegion(
    std::vector<p_Read>& all_reads, int max_num_reads, ReadReader& reader, const Region& region,
    int avr_fragment_length);

/**
 * Low-level read extraction for mapped reads in target region
 * @param read_pairs Container for extracted reads
 * @param reader Reader that will provide the reads
 * @param region Region to check if a read is in
 */
void extractMappedReadsFromRegion(ReadPairs& read_pairs, int max_num_reads, ReadReader& reader, const Region& region);

/**
 * return true if this aligned read or its mate overlaps >= 1 base with the target region
 */
bool isReadOrItsMateInRegion(Read& read, const Region& region);

/**
 * Lower-level read extraction interface -- recover mates
 * @param reader ReadReader to find mates in
 * @param read_pairs ReadPairs data structure to update
 */
void recoverMissingMates(ReadReader& reader, ReadPairs& read_pairs);
};
