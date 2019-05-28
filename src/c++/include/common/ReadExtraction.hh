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
    BamReader& reader, std::list<Region> const& target_regions, int max_num_reads, unsigned longest_alt_insertion,
    std::vector<p_Read>& all_reads, int avr_fragment_length = 333);

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
    std::list<Region> const& target_regions, int max_num_reads, unsigned longest_alt_insertion,
    std::vector<p_Read>& all_reads, int avr_fragment_length = 333);

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
    unsigned longest_alt_insertion, int avr_fragment_length);

/**
 * Low-level read extraction for mapped reads in target region
 * @param read_pairs Container for extracted reads
 * @param reader Reader that will provide the reads
 * @param region Region to check if a read is in
 * @return average read length
 */
int extractMappedReadsFromRegion(ReadPairs& read_pairs, int max_num_reads, ReadReader& reader, const Region& region);

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
