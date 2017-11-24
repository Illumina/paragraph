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
#include <vector>

namespace common
{

/**
 * retrieve reads from target regions
 * @param bam_path path to BAM file
 * @param reference_path path to FASTA reference file
 * @param target_regions list of target regions
 * @param max_reads maximum number of reads per target region to retrieve
 * @param all_reads output vector for the reads we retrieve
 */
void extractReads(
    const std::string& bam_path, const std::string& reference_path, std::list<Region> const& target_regions,
    int max_reads, std::vector<p_Read>& all_reads);

/**
 * Lower-level read extraction interface
 * @param max_reads maximum number of reads to load
 * @param reader Reader that will provide the reads
 * @param read_pairs output ReadPairs structure
 */
void extractReads(int max_reads, ReadReader& reader, ReadPairs& read_pairs);

/**
 * Lower-level read extraction interface -- recover mates
 * @param reader ReadReader to find mates in
 * @param read_pairs ReadPairs data structure to update
 */
void recoverMissingMates(ReadReader& reader, ReadPairs& read_pairs);
};
