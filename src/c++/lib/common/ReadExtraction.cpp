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

#include "common/ReadExtraction.hh"
#include "common/Error.hh"

#include <cstdlib>
#include <list>

namespace common
{

/**
 * Wrapper function to retrieve reads from target regions
 * @param bam_path path to BAM file
 * @param target_regions list of target regions
 * @param max_reads maximum number of reads per target region to retrieve
 * @param all_reads output vector for the reads we retrieve
 */
void extractReads(int max_reads, ReadReader& reader, ReadPairs& read_pairs)
{
    Read read;
    while ((read_pairs.num_reads() != max_reads) && reader.getAlign(read))
    {
        read_pairs.add(read);
    }
}

/**
 * Lower-level read extraction interface
 * @param max_reads maximum number of reads to load
 * @param reader Reader that will provide the reads
 * @param read_pairs output ReadPairs structure
 */
void recoverMissingMates(ReadReader& reader, ReadPairs& read_pairs)
{
    for (const auto& kv : read_pairs)
    {
        const ReadPair& read_pair = kv.second;

        // If a mate is missing, try to recover it.
        if (!read_pair.first_mate().is_initialized() || !read_pair.second_mate().is_initialized())
        {
            // Get references for existing Align and the one that
            // needs to be filled in and then process them below. This
            // is done to avoid code duplication.
            const Read& initialized_read
                = read_pair.first_mate().is_initialized() ? read_pair.first_mate() : read_pair.second_mate();

            // Do not recover nearby mates.
            const int kMaxNormalDistanceBetweenMates = 1000;
            if ((initialized_read.chrom_id() == initialized_read.mate_chrom_id())
                && (std::abs(initialized_read.pos() - initialized_read.mate_pos()) < kMaxNormalDistanceBetweenMates))
            {
                continue;
            }

            Read missing_read;
            reader.getAlignedMate(initialized_read, missing_read);
            if (missing_read.is_initialized())
            {
                read_pairs.add(missing_read);
            }
        }
    }
}

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
    int max_reads, std::vector<p_Read>& all_reads)
{
    auto logger = LOG();
    logger->info("Retrieving reads from {}", bam_path);

    for (const auto& region : target_regions)
    {
        logger->info("[Initializing BAM reader for region {}]", (std::string)region);
        BamReader reader(bam_path, reference_path);
        reader.setRegion(region);
        logger->info("[Done initializing BAM reader]");

        ReadPairs read_pairs;
        // limit the maximum number of reads we read
        extractReads(max_reads, reader, read_pairs);
        if (max_reads == read_pairs.num_reads())
        {
            logger->warn("Reached maximum number of reads ({}) in region {}", max_reads, (std::string)region);
        }
        else
        {
            const int num_reads_original = read_pairs.num_reads();
            recoverMissingMates(reader, read_pairs);
            const int num_reads_recovered = read_pairs.num_reads() - num_reads_original;
            logger->info("[Recovered {} + {} additional reads]", num_reads_original, num_reads_recovered);
        }

        read_pairs.getReads(all_reads);
    }
    logger->info("Done retrieving reads from {}", bam_path);
}
}
