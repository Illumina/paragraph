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

#include "idxdepth/IndexBinning.hh"

#include <functional>
#include <htslib/sam.h>
#include <memory>

extern "C" {
#include <htslib/hts.h>
#include <htslib/sam.h>
};

#include "common/Error.hh"
#include "common/StringUtil.hh"

namespace idxdepth
{

typedef std::vector<IndexBin> PerChrBins;

/**
 * Get index bins
 * @param bam BAM / CRAM file
 * @param output output vector for bins
 */
void getIndexBins(std::string const& bam_path, std::vector<IndexBin>& output)
{
    auto logger = LOG();
    // Open alignment file for reading.
    auto hts_file_ptr
        = std::unique_ptr<htsFile, std::function<void(htsFile*)>>{ sam_open(bam_path.c_str(), "r"), hts_close };

    if (!hts_file_ptr)
    {
        error("ERROR: Failed to open %s", bam_path.c_str());
    }

    if (hts_file_ptr->format.format == bam)
    {
        logger->debug("\t[Input {} is a BAM file]", bam_path);
        error("BAM index reading not available yet.");
    }
    else if (hts_file_ptr->format.format == cram)
    {
        logger->debug("\t[Input {} is a CRAM file]", bam_path);
    }
    else
    {
        error("ERROR: Unknown alignment file format.");
    }

    // Read BAM header
    auto header = std::unique_ptr<bam_hdr_t, std::function<void(bam_hdr_t*)>>{ sam_hdr_read(hts_file_ptr.get()),
                                                                               bam_hdr_destroy };
    if (!header)
    {
        error("ERROR: Failed to read header of %s", bam_path.c_str());
    }

    // 16kb binning
    std::vector<PerChrBins> bins;
    bins.resize(static_cast<unsigned long>(header->n_targets));
    for (int i = 0; i < header->n_targets; ++i)
    {
        bins[i].resize((header->target_len[i] >> 14) + 1);
    }

    // gets bin index given a reference ID and start position
    auto CanonicalBin = [&bins](int rid, uint64_t start) -> PerChrBins::iterator {
        assert(rid < (int)bins.size());
        const uint64_t canonical_start = start >> 14;
        assert(canonical_start < (uint64_t)bins[rid].size());
        auto it = bins[rid].begin();
        std::advance(it, (long)canonical_start);
        return it;
    };

    auto UpdateBin
        = [&bins](int rid, PerChrBins::iterator bin_it, uint64_t idx_start, uint64_t idx_end, size_t bytes) -> bool {
        assert(rid < (int)bins.size());
        assert(idx_start <= idx_end);
        if (bin_it == bins[rid].end())
        {
            return false;
        }

        // calculate bin start / end
        uint64_t bin_start = static_cast<size_t>(std::distance(bins[rid].begin(), bin_it));
        bin_start <<= 14;
        uint64_t bin_end = bin_start + 16383;

        bin_it->start = bin_start;
        bin_it->end = bin_end;

        // if the index ends before this bin starts or begins after this bin ends,
        // its bytes do not contribute to this bin's bytes
        if (idx_start > bin_end || idx_end < bin_start)
        {
            return false;
        }

        bin_it->slices++;
        bin_it->overlapping_bytes += bytes;
        uint64_t overlap = std::min(bin_end, idx_end) - std::max(bin_start, idx_start) + 1;

        // adjusted bytes is the adjusted by the percent the index entry overlaps by the bin and
        // normalized to report the average number of bytes per position
        bin_it->adjusted_bytes += ((double)bytes * overlap) / (idx_end - idx_start + 1) / 16384;

        return true;
    };
    size_t unaligned_bytes = 0;

    auto indexfile
        = std::unique_ptr<htsFile, std::function<void(htsFile*)>>{ hts_open((bam_path + ".crai").c_str(), "rz"),
                                                                   hts_close };
    assert(indexfile);

    int res = 1;
    while (res > 0)
    {
        kstring_t str;
        str.l = 0;
        str.m = 0;
        str.s = nullptr;
        res = hts_getline(indexfile.get(), '\n', &str);
        if (res > 0)
        {
            std::vector<std::string> fields;
            common::stringutil::split(str.s, fields, "\t");
            // CRAM index fields are:
            // (0) sequence ID, i.e. contig ID; (1) alignment start; (2) alignment span;
            // (3) container start byte offset in file; (4) slice byte offset in container; (5) bytes in slice

            assert(fields.size() == 6);

            const int rid = std::stoi(fields[0]);
            const size_t bytes = std::stoull(fields[5]);

            if (rid < 0)
            {
                unaligned_bytes += bytes;
                continue;
            }

            const uint64_t start = std::stoull(fields[1]);
            const uint64_t end = start + std::stoull(fields[2]) - 1;

            auto start_slice_it = CanonicalBin(rid, start);
            while (UpdateBin(rid, start_slice_it, start, end, bytes))
            {
                ++start_slice_it;
            }
        }

        if (str.s != nullptr)
        {
            free(str.s);
        }
    }

    int bin_id = 0;
    double mean_adjusted_bytes = 0;
    for (int rid = 0; rid < header->n_targets; ++rid)
    {
        for (auto& input_bin : bins[rid])
        {
            input_bin.bin_id = bin_id++;
            mean_adjusted_bytes += input_bin.adjusted_bytes;
        }
    }
    if (bin_id > 0)
    {
        mean_adjusted_bytes /= bin_id;
    }
    for (int rid = 0; rid < header->n_targets; ++rid)
    {

        for (auto const& input_bin : bins[rid])
        {
            output.push_back(input_bin);
            output.back().chrom = header->target_name[rid];
            output.back().normalized_depth = output.back().adjusted_bytes / mean_adjusted_bytes;
        }
    }
}
}
