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

#include "common/BamReader.hh"

#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <htslib/sam.h>

extern "C" {
#include <htslib/hts.h>
#include <htslib/sam.h>
};

#include "common/Error.hh"
#include "common/Fasta.hh"
#include "common/ReadPileup.hh"
#include "common/StringUtil.hh"
#include "spdlog/spdlog.h"

#include "htslib/hts.h"
#include "htslib/sam.h"

using std::string;
using std::vector;

namespace common
{

static inline void decodeHtsBases(bam1_t* hts_align_ptr, string& bases)
{
    uint8_t* hts_seq_ptr = bam_get_seq(hts_align_ptr);
    const int32_t read_len = hts_align_ptr->core.l_qseq;
    bases.resize((unsigned long)read_len);

    for (int32_t i = 0; i < read_len; ++i)
    {
        bases[i] = seq_nt16_str[bam_seqi(hts_seq_ptr, i)];
    }
}

static inline void decodeHtsQuals(bam1_t* hts_align_ptr, string& quals)
{
    uint8_t* hts_quals_ptr = bam_get_qual(hts_align_ptr);
    const int32_t read_len = hts_align_ptr->core.l_qseq;
    quals.resize((unsigned long)read_len);

    uint8_t* test_hts_quals_ptr = hts_quals_ptr;

    for (int32_t i = 0; i < read_len; ++i)
    {
        quals[i] = static_cast<char>(33 + test_hts_quals_ptr[i]);
    }
}

/**
 * Decode BAM alignment from HTSLib struct.
 *
 * @param hts_align_ptr a bam1_t * to initialize from.
 * Passed as void* to avoid dependency on htslib headers
 *
 * TODO we could make this a constructor
 */
static inline void decodeHtsAlign(void* _hts_align_ptr, Read& read)
{
    auto* hts_align_ptr = (bam1_t*)_hts_align_ptr;
    const string fragment_id = bam_get_qname(hts_align_ptr);
    string bases, quals;
    decodeHtsBases(hts_align_ptr, bases);
    decodeHtsQuals(hts_align_ptr, quals);
    read.set_fragment_id(fragment_id);
    read.set_bases(bases);
    read.set_quals(quals);

    const auto& flag = hts_align_ptr->core.flag;
    read.set_is_mapped((flag & BamReader::kIsMapped) == 0);
    read.set_is_first_mate((flag & BamReader::kIsFirstMate) != 0);
    read.set_is_mate_mapped((flag & BamReader::kIsMateMapped) == 0);
    read.set_is_reverse_strand(bam_is_rev(hts_align_ptr));
    read.set_is_mate_reverse_strand(bam_is_mrev(hts_align_ptr));

    read.set_chrom_id(hts_align_ptr->core.tid);
    read.set_pos(hts_align_ptr->core.pos);
    read.set_mapq(hts_align_ptr->core.qual);
    read.set_mate_chrom_id(hts_align_ptr->core.mtid);
    read.set_mate_pos(hts_align_ptr->core.mpos);
}

struct BamReader::BamReaderImpl
{
    BamReaderImpl() = default;
    BamReaderImpl(BamReaderImpl const&) = delete;
    BamReaderImpl(BamReaderImpl&&) = delete;
    ~BamReaderImpl() { close(); }
    BamReaderImpl& operator=(BamReaderImpl const&) = delete;
    BamReaderImpl& operator=(BamReaderImpl&&) = delete;

    void open(const std::string& path, const std::string& reference)
    {
        auto logger = LOG();
        if (hts_file_ptr_ != nullptr || hts_bam_hdr_ptr_ != nullptr || hts_idx_ptr_ != nullptr
            || hts_bam_align_ptr_ != nullptr)
        {
            close();
        }

        // Open alignment file file for reading.
        hts_file_ptr_ = sam_open(path.c_str(), "r");

        if (hts_file_ptr_ == nullptr)
        {
            error("ERROR: Failed to open %s", path.c_str());
        }

        if (hts_file_ptr_->format.format == bam)
        {
            logger->debug("\t[Input {} is a BAM file]", path);
        }
        else if (hts_file_ptr_->format.format == cram)
        {
            logger->debug("\t[Input {} is a CRAM file]", path);
        }
        else
        {
            error("ERROR: Unknown alignment file format.");
        }

        assert(hts_set_fai_filename(hts_file_ptr_, (reference + ".fai").c_str()) == 0);

        // Read BAM header
        hts_bam_hdr_ptr_ = sam_hdr_read(hts_file_ptr_);

        if (hts_bam_hdr_ptr_ == nullptr)
        {
            error("ERROR: Failed to read header of %s", path.c_str());
        }

        for (int n = 0; n < hts_bam_hdr_ptr_->n_targets; ++n)
        {
            header_contig_map[hts_bam_hdr_ptr_->target_name[n]] = n;
        }

        // Load the index
        hts_idx_ptr_ = sam_index_load(hts_file_ptr_, path.c_str());

        if (hts_idx_ptr_ == nullptr)
        {
            error("ERROR: Failed to read index of %s", path.c_str());
        }

        file_path = path;
        at_file_end_ = false;
    }

    void close()
    {
        if (hts_itr_ptr_ != nullptr)
        {
            hts_itr_destroy(hts_itr_ptr_);
            hts_itr_ptr_ = nullptr;
        }

        if (hts_bam_align_ptr_ != nullptr)
        {
            bam_destroy1(hts_bam_align_ptr_);
            hts_bam_align_ptr_ = nullptr;
        }

        if (hts_idx_ptr_ != nullptr)
        {
            hts_idx_destroy(hts_idx_ptr_);
            hts_idx_ptr_ = nullptr;
        }

        if (hts_bam_hdr_ptr_ != nullptr)
        {
            bam_hdr_destroy(hts_bam_hdr_ptr_);
            hts_bam_hdr_ptr_ = nullptr;
        }

        if (hts_file_ptr_ != nullptr)
        {
            sam_close(hts_file_ptr_);
            hts_file_ptr_ = nullptr;
        }
    }

    std::string file_path;
    std::string reference_path;

    // Pointer to the input BAM/CRAM file itself.
    htsFile* hts_file_ptr_ = nullptr;
    // A pointer to BAM header.
    bam_hdr_t* hts_bam_hdr_ptr_ = nullptr;
    // A pointer to BAM index.
    hts_idx_t* hts_idx_ptr_ = nullptr;
    // Pointer to the target region.
    hts_itr_t* hts_itr_ptr_ = nullptr;
    // A pointer to an alignment in the BAM file.
    bam1_t* hts_bam_align_ptr_ = nullptr;
    bool at_file_end_ = false;

    std::unordered_map<std::string, int> header_contig_map;
};

/**
 * Create BAM / CRAM reader. Reference is required to match the
 * reference FASTA file for the BAM/CRAM
 * @param path file path
 * @param reference path to FASTA reference
 */
BamReader::BamReader(const std::string& path, const std::string& reference)
    : _impl(new BamReaderImpl())
{
    assertFileExists(path);
    assertFileExists(reference);
    assertFileExists(reference + ".fai");
    _impl->file_path = path;
    _impl->reference_path = reference;
}

BamReader::BamReader(BamReader&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

BamReader& BamReader::operator=(BamReader&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

BamReader::~BamReader() = default;

void BamReader::setRegion(const std::string& region_encoding)
{
    _impl->open(_impl->file_path, _impl->reference_path);
    _impl->hts_itr_ptr_ = sam_itr_querys(_impl->hts_idx_ptr_, _impl->hts_bam_hdr_ptr_, region_encoding.c_str());
    if (_impl->hts_itr_ptr_ == nullptr)
    {
        error("Failed to jump to %s in %s", region_encoding.c_str(), _impl->file_path.c_str());
    }
}

bool BamReader::getAlign(Read& read)
{
    if (_impl->hts_file_ptr_ == nullptr)
    {
        throw std::logic_error("Error: BAM file is not open.");
    }

    if (_impl->at_file_end_)
    {
        return false;
    }

    if (_impl->hts_bam_align_ptr_ == nullptr)
    {
        _impl->hts_bam_align_ptr_ = bam_init1();
    }

    int read_ret = 0;
    assert(_impl->hts_itr_ptr_);
    read_ret = SkipToNextGoodAlign();

    if (read_ret == -1)
    {
        _impl->at_file_end_ = true; // EOF
        return false;
    }

    if (read_ret < -1)
    {
        error("ERROR: Failed to extract read from BAM.");
    }

    decodeHtsAlign(_impl->hts_bam_align_ptr_, read);

    return true;
}

int BamReader::SkipToNextGoodAlign()
{
    bool is_primary_align = false;
    int return_value = 0;

    while (!is_primary_align)
    {
        return_value = sam_itr_next(_impl->hts_file_ptr_, _impl->hts_itr_ptr_, _impl->hts_bam_align_ptr_);
        if (return_value < 0)
        {
            // low-level reading failed so report the return code.
            return return_value;
        }
        const auto is_supplementary
            = static_cast<const bool>(_impl->hts_bam_align_ptr_->core.flag & kSupplementaryAlign);
        const auto is_secondary = static_cast<const bool>(_impl->hts_bam_align_ptr_->core.flag & kSecondaryAlign);
        is_primary_align = (!is_supplementary) && (!is_secondary);
    }
    return return_value;
}

// Try to get an aligned mate.
bool BamReader::getAlignedMate(const Read& read, Read& mate)
{
    int32_t tid = 0;
    int32_t beg = 0;
    int32_t end = 0;

    if (read.is_mate_mapped())
    {
        tid = read.mate_chrom_id();
        beg = read.mate_pos();
        end = read.mate_pos() + 1;
    }
    else
    {
        tid = read.chrom_id();
        beg = read.pos();
        end = read.pos() + 1;
    }

    hts_itr_t* iter;
    iter = sam_itr_queryi(_impl->hts_idx_ptr_, tid, beg, end);
    if (iter == nullptr)
    {
        return false;
    }
    while (sam_itr_next(_impl->hts_file_ptr_, iter, _impl->hts_bam_align_ptr_) >= 0)
    {
        decodeHtsAlign(_impl->hts_bam_align_ptr_, mate);
        if ((mate.fragment_id() == read.fragment_id()) && (mate.is_first_mate() != read.is_first_mate()))
        {
            hts_itr_destroy(iter);
            return true;
        }
    }
    hts_itr_destroy(iter);
    return false;
}

std::unique_ptr<DepthInfo> BamReader::estimateDepth(std::string const& region)
{
    auto logger = LOG();

    // max change of DP for convergence
    static const double kDPAccuracy = 0.05;
    // Minimum region length to estimate depth over
    static const auto kIntervalLength = 2000000;

    // output variable that will be returned
    std::unique_ptr<DepthInfo> dp_info{ new DepthInfo };

    // file needs to be open already so we have a path
    // we will open a separate handle below s.t. we can
    // run this function in parallel with itself + read CRAM files
    // which have side effects and cannot have multiple iterators
    // opened for the same file / index
    assert(!_impl->file_path.empty());
    _impl->open(_impl->file_path, _impl->reference_path);

    std::string chr;
    int64_t region_start = -1;
    int64_t region_end = -1;
    common::stringutil::parsePos(region, chr, region_start, region_end);

    assert(_impl->header_contig_map.count(chr) != 0);

    const int tid = _impl->header_contig_map[chr];
    {
        if (region_start < 0)
        {
            region_start = 0;
        }
        if (region_end < 0)
        {
            region_end = _impl->hts_bam_hdr_ptr_->target_len[tid];
        }
    }

    if (region_end < region_start)
    {
        region_end = region_start;
    }

    // walk through chromosome in intervals; remove intervals with no reads
    std::vector<std::pair<int64_t, int64_t>> estimation_intervals{ { region_start, region_end } };
    int64_t max_interval_size = region_end - region_start + 1;
    while (max_interval_size > kIntervalLength && estimation_intervals.size() < 20)
    {
        std::vector<std::pair<int64_t, int64_t>> new_intervals;
        max_interval_size = 0;
        for (const auto& iv : estimation_intervals)
        {
            const auto iv_length = iv.second - iv.first + 1;
            if (iv_length > kIntervalLength)
            {
                new_intervals.emplace_back(iv.first, iv.first + iv_length / 2);
                if (iv_length > 1)
                {
                    new_intervals.emplace_back(iv.first + iv_length / 2 + 1, iv.second);
                }
                max_interval_size = std::max(max_interval_size, (iv_length + 1) / 2);
            }
            else
            {
                new_intervals.push_back(iv);
                max_interval_size = std::max(max_interval_size, iv_length);
            }
        }
        estimation_intervals = new_intervals;
    }

    std::vector<bool> intervals_empty;
    intervals_empty.resize(estimation_intervals.size(), false);

    std::unique_ptr<bam1_t, std::function<void(bam1_t*)>> hts_align_ptr{ bam_init1(), bam_destroy1 };

    using namespace boost::accumulators;
    accumulator_set<double, features<tag::count, tag::median, tag::variance>> read_length_accumulator;
    accumulator_set<double, features<tag::median, tag::variance>> depth_accumulator;

    bool all_finished = false;
    bool converged = false;
    double convergence_previous_depth = std::numeric_limits<double>::max();

    int cycle = 0;

    auto check_convergence = [&]() {
        if (cycle == 0)
        {
            // look at each segment
            return;
        }
        // check every ~1M reads
        const double current_depth = median(depth_accumulator);
        if (fabs(current_depth - convergence_previous_depth) < kDPAccuracy)
        {
#ifdef _DEBUG
            logger->debug("Converged for {} DP: {}", chr, current_depth);
#endif
            converged = true;
        }
        else
        {
#ifdef _DEBUG
            logger->debug(
                "Not converged yet for {} previous: {} current: {}", chr, convergence_previous_depth, current_depth);
#endif
            convergence_previous_depth = current_depth;
        }
    };

    BamReaderImpl tmp_impl;
    tmp_impl.file_path = _impl->file_path;
    tmp_impl.reference_path = _impl->reference_path;
    tmp_impl.open(tmp_impl.file_path, tmp_impl.reference_path);
    while (!(converged || all_finished) && cycle < 10)
    {
        for (size_t interval_ptr = 0; interval_ptr < estimation_intervals.size(); ++interval_ptr)
        {
            if (intervals_empty[interval_ptr])
            {
                continue;
            }
            const int64_t start = estimation_intervals[interval_ptr].first;
            const int64_t end = estimation_intervals[interval_ptr].second;
#ifdef _DEBUG
            logger->debug("Starting interval {}:{}-{}", chr, start, end);
#endif
            std::unique_ptr<hts_itr_t, std::function<void(hts_itr_t*)>> hts_itr_ptr{
                sam_itr_queryi(tmp_impl.hts_idx_ptr_, tid, static_cast<int>(start), static_cast<int>(end)),
                hts_itr_destroy
            };
            if (!hts_itr_ptr)
            {
                intervals_empty[interval_ptr] = true;
                continue;
            }

            int64_t last_pos = start;

            int return_value;
            int any_reads = 0;

            ReadPileup pileup;

            while ((return_value = sam_itr_next(tmp_impl.hts_file_ptr_, hts_itr_ptr.get(), hts_align_ptr.get())) >= 0)
            {
                const auto is_supplementary = static_cast<const bool>(hts_align_ptr->core.flag & kSupplementaryAlign);
                const auto is_secondary = static_cast<const bool>(hts_align_ptr->core.flag & kSecondaryAlign);
                const auto is_primary_align = (!is_supplementary) && (!is_secondary);

                if (hts_align_ptr->core.tid != tid || !is_primary_align || bam_get_qual(hts_align_ptr.get()) == 0
                    || hts_align_ptr->core.pos + hts_align_ptr->core.l_qseq < start)
                {
                    continue;
                }

                // update read length
                any_reads++;
                read_length_accumulator(hts_align_ptr->core.l_qseq);

                // update pileup
                pileup.addRead(SimpleRead::make(hts_align_ptr->core.pos, hts_align_ptr->core.l_qseq));

                last_pos = static_cast<int64_t>(hts_align_ptr->core.pos);

                // change intervals every 10kb / 40k reads
                if (last_pos - start > 10000 && any_reads > 40000)
                {
                    break;
                }
            }

            int64_t last_pileup_pos = start;
            const size_t current_read_length = median(read_length_accumulator);
            while (last_pileup_pos <= last_pos)
            {
                size_t dp = 0;
                pileup.pileup(last_pileup_pos, [&dp](ReadPileupInfo*) { ++dp; });
                depth_accumulator(dp);
                last_pileup_pos += std::max((size_t)1, current_read_length / 2);
            }

            if (return_value < 0)
            {
                intervals_empty[interval_ptr] = true;
            }
            else
            {
                estimation_intervals[interval_ptr].first = last_pos;
                estimation_intervals[interval_ptr].second
                    = std::max(last_pos, estimation_intervals[interval_ptr].second);
            }
            if (any_reads > 10000 && cycle > 0)
            {
                check_convergence();
                if (converged)
                {
                    break;
                }
            }
        }
        all_finished = std::all_of(intervals_empty.begin(), intervals_empty.end(), [](bool b) { return b; });
        ++cycle;
        check_convergence();
    }

    dp_info->read_length_unique = fabs(variance(read_length_accumulator)) < std::numeric_limits<double>::epsilon();
    dp_info->read_length = static_cast<size_t>(median(read_length_accumulator));
    dp_info->read_count = static_cast<size_t>(count(read_length_accumulator));
    dp_info->depth_median = round(median(depth_accumulator) * 100) / 100;
    dp_info->depth_variance = round(variance(depth_accumulator) * 100 / 100);

    if (!converged && !all_finished)
    {
        logger->warn("DP estimation for region {} did not converge.", region);
    }

    return dp_info;
}
}