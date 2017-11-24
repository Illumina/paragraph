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

#include <string>
#include <vector>

#include "common/DepthEstimator.hh"
#include "common/Read.hh"
#include "common/ReadReader.hh"

namespace common
{

/**
 * @brief BAM reader class
 *
 * Extracts alignments from BAM/CRAM file contained in a specified region.
 * Example:
 *   reader.init("crams/LP6008113-DNA_A01.cram", "chr2:44872230-44875849");
 *   BamAlignment al;
 *   while (reader.getAlign(al)) {
 *     cerr << al.bases << endl;
 *   }
 */
class BamReader : public ReadReader, public DepthEstimator
{
public:
    enum
    {
        kSupplementaryAlign = 0x800,
        kSecondaryAlign = 0x100,
        kIsMapped = 0x0004,
        kIsFirstMate = 0x0040,
        kIsMateMapped = 0x0008
    };

    /**
     * Create BAM / CRAM reader. Reference is required to match the
     * reference FASTA file for the BAM/CRAM
     * @param path file path
     * @param reference path to FASTA reference
     */
    explicit BamReader(const std::string& path, const std::string& reference);

    BamReader(BamReader&&) noexcept;
    BamReader& operator=(BamReader&&) noexcept;
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;

    ~BamReader() override;

    // Initialize BAM file for reading at a given region.
    void setRegion(const std::string& region_encoding = "");
    // Reads the next alignment from the region; returns false if the region
    // is exhausted.
    bool getAlign(Read& align) override;

    bool getAlignedMate(const Read& read, Read& mate) override;

    std::unique_ptr<DepthInfo> estimateDepth(std::string const& region) override;

protected:
    int SkipToNextGoodAlign();

private:
    struct BamReaderImpl;
    std::unique_ptr<BamReaderImpl> _impl;
};
}
