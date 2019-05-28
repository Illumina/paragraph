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
     * @param index_path bam/cram index file path (pass an empty string to use default location)
     * @param reference path to FASTA reference
     */
    explicit BamReader(const std::string& path, const std::string& index_path, const std::string& reference);

    BamReader(BamReader&&) noexcept;
    BamReader& operator=(BamReader&&) noexcept;
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;

    ~BamReader() override;

    /**
     * Initialize BAM file for reading at a given region.
     * return true if find the region
     */
    void setRegion(const std::string& region_encoding = "") override;

    /**
     *  Reads the next alignment from the region
     * return false if the region is exhausted.
     */
    bool getAlign(Read& align) override;

    /**
     *  fetch information for the mate of this aligned read
     * return true if found
     */
    bool getAlignedMate(const Read& read, Read& mate) override;

    /**
     * estimate depth on a given region
     */
    std::unique_ptr<DepthInfo> estimateDepth(std::string const& region) override;

protected:
    int SkipToNextGoodAlign();

private:
    struct BamReaderImpl;
    std::unique_ptr<BamReaderImpl> _impl;
};
}
