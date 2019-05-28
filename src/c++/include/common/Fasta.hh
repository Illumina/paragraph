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
 *  \brief Wrapper for htslib faidx
 *
 * \file Fasta.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Reference.hh"
#include <list>
#include <string>

namespace common
{

struct FastaFileImpl;

class FastaFile : public Reference
{
public:
    FastaFile();

    explicit FastaFile(std::string const& filename);

    ~FastaFile();

    FastaFile(FastaFile const&);

    FastaFile& operator=(FastaFile const&);

    std::string getFilename() const;

    std::string query(std::string const& location) const;

    std::string query(const char* chr, int64_t start, int64_t end) const;

    /**
     * return the size of a contig. Sum of all contig sizes if contig == ""
     * @param contig name of contig, or empty
     * @return number of bases in the contig
     */
    size_t contigSize(std::string const& contig = "") const;

    /**
     * return the non-N padded size of a contig. This is calculated
     * as the size of the contig minus any N's at the beginning or
     * at the end.
     *
     * Returns the sum of all contig sizes if contig == ""
     * @param contig name of contig, or empty
     * @return number of bases in the contig
     */
    size_t contigNonNSize(std::string const& contig = "") const;

    /**
     * Return the first position that is not an N character
     */
    size_t contigNonNStart(std::string const& contig) const;

    /**
     * @return all contig names
     */
    std::list<std::string> getContigNames() const;

private:
    FastaFileImpl* _impl;
};
}
