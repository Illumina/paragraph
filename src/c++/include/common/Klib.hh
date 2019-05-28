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
 * \brief klib Smith Waterman wrapped. See external/klib.
 *
 * \file Klib.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "common/Alignment.hh"

namespace common
{
struct KlibAlignmentImpl;

class KlibAlignment : public Alignment
{
public:
    KlibAlignment();
    KlibAlignment(KlibAlignment&& that);
    KlibAlignment(const KlibAlignment& that) = delete;

    ~KlibAlignment();

    void setParameters(AlignmentParameters const& ap);

    void getParameters(AlignmentParameters& ap);

    /**
     * @brief set target sequence
     */
    void setRef(const char* seq);

    /*
     * @brief set query sequence
     */
    void setQuery(const char* seq);

    /**
     * @brief Get the alignment score
     */
    int getScore();

    /**
     * @brief Get a cigar string
     */
    void getCigar(int& r0, int& r1, int& a0, int& a1, std::string& cig);

    void getCigar(int& r0, int& r1, int& a0, int& a1, int& n_cigar, uint32_t*& cigar);

protected:
    virtual void update();

    KlibAlignmentImpl* _impl;
};
}
