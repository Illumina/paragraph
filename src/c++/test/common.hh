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

/**
 *
 * \file common.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include "gtest/gtest.h"

#include <cstdlib>
#include <list>
#include <string>
#include <unistd.h>

class GTestEnvironment : public testing::Environment
{
public:
    explicit GTestEnvironment(const char* _bpath)
        : bpath(_bpath)
    {

        std::list<std::string> paths_to_try_hg19
            = { "/illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" };

        const char* hg19_var = std::getenv("HG19");
        if (hg19_var != nullptr)
        {
            paths_to_try_hg19.push_front(std::string(hg19_var));
        }
        for (const auto& p : paths_to_try_hg19)
        {
            if (access(p.c_str(), F_OK) != -1)
            {
                hg19path = p;
                break;
            }
        }

        std::list<std::string> paths_to_try_hg38
            = { "/illumina/sync/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa" };

        const char* hg38_var = std::getenv("HG38");
        if (hg38_var != nullptr)
        {
            paths_to_try_hg38.push_front(std::string(hg38_var));
        }
        for (const auto& p : paths_to_try_hg38)
        {
            if (access(p.c_str(), F_OK) != -1)
            {
                hg38path = p;
                break;
            }
        }
    }
    ~GTestEnvironment() = default;

    std::string getBasePath() const { return bpath; }

    std::string getHG19Path() const { return hg19path; }

    std::string getHG38Path() const { return hg38path; }

protected:
    std::string bpath;
    std::string hg19path;
    std::string hg38path;
};

extern GTestEnvironment* g_testenv;
