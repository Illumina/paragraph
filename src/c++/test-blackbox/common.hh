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