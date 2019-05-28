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

namespace idxdepth
{

struct IndexBin
{
    int bin_id = 0;
    std::string chrom;
    int64_t start = -1;
    int64_t end = -1;
    size_t slices = 0; // number of times this index bin has been encountered (can be >1 for CRAM indices)
    size_t overlapping_bytes = 0; // total byte count from overlapping index slices
    double adjusted_bytes = 0; // adjusted byte count (after adjusting for slice vs bin size)
    double normalized_depth = 0;
};

/**
 * Get index bins
 * @param bam BAM / CRAM file
 * @param output output vector for bins
 */
void getIndexBins(std::string const& bam, std::vector<IndexBin>& output);
}