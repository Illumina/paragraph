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
 * \brief Genotyping helpers
 *
 * \file Genotype.cpp
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/Genotype.hh"
#include "common/Error.hh"

using std::string;
using std::vector;

namespace genotyping
{
string Genotype::genotypeString() const
{
    string gt_str;
    if (gt.empty())
    {
        gt_str = "./.";
    }
    else
    {
        gt_str = std::to_string(gt[0]) + "/" + std::to_string(gt[1]);
    }
    return gt_str;
}

string Genotype::genotypeString(vector<string>& allele_names) const
{
    string gt_str;
    if (gt.empty())
    {
        gt_str = "./.";
    }
    else
    {
        gt_str = allele_names[gt[0]] + "/" + allele_names[gt[1]];
    }
    return gt_str;
}

Genotype::operator std::string() const { return genotypeString(); }

void Genotype::recode(vector<vector<uint64_t>>& alternative_codes)
{
    vector<uint64_t> hap1 = alternative_codes[gt[0]];
    vector<uint64_t> hap2 = alternative_codes[gt[1]];

    if (hap1[0] <= hap2[0])
    {
        gt = { hap1[0], hap2[0] };
    }
    else
    {
        gt = { hap2[0], hap1[0] };
    }

    // figure out other possibilities
    for (auto& equal_gt : equivalent_gt)
    {
        hap1.insert(hap1.end(), alternative_codes[equal_gt[0]].begin(), alternative_codes[equal_gt[0]].end());
        hap2.insert(hap2.end(), alternative_codes[equal_gt[1]].begin(), alternative_codes[equal_gt[1]].end());
    }
    std::map<GenotypeVector, bool> possible_gt_hash;
    for (auto& h1 : hap1)
    {
        for (auto& h2 : hap2)
        {
            GenotypeVector current_gt;
            if (h1 <= h2)
            {
                current_gt = { h1, h2 };
            }
            else
            {
                current_gt = { h2, h1 };
            }
            if (possible_gt_hash.find(current_gt) == possible_gt_hash.end())
            {
                possible_gt_hash[current_gt] = true;
            }
        }
    }

    // implement equivalent gt
    equivalent_gt.clear();
    for (auto& possible_gt : possible_gt_hash)
    {
        if (possible_gt.first != gt)
        {
            equivalent_gt.push_back(possible_gt.first);
        }
    }

    // implement gl
    for (auto& element : gl_name)
    {
        uint64_t gl0 = alternative_codes[element[0]][0];
        uint64_t gl1 = alternative_codes[element[1]][0];
        if (gl0 <= gl1)
        {
            element = { gl0, gl1 };
        }
        else
        {
            element = { gl1, gl0 };
        }
    }
}

// output allele name instead of indexes
Json::Value Genotype::toJson(vector<string>& allele_names) const
{
    Json::Value json_result;
    json_result["GT"] = genotypeString(allele_names);

    if (!equivalent_gt.empty())
    {
        json_result["equivalentGT"] = Json::arrayValue;
        for (auto& equal_gt : equivalent_gt)
        {
            string equal_gt_str = allele_names[equal_gt[0]] + "/" + allele_names[equal_gt[1]];
            json_result["equivalentGT"].append(equal_gt_str);
        }
    }

    if (!gl.empty())
    {
        json_result["GL"] = Json::Value();
        for (size_t i = 0; i < gl.size(); ++i)
        {
            string gl_str = allele_names[gl_name[i][0]] + "/" + allele_names[gl_name[i][1]];
            json_result["GL"][gl_str] = gl[i];
        }
    }

    if (!allele_fractions.empty())
    {
        json_result["allele_fractions"] = Json::arrayValue;
        for (auto& allele_fraction : allele_fractions)
        {
            json_result["allele_fractions"].append(allele_fraction);
        }
    }
    if (!filter.empty())
    {
        json_result["Filter"] = filter;
    }
    if (!gt.empty())
    {
        json_result["DP"] = num_reads;
    }
    return json_result;
}
}
