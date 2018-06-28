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
#include "common/StringUtil.hh"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <algorithm>
#include <cassert>
#include <functional>

using boost::adaptors::transformed;
using boost::algorithm::join;
using std::string;
using std::vector;

namespace genotyping
{

GenotypeVector genotypeVectorFromString(std::string& gt_str)
{
    GenotypeVector gv;
    size_t prev_pos = 0;
    size_t current_pos = gt_str.find("/");
    while (current_pos != std::string::npos)
    {
        string allele_str = gt_str.substr(prev_pos, current_pos - prev_pos);
        assert(std::all_of(allele_str.begin(), allele_str.end(), ::isdigit));
        gv.push_back(stoi(allele_str));
        prev_pos = current_pos + 1;
        current_pos = gt_str.find("/", prev_pos);
    }
    return gv;
}

string Genotype::toString() const
{
    string gt_str;
    if (gt.empty())
    {
        gt_str = "."; // TODO: change this for haploid
    }
    else
    {
        auto str = static_cast<std::string (*)(int)>(std::to_string);
        gt_str = join(gt | transformed(str), "/");
    }
    return gt_str;
}

string Genotype::toString(vector<string> const& allele_names) const
{
    string gt_str;
    if (gt.empty())
    {
        gt_str = "./.";
    }
    else
    {
        auto al = [&allele_names](int g) { return allele_names[g]; };
        gt_str = join(gt | transformed(al), "/");
    }
    return gt_str;
}

Genotype::operator std::string() const { return toString(); }

/**
 * relabel genotypes
 */
void Genotype::relabel(std::vector<uint64_t> const& new_labels)
{
    // remap genotypes
    for (auto& g : gt)
    {
        assert(new_labels.size() > g);
        g = new_labels[g];
    }
    std::sort(gt.begin(), gt.end());
    // remap labels for GL
    for (auto& l : gl_name)
    {
        for (auto& g : l)
        {

            assert(new_labels.size() > g);
            g = new_labels[g];
        }
        std::sort(l.begin(), l.end());
    }

    // relabel allele fractions
    std::vector<double> new_allele_fractions;
    new_allele_fractions.resize(new_labels.size(), 0ull);
    for (size_t g = 0; g < allele_fractions.size(); ++g)
    {
        assert(new_labels.size() > g);
        const uint64_t t = new_labels[g];
        new_allele_fractions[t] = allele_fractions[g];
    }

    allele_fractions = new_allele_fractions;
}

// output allele name instead of indexes
Json::Value Genotype::toJson(vector<string> const& allele_names) const
{
    Json::Value json_result;
    json_result["GT"] = toString(allele_names);

    if (!gl.empty())
    {
        json_result["GL"] = Json::Value();
        for (size_t i = 0; i < gl.size(); ++i)
        {
            string gl_str = allele_names[gl_name[i][0]];
            for (size_t j = 1; j < gl_name[i].size(); j++)
            {
                gl_str += "/" + allele_names[gl_name[i][j]];
            }
            json_result["GL"][gl_str] = gl[i];
        }
    }

    if (!allele_fractions.empty())
    {
        json_result["allele_fractions"] = Json::objectValue;

        assert(allele_names.size() == allele_fractions.size());
        for (size_t allele = 0; allele < allele_fractions.size(); ++allele)
        {
            const std::string& allele_name = allele_names[allele];
            const double& allele_fraction = allele_fractions[allele];
            json_result["allele_fractions"][allele_name] = allele_fraction;
        }
    }
    if (!filter.empty())
    {
        json_result["filter"] = filter;
    }
    if (!gt.empty())
    {
        json_result["num_reads"] = num_reads;
        if (coverage_test_pvalue != -1)
        {
            json_result["coverage_test_pvalue"] = coverage_test_pvalue;
        }
    }
    return json_result;
}
}
