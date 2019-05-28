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

string Genotype::toString(const vector<string>* p_allele_names) const
{
    string gt_str;
    if (gt.empty())
    {
        gt_str = ".";
    }
    else
    {
        if (p_allele_names == NULL)
        {
            auto al = static_cast<std::string (*)(int)>(std::to_string);
            gt_str = join(gt | transformed(al), "/");
        }
        else
        {
            auto al = [p_allele_names](int g) { return (*p_allele_names)[g]; };
            gt_str = join(gt | transformed(al), "/");
        }
    }
    return gt_str;
}

std::string Genotype::filterString() const
{
    string filter_string = boost::algorithm::join(filters, ";");
    return filter_string;
}

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
    json_result["GT"] = toString(&allele_names);

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

    if (gq != -1)
    {
        json_result["GQ"] = gq;
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
    if (!filters.empty())
    {
        json_result["filters"] = Json::arrayValue;
        for (auto f : filters)
        {
            json_result["filters"].append(f);
        }
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

Genotype::operator std::string() const { return toString(); }
}
