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

#include "genotyping/VariantGenotype.hh"
#include "common/Error.hh"

using std::string;
using std::vector;

namespace genotyping
{

void VariantGenotype::addInfoFromSingleBreakpoint(const Genotype& breakpoint_genotype)
{
    if (breakpoint_genotype.gt.empty())
    {
        num_missing_++;
        return;
    }
    for (auto& possible_gt : breakpoint_genotype.gl_name)
    {
        if (gt_vote_info_.find(possible_gt) == gt_vote_info_.end())
        {
            gt_vote_info_[possible_gt] = genotypeVoteInfo();
        }
        if (possible_gt == breakpoint_genotype.gt)
        {
            gt_vote_info_[possible_gt].num_support++;
        }
        else
        {
            gt_vote_info_[possible_gt].num_unsupport++;
        }
    }
    for (auto& equal_gt : breakpoint_genotype.equivalent_gt)
    {
        if (gt_vote_info_.find(equal_gt) == gt_vote_info_.end())
        {
            gt_vote_info_[equal_gt] = genotypeVoteInfo();
        }
        gt_vote_info_[equal_gt].num_support++;
    }
}

void VariantGenotype::genotype(GenotypeVector reference_homozygote)
{
    vector<GenotypeVector> best_genotype_vec;
    int max_num_support = -1;
    for (auto& current_gt : gt_vote_info_)
    {
        if (current_gt.second.num_support > max_num_support)
        {
            best_genotype_vec.clear();
            best_genotype_vec.push_back(current_gt.first);
            max_num_support = current_gt.second.num_support;
        }
        else if (current_gt.second.num_support == max_num_support)
        {
            best_genotype_vec.push_back(current_gt.first);
        }
    }
    if (max_num_support <= 0)
    {
        gt_ = {};
        filter_ = "ZERO_COUNT";
    }
    else
    {
        if (best_genotype_vec.size() > 1)
        {
            genotypeFromDecisionTree(reference_homozygote, best_genotype_vec);
            filter_ = "CONFLICTS";
        }
        else
        {
            gt_ = best_genotype_vec[0];
            filter_ = "PASS";
        }
    }
}

GenotypeVector VariantGenotype::genotypeFromDecisionTree(
    GenotypeVector& reference_homozygote, vector<GenotypeVector>& best_genotype_vec)
{
    GenotypeVector best_genotype;
    if (best_genotype_vec.size() == 2)
    {
        if (best_genotype_vec[0] == reference_homozygote)
        {
            best_genotype = best_genotype_vec[1];
        }
        else if (best_genotype_vec[1] == reference_homozygote)
        {
            best_genotype = best_genotype_vec[0];
        }
    }
    else
    {
        best_genotype = {};
    }
    return best_genotype;
}

string VariantGenotype::genotypeString() const
{
    string gt_str;
    if (gt_.empty())
    {
        gt_str = "./.";
    }
    else
    {
        gt_str = std::to_string(gt_[0]) + "/" + std::to_string(gt_[1]);
    }
    return gt_str;
}

string VariantGenotype::genotypeString(vector<string>& allele_names) const
{
    string gt_str;
    if (gt_.empty())
    {
        gt_str = "./.";
    }
    else
    {
        gt_str = allele_names[gt_[0]] + "/" + allele_names[gt_[1]];
    }
    return gt_str;
}

VariantGenotype::operator std::string() const { return genotypeString(); }

// output allele name instead of indexes
Json::Value VariantGenotype::toJson(vector<string> allele_names) const
{
    Json::Value json_result;

    json_result["GT"] = genotypeString(allele_names);

    if (!equivalent_gt_.empty())
    {
        json_result["equivalentGT"] = Json::Value();
        for (int i = 0; i < (int)equivalent_gt_.size(); i++)
        {
            string equal_gt_str = allele_names[equivalent_gt_[i][0]] + "/" + allele_names[equivalent_gt_[i][1]];
            json_result["equivalentGT"][i] = equal_gt_str;
        }
    }
    if (!filter_.empty())
    {
        json_result["Filter"] = filter_;
    }

    if (!gt_.empty())
    {
        const genotypeVoteInfo best_vote = gt_vote_info_.at(gt_);
        if (best_vote.num_support > 0)
        {
            json_result["num_support"] = best_vote.num_support;
        }
        if (best_vote.num_unsupport > 0)
        {
            json_result["num_unsupport"] = best_vote.num_unsupport;
        }
    }
    if (num_missing_ > 0)
    {
        json_result["num_missing"] = num_missing_;
    }

    if (!gt_vote_info_.empty())
    {
        for (auto& vote_info : gt_vote_info_)
        {
            if (vote_info.first == gt_)
            {
                continue;
            }
            string current_gt = allele_names[vote_info.first[0]] + "/" + allele_names[vote_info.first[1]];
            json_result["otherGT"][current_gt] = Json::Value();
            if (vote_info.second.num_support > 0)
            {
                json_result["otherGT"][current_gt]["num_support"] = vote_info.second.num_support;
            }
            if (vote_info.second.num_unsupport > 0)
            {
                json_result["otherGT"][current_gt]["num_unsupport"] = vote_info.second.num_unsupport;
            }
        }
    }

    return json_result;
}
}
