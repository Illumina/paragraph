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

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <set>
#include <string>
#include <unordered_map>

#include "common/Error.hh"
#include "genotyping/CombinedGenotype.hh"

using std::string;
using std::vector;

namespace genotyping
{

Genotype combinedGenotype(GenotypeSet const& genotypes)
{
    Genotype result;

    std::set<std::string> filters;

    // combine GLs
    struct GLInfo
    {
        explicit GLInfo(GenotypeVector gt_, double gl_)
            : gt(gt_)
            , gl(gl_)
        {
        }
        explicit GLInfo(GLInfo const& rhs) = default;
        GenotypeVector gt;
        double gl = std::numeric_limits<double>::min();
    };
    std::unordered_map<std::string, GLInfo> GLs;
    result.num_reads = 0;
    result.allele_fractions.clear();

    // count votes for each GT
    std::map<GenotypeVector, std::list<Genotype>> votes_for;
    for (auto& bp : genotypes)
    {
        if (bp.gt.empty())
        {
            filters.insert("MISSING");
            continue;
        }
        // sort genotype -- this will need to handle phasing
        // at some point
        Genotype sorted_bp = bp;
        std::sort(sorted_bp.gt.begin(), sorted_bp.gt.end());
        auto v_it = votes_for.find(sorted_bp.gt);
        if (v_it == votes_for.end())
        {
            votes_for[sorted_bp.gt] = { bp };
        }
        else
        {
            v_it->second.push_back(bp);
        }
        result.num_reads += bp.num_reads;
        if (bp.allele_fractions.size() > result.allele_fractions.size())
        {
            result.allele_fractions.resize(bp.allele_fractions.size(), 0);
        }
        for (size_t i = 0; i < bp.allele_fractions.size(); ++i)
        {
            result.allele_fractions[i] += bp.num_reads * bp.allele_fractions[i];
        }

        assert(bp.gl.size() == bp.gl_name.size());
        for (size_t i = 0; i < bp.gl.size(); ++i)
        {
            using boost::adaptors::transformed;
            using boost::algorithm::join;
            auto str = static_cast<std::string (*)(int)>(std::to_string);
            auto sorted_gl_name = bp.gl_name[i];
            std::sort(sorted_gl_name.begin(), sorted_gl_name.end());
            const std::string gl_gt_as_string = join(sorted_gl_name | transformed(str), "|");
            auto gl_it = GLs.find(gl_gt_as_string);
            if (gl_it == GLs.end())
            {
                const std::unordered_map<std::string, GLInfo>::value_type entry{ gl_gt_as_string,
                                                                                 GLInfo(sorted_gl_name, bp.gl[i]) };
                GLs.insert(entry);
            }
            else
            {
                gl_it->second.gl = std::max(gl_it->second.gl, bp.gl[i]);
            }
        }
    }

    for (auto& af : result.allele_fractions)
    {
        af /= result.num_reads;
    }

    result.gl.reserve(GLs.size());
    result.gl_name.reserve(GLs.size());
    for (auto const& gl : GLs)
    {
        result.gl.push_back(gl.second.gl);
        result.gl_name.push_back(gl.second.gt);
    }

    if (votes_for.size() > 1)
    {
        filters.insert("CONFLICT");
        // output most popular GT that is not reference
        size_t best_bp_votes = 0;
        size_t second_best_bp_votes = 0;
        GenotypeVector best_bp;
        size_t ref_ploidy = 0;
        for (auto& bp_votes : votes_for)
        {
            const bool this_bp_is_missing = bp_votes.first.empty();
            if (this_bp_is_missing)
            {
                continue;
            }

            const bool this_bp_is_ref
                = std::all_of(bp_votes.first.begin(), bp_votes.first.end(), [](int32_t g) { return g == 0; });
            if (this_bp_is_ref)
            {
                ref_ploidy = bp_votes.first.size();
                continue;
            }

            const size_t this_bp_votes = bp_votes.second.size();
            // compare vote counts for all cases where we called an allele
            if (this_bp_votes > best_bp_votes)
            {
                second_best_bp_votes = std::max(second_best_bp_votes, best_bp_votes);
                best_bp_votes = this_bp_votes;
                best_bp = bp_votes.first;
            }
            else
            {
                second_best_bp_votes = std::max(second_best_bp_votes, best_bp_votes);
            }
        }
        // -> implies also best_bp_votes > 0
        if (best_bp_votes > second_best_bp_votes)
        {
            result.gt = best_bp;
        }
        else if (best_bp_votes == 0 && ref_ploidy > 0)
        {
            result.gt.empty();
            result.gt.resize(ref_ploidy, 0);
        }
        else
        {
            result.gt.empty();
        }
    }
    else if (votes_for.size() == 1)
    {
        // everyone voted for the same GT
        result.gt = votes_for.begin()->first;
    }
    else
    {
        filters.insert("NOGT");
    }

    if (filters.empty())
    {
        result.filter = "PASS";
    }
    else
    {
        result.filter = boost::algorithm::join(filters, ";");
    }

    return result;
}
}
