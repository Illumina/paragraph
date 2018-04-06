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

Genotype
combinedGenotype(GenotypeSet const& genotypes, const BreakpointGenotyper* p_genotyper, double depth, int read_length)
{
    Genotype result;

    size_t num_pass_genotypes = countUniqGenotypes(genotypes, true);
    if (num_pass_genotypes == 0)
    {
        int num_fail_genotypes = countUniqGenotypes(genotypes, false);
        if (num_fail_genotypes == 0)
        {
            result.filter = "MISSING";
        }
        else if (num_fail_genotypes == 1)
        {
            result = reportConsensusGenotypes(genotypes, false);
        }
        else
        {
            result = genotypeByTotalCounts(genotypes, false, p_genotyper, depth, read_length);
        }
    }
    else if (num_pass_genotypes == 1)
    {
        result = reportConsensusGenotypes(genotypes, true);
    }
    else
    {
        result = genotypeByTotalCounts(genotypes, true, p_genotyper, depth, read_length);
    }

    if (result.filter.empty())
    {
        result.filter = "PASS";
    }

    return result;
}

size_t countUniqGenotypes(GenotypeSet const& genotypes, bool pass_only)
{
    std::set<GenotypeVector> votes_for;
    for (auto& bp : genotypes)
    {
        if (bp.gt.empty())
        {
            continue;
        }

        if (pass_only)
        {
            if (!bp.filter.empty())
            {
                continue;
            }
        }

        Genotype sorted_bp = bp;
        std::sort(sorted_bp.gt.begin(), sorted_bp.gt.end());
        votes_for.insert(sorted_bp.gt);
    }

    return votes_for.size();
}

Genotype reportConsensusGenotypes(GenotypeSet const& genotypes, bool pass_only)
{
    Genotype result;
    std::set<string> filters;

    if (!pass_only)
    {
        filters.insert("ALL_BAD_BP");
    }

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

    for (auto& bp : genotypes)
    {

        if (bp.gt.empty())
        {
            filters.insert("MISSING");
            continue;
        }

        if (pass_only)
        {
            if (!bp.filter.empty())
            {
                filters.insert("EXIST_BAD_BP");
                continue;
            }
        }

        if (result.gt.empty())
        {
            // sort genotype -- this will need to handle phasing at some point
            Genotype sorted_bp = bp;
            std::sort(sorted_bp.gt.begin(), sorted_bp.gt.end());
            result.gt = sorted_bp.gt;
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

    result.filter = boost::algorithm::join(filters, ";");

    return result;
}

Genotype genotypeByTotalCounts(
    GenotypeSet const& genotypes, bool use_pass_only, const BreakpointGenotyper* p_genotyper, double depth,
    int read_length)
{
    assert(p_genotyper != NULL && depth > 0 && read_length > 0);

    std::set<string> filters;
    filters.insert("CONFLICT");
    if (!use_pass_only)
    {
        filters.insert("ALL_BAD_BP");
    }

    // get total count vector
    std::vector<int> sum_counts;
    int num_bp = 0;
    for (auto const& bp : genotypes)
    {
        if (use_pass_only)
        {
            if (!bp.filter.empty())
            {
                filters.insert("EXIST_BAD_BP");
                continue;
            }
        }

        if (bp.num_reads == 0)
        {
            filters.insert("MISSING");
            continue;
        }

        if (sum_counts.empty())
        {
            sum_counts.resize(bp.allele_fractions.size(), 0);
        }
        else
        {
            assert(sum_counts.size() == bp.allele_fractions.size());
        }
        size_t allele_index = 0;
        for (auto& af : bp.allele_fractions)
        {
            sum_counts[allele_index] += std::round(af * bp.num_reads);
            allele_index++;
        }
        num_bp++;
    }

    // genotype with mean
    for (auto& s : sum_counts)
    {
        s = std::round((double)s / num_bp);
    }
    const auto input_counts = sum_counts;
    const auto sum_gt = p_genotyper->genotype(depth, read_length, input_counts);
    Genotype result = sum_gt;
    result.filter = boost::algorithm::join(filters, ";");
    return result;
}
}
