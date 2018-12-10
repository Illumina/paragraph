//
// Copyright (c) 2016 Illumina, Inc.
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

#include "genotyping/GraphBreakpointGenotyper.hh"
#include "genotyping/BreakpointGenotyper.hh"
#include "genotyping/CombinedGenotype.hh"

#include "common/Error.hh"
#include "common/JsonHelpers.hh"

using std::map;
using std::string;
using std::vector;

namespace genotyping
{

void GraphBreakpointGenotyper::setParameters(const string& genotyping_parameter_path)
{
    p_genotype_parameter
        = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleleNames(), female_ploidy_));
    p_male_genotype_parameter
        = std::unique_ptr<GenotypingParameters>(new GenotypingParameters(alleleNames(), male_ploidy_));
    if (!genotyping_parameter_path.empty())
    {
        Json::Value param_json = common::getJSON(genotyping_parameter_path);
        p_genotype_parameter->setFromJson(param_json);
    }
}

void GraphBreakpointGenotyper::runGenotyping()
{
    auto breakpoint_names = breakpointNames();

    // genotype all breakpoints
    auto const& allelenames = alleleNames();
    BreakpointGenotyper genotyper(p_genotype_parameter);
    BreakpointGenotyper male_genotyper(p_male_genotype_parameter);
    for (const auto& breakpointname : breakpoint_names)
    {
        size_t sample_index = 0;
        for (const auto& samplename : sampleNames())
        {
            auto const& depth_readlength = getDepthAndReadlength(sample_index);
            vector<int32_t> counts;
            counts.reserve(allelenames.size());
            for (const auto& e : allelenames)
            {
                counts.push_back(getCount(sample_index, breakpointname, e));
            }
            auto sample_ploidy = getSamplePloidy(sample_index);
            double expected_depth = depth_readlength.first * ((double)sample_ploidy / female_ploidy_);
            double depth_sd = getDepthSD(sample_index);
            const BreakpointGenotyperParameter b_param(
                expected_depth, depth_readlength.second, depth_sd, p_genotype_parameter->usePoissonDepth());

            if (getSamplePloidy(sample_index) == male_ploidy_)
            {
                const auto gt = male_genotyper.genotype(b_param, counts);
                setGenotype(samplename, breakpointname, gt);
            }
            else // treat unknown as female
            {
                const auto gt = genotyper.genotype(b_param, counts);
                setGenotype(samplename, breakpointname, gt);
            }
            ++sample_index;
        }
    }

    // compute combined genotype
    size_t sample_index = 0;
    for (const auto& samplename : sampleNames())
    {
        GenotypeSet all_breakpoint_gts;
        for (const auto& breakpointname : breakpoint_names)
        {
            all_breakpoint_gts.add(allelenames, getGenotype(samplename, breakpointname));
        }
        auto const& depth_readlength = getDepthAndReadlength(sample_index);
        auto depth_sd = getDepthSD(sample_index);
        const BreakpointGenotyperParameter b_param(
            depth_readlength.first, depth_readlength.second, depth_sd, p_genotype_parameter->usePoissonDepth());
        setGenotype(samplename, "", combinedGenotype(all_breakpoint_gts, &b_param, &genotyper));
        ++sample_index;
    }
}

unsigned int GraphBreakpointGenotyper::getSamplePloidy(size_t sample_index)
{
    if (getSampleSex(sample_index) == SampleInfo::MALE)
    {
        return male_ploidy_;
    }
    else // treat unknown as female
    {
        return female_ploidy_;
    }
}
}
