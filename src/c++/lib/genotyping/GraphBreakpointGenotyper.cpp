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
