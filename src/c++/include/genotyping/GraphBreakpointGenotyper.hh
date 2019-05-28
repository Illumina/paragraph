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
 * Genotyper for the graph-represented whole variants with multiple breakpoints.
 * Most likely genotype is voted from all breakpoint genotypes.
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <cstdint>

#include "GraphGenotyper.hh"
#include "genotyping/Genotype.hh"
#include "genotyping/GenotypingParameters.hh"

namespace genotyping
{
class GraphBreakpointGenotyper : public GraphGenotyper
{
public:
    GraphBreakpointGenotyper(unsigned int male_ploidy = 2, unsigned int female_ploidy = 2)
        : male_ploidy_(male_ploidy)
        , female_ploidy_(female_ploidy){};

    /**
     * set genotyping parameters from JSON
     */
    void setParameters(const std::string& genotyping_parameter_path) override;

protected:
    /**
     * Implemented by derived classes which implement genotyping
     */
    void runGenotyping() override;

private:
    /**
     * get ploidy according to its sex
     */
    unsigned int getSamplePloidy(size_t sample_index);

    /**
     * genotyping parameters
     */
    std::unique_ptr<GenotypingParameters> p_genotype_parameter;
    std::unique_ptr<GenotypingParameters> p_male_genotype_parameter;

    unsigned int male_ploidy_;
    unsigned int female_ploidy_;
};
}
