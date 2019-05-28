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
 * \brief Implementation for NaiveGaussian class
 *
 * \file NaiveGaussian.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#include "statistics/NaiveGaussian.hh"
#include "statistics/Basics.hh"

namespace statistics
{

/**
 * Fit the observed data with a naive Gaussian (i.e. simple mean and variance)
 */
void NaiveGaussian::fit(const std::vector<double>& nums)
{
    std::pair<double, double> mean_var = statistics::basics::one_pass_mean_var(nums);

    mean(mean_var.first);
    variance(mean_var.second);
}
}
