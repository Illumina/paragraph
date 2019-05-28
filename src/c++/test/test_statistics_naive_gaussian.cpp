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
 * \file test_statistics_gaussian_fit.cpp
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#include "common.hh"
#include "gtest/gtest.h"

#include "boost/math/distributions/normal.hpp"
#include "statistics/Basics.hh"
#include "statistics/NaiveGaussian.hh"

typedef boost::math::normal Gaussian;

using namespace statistics;

TEST(GaussianFitStatistics, NaiveGaussian)
{
    std::vector<double> numbers{ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    NaiveGaussian res;
    res.fit(numbers);

    ASSERT_NEAR(res.mean(), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(res.variance(), 0.14, ABS_ERROR_TOL);
    ASSERT_NEAR(res.standard_deviation(), pow(0.14, 0.5), ABS_ERROR_TOL);

    // Test reset num
    std::vector<double> numbers2{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    res.fit(numbers2);

    ASSERT_NEAR(res.mean(), 4., ABS_ERROR_TOL);
    ASSERT_NEAR(res.variance(), 7.5, ABS_ERROR_TOL);
}