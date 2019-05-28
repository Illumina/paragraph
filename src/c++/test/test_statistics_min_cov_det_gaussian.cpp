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
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "boost/math/distributions/normal.hpp"
#include "statistics/Basics.hh"
#include "statistics/MinCovDetGaussian.hh"

using namespace statistics;
using namespace testing;

TEST(GaussianFitStatistics, MinCovDetGaussian)
{
    std::vector<double> numbers{ 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 10.0 };
    std::vector<double> support(numbers.cbegin(), numbers.cend() - 1);

    MinCovDetGaussian res;

    ASSERT_EQ(res.raw_support_fraction(), 0.5);

    res.fit(numbers);

    ASSERT_NEAR(res.raw_mean(), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(res.mean(), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(res.raw_variance(), 0.04, ABS_ERROR_TOL);
    ASSERT_NEAR(res.variance(), 0.1, ABS_ERROR_TOL);
    ASSERT_THAT(res.support(), ElementsAreArray(support));

    std::vector<double> numbers2{ 9.8, 7.5, 6.4, 8.5, 5.5, 1.1, 7.4, 8.9 };
    std::vector<double> support2{ 9.8, 7.5, 6.4, 8.5, 5.5, 7.4, 8.9 };

    res.fit(numbers2);

    ASSERT_NEAR(res.raw_mean(), 7.65, ABS_ERROR_TOL);
    ASSERT_NEAR(res.mean(), 7.71428571, ABS_ERROR_TOL);
    ASSERT_NEAR(res.raw_variance(), 0.7784, ABS_ERROR_TOL);
    ASSERT_NEAR(res.variance(), 1.87836735, ABS_ERROR_TOL);
    ASSERT_THAT(res.support(), ElementsAreArray(support2));
}