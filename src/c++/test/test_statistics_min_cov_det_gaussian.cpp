// -*- mode: c++; indent-tabs-mode: nil; -*-
//
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