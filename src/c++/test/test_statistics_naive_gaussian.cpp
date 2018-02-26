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