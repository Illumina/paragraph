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
 * \file test_statistics_basics.cpp
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#include "common.hh"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <vector>

#include "statistics/Basics.hh"

using namespace statistics::basics;
using namespace testing;

class BasicStatistics : public ::testing::Test
{
public:
    std::vector<double> target_zscores;
    std::vector<double> target_zscores_uint;

    virtual void SetUp()
    {
        std::vector<double> obs_d{ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

        std::transform(obs_d.cbegin(), obs_d.cend(), std::back_inserter(target_zscores), [](double x) {
            return (x - 1.5) / std::sqrt(0.14);
        });

        std::vector<double> obs_u{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };

        std::transform(obs_u.begin(), obs_u.end(), std::back_inserter(target_zscores_uint), [](unsigned x) {
            return (static_cast<double>(x) - 4) / std::sqrt(7.5);
        });
    }
};

TEST_F(BasicStatistics, ZeroVectorUint)
{
    std::vector<unsigned> numbers{ 0, 0, 0, 0 };

    ASSERT_EQ(mean(numbers), 0.);
    ASSERT_EQ(var(numbers), 0.);

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_EQ(mean_var.first, 0.);
    ASSERT_EQ(mean_var.second, 0.);

    ASSERT_EQ(mean_var.first, mean(numbers));
    ASSERT_EQ(mean_var.second, var(numbers));

    ASSERT_THROW(zscore(numbers, 0, 0), std::runtime_error);
}

TEST_F(BasicStatistics, ZeroVectorDouble)
{
    std::vector<double> numbers{ 0., 0., 0., 0. };

    ASSERT_EQ(mean(numbers), 0.);
    ASSERT_EQ(var(numbers), 0.);

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_EQ(mean_var.first, 0.);
    ASSERT_EQ(mean_var.second, 0.);

    ASSERT_EQ(mean_var.first, mean(numbers));
    ASSERT_EQ(mean_var.second, var(numbers));
}

TEST_F(BasicStatistics, DoubleVector)
{
    std::vector<double> numbers{ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    ASSERT_NEAR(mean(numbers), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(var(numbers), 0.14, ABS_ERROR_TOL);
    ASSERT_EQ(median(numbers), 1.5);

    std::vector<double> zscores = zscore(numbers, 1.5, 0.14);
    ASSERT_THAT(zscores, ElementsAreArray(target_zscores));

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_NEAR(mean_var.first, mean(numbers), ABS_ERROR_TOL);
    ASSERT_NEAR(mean_var.second, var(numbers), ABS_ERROR_TOL);

    std::vector<double> numbers2{ 1.0, 1.2, 1.4, 1.6, 1.8 };
    ASSERT_EQ(median(numbers2), 1.4);
}

TEST_F(BasicStatistics, UintVector)
{
    std::vector<unsigned int> numbers{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };

    ASSERT_NEAR(mean(numbers), 4., ABS_ERROR_TOL);
    ASSERT_NEAR(var(numbers), 7.5, ABS_ERROR_TOL);
    ASSERT_EQ(median(numbers), 4u);

    std::vector<double> zscores = zscore(numbers, 4, 7.5);
    ASSERT_THAT(zscores, ElementsAreArray(target_zscores_uint));

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_NEAR(mean_var.first, mean(numbers), ABS_ERROR_TOL);
    ASSERT_NEAR(mean_var.second, var(numbers), ABS_ERROR_TOL);
}

TEST_F(BasicStatistics, SingleElementVector)
{
    std::vector<unsigned int> num{ 1 };

    ASSERT_EQ(mean(num), 1);
    ASSERT_EQ(median(num), 1u);
    ASSERT_EQ(std::isnan(var(num)), true);
    ASSERT_THROW(one_pass_mean_var(num), std::runtime_error);
}

TEST_F(BasicStatistics, DoubleArray)
{
    std::array<double, 6> numbers = { 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    ASSERT_NEAR(mean(numbers), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(var(numbers), 0.14, ABS_ERROR_TOL);
    ASSERT_EQ(median(numbers), 1.5);

    std::vector<double> zscores = zscore(numbers, 1.5, 0.14);
    ASSERT_THAT(zscores, ElementsAreArray(target_zscores));

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_NEAR(mean_var.first, mean(numbers), ABS_ERROR_TOL);
    ASSERT_NEAR(mean_var.second, var(numbers), ABS_ERROR_TOL);
}

TEST_F(BasicStatistics, DoubleDeque)
{
    std::deque<double> numbers = { 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    ASSERT_NEAR(mean(numbers), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(var(numbers), 0.14, ABS_ERROR_TOL);
    ASSERT_EQ(median(numbers), 1.5);

    std::vector<double> zscores = zscore(numbers, 1.5, 0.14);
    ASSERT_THAT(zscores, ElementsAreArray(target_zscores));

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_NEAR(mean_var.first, mean(numbers), ABS_ERROR_TOL);
    ASSERT_NEAR(mean_var.second, var(numbers), ABS_ERROR_TOL);
}

TEST_F(BasicStatistics, DoubleList)
{
    std::list<double> numbers = { 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    ASSERT_NEAR(mean(numbers), 1.5, ABS_ERROR_TOL);
    ASSERT_NEAR(var(numbers), 0.14, ABS_ERROR_TOL);
    ASSERT_EQ(median(numbers), 1.5);

    std::vector<double> zscores = zscore(numbers, 1.5, 0.14);
    ASSERT_THAT(zscores, ElementsAreArray(target_zscores));

    std::pair<double, double> mean_var = one_pass_mean_var(numbers);

    ASSERT_NEAR(mean_var.first, mean(numbers), ABS_ERROR_TOL);
    ASSERT_NEAR(mean_var.second, var(numbers), ABS_ERROR_TOL);
}

TEST_F(BasicStatistics, MinVectorElement)
{
    std::vector<double> d_numbers = { 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };

    ASSERT_EQ(*min_element_indices(d_numbers).cbegin(), 0u);

    d_numbers = { 1.2, 1.4, 0.0, 1.6, 1.8, 2.0 };
    ASSERT_EQ(*min_element_indices(d_numbers).cbegin(), 2u);

    std::vector<unsigned> unumbers = { 1, 2, 3, 4, 0, 5, 6, 7, 8 };
    ASSERT_EQ(*min_element_indices(unumbers).cbegin(), 4u);
}