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
 * \brief Basic statistical functions
 *
 * \file basics.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include "common/Error.hh"

namespace statistics
{
namespace basics
{

    /**
     * Compute the median of an STL container of numeric types
     * @tparam T an STL container of a numeric type
     * @param nums an STL container of a numeric type
     * @return the median of the container
     */
    template <typename T> double median(const T& nums)
    {
        if (nums.size() == 1)
        {
            return static_cast<double>(*(nums.cbegin()));
        }

        // this circumvents the fact that std::list does not use
        // std::sort, but has its own member sort function
        std::vector<typename T::value_type> nums_copy(nums.cbegin(), nums.cend());

        // this could be reimplemented with quickselect
        std::sort(nums_copy.begin(), nums_copy.end());

        std::size_t size = nums_copy.size();

        if (size % 2 != 0)
        {
            auto mid_it = std::next(nums_copy.cbegin(), size / 2);
            return static_cast<double>(*(mid_it));
        }

        auto mid_it = std::next(nums_copy.cbegin(), (size / 2) - 1);
        auto m = static_cast<double>(*(mid_it));
        m += static_cast<double>(*(++mid_it));
        return m / 2;
    }

    /**
     * Compute the mean of an STL container of numeric types
     * @tparam T an STL container of a numeric type
     * @param nums an STL container of a numeric type
     * @return the mean of the container
     */
    template <typename T> double mean(const T& nums)
    {
        return std::accumulate(nums.cbegin(), nums.cend(), 0.0) / nums.size();
    }

    /**
     * Compute the sample variance of an STL container of numeric types
     * @tparam T an STL container of a numeric type
     * @param nums a container of a numeric type
     * @return the variance of the container
     */
    template <typename T> double var(const T& nums)
    {
        double n_mean = mean(nums);

        double var = 0;

        for (auto it = nums.cbegin(); it != nums.cend(); ++it)
        {
            var += pow(*it - n_mean, 2);
        }

        return var / (nums.size() - 1);
    }

    /**
     * Compute the mean and variance of an STL container of numeric types in a single pass
     * Taken from: http://mathworld.wolfram.com/SampleVarianceComputation.html
     * @tparam an STL container of a numeric type
     * @param nums an STL container of a numeric type
     * @return a pair of doubles (mean, variance)
     */
    template <typename T> std::pair<double, double> one_pass_mean_var(const T& nums)
    {
        assert(nums.size() > 1);

        auto it = nums.cbegin();

        auto mean = static_cast<double>(*(it++));
        double var = 0;

        for (double count = 2; it != nums.end(); ++it, ++count)
        {
            double last_mean = mean;

            mean = last_mean + ((static_cast<double>(*it) - last_mean) / count);

            double v1 = var * (1 - (1 / (count - 1)));
            double v2 = count * (pow(mean - last_mean, 2));
            var = v1 + v2;
        }

        return std::make_pair(mean, var);
    }

    /**
     * Calculate the z-scores for a set of numbers compared to a normal distribution
     * defined by mean and variance
     * @tparam T an STL container of a numeric type
     * @param nums an STL container of a numeric type
     * @param mean the distribution mean
     * @param variance the distribution variance
     * @return a vector of doubles
     */
    template <typename T> std::vector<double> zscore(const T& nums, double mean, double variance)
    {
        assert(variance > 0);
        std::vector<double> zscores;
        std::transform(
            nums.cbegin(), nums.cend(), std::back_inserter(zscores), [mean, variance](typename T::value_type x) {
                return (static_cast<double>(x) - mean) / std::sqrt(variance);
            });

        return zscores;
    }

    /**
     * Return all indices with the minimum element in an STL container of numeric types
     * @tparam T an STL container of a numeric type
     * @param nums an STL container of a numeric type
     * @return unsigned integer index of minimum value in container
     */
    template <typename T> std::vector<unsigned> min_element_indices(const T& nums)
    {
        auto min = *std::min_element(nums.cbegin(), nums.cend());

        std::vector<unsigned> min_element_indices;

        auto it = nums.cbegin();
        while ((it = std::find(it, nums.end(), min)) != nums.end())
        {
            auto index = static_cast<unsigned>(std::distance(nums.cbegin(), it));
            min_element_indices.emplace_back(index);
            ++it;
        }

        assert(!min_element_indices.empty());
        return min_element_indices;
    }
}
}
