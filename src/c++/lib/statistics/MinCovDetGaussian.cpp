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
 * \brief Implementation for MinCovDetGaussian class
 *
 * \file MinCovDetGaussian.cpp
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <vector>

#include "boost/math/distributions/chi_squared.hpp"

#include "statistics/Basics.hh"
#include "statistics/MinCovDetGaussian.hh"

#include "common/Error.hh"

namespace statistics
{

/**
 * Set the minimum allowable support fraction for a vector of doubles
 * @param nums a vector of doubles
 */
double min_support_fraction(const std::vector<double>& nums) { return (1 / static_cast<double>(nums.size())) + 0.5; }

/**
 * Get support count from support fraction
 * @param nums a vector of doubles
 * @param support_frac the support fraction to fit the MCD
 * @return an unsigned integer
 */
unsigned get_n_support(const std::vector<double>& nums, double support_frac)
{
    return static_cast<unsigned>(ceil(nums.size() * support_frac));
}

struct MinCovDetGaussian::MinCovDetGaussianImpl
{
    double raw_support_frac;
    double raw_mean, raw_variance;
    std::vector<double> support;

    MinCovDetGaussianImpl()
    {
        reset();
        raw_support_frac = std::numeric_limits<double>::quiet_NaN();
    }

    void reset()
    {
        raw_mean = std::numeric_limits<double>::quiet_NaN();
        raw_variance = std::numeric_limits<double>::quiet_NaN();
        support.clear();
    }
};

/**
 * Minimum Covariance Determinant fit constructor using default
 * raw support fraction (0.5)
 */
MinCovDetGaussian::MinCovDetGaussian()
    : AbstractGaussian()
    , _impl(new MinCovDetGaussianImpl())
{
    _impl->raw_support_frac = 0.5;
}

/**
 * Minimum Covariance Determinant fit constructor using custom raw support fraction
 * @param raw_support_frac double between 0.5 and 1
 */
MinCovDetGaussian::MinCovDetGaussian(double raw_support_frac)
    : AbstractGaussian()
    , _impl(new MinCovDetGaussianImpl())
{
    raw_support_fraction(raw_support_frac);
}

/**
 * MinCovDetGaussian move constructor
 * @param rhs another MinCovDetGaussian instance
 */
MinCovDetGaussian::MinCovDetGaussian(MinCovDetGaussian&& rhs) noexcept
    : AbstractGaussian(std::move(rhs))
    , _impl(std::move(rhs._impl))
{
}

/**
 * assignment operator overload for moving a MinCovDetGaussian instance
 * @param rhs another MinCovDetGaussian instance
 * @return a pointer to a MinCovDetGaussian instance
 */
MinCovDetGaussian& MinCovDetGaussian::operator=(MinCovDetGaussian&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    AbstractGaussian::operator=(std::move(rhs));
    return *this;
}

MinCovDetGaussian::~MinCovDetGaussian() = default;

/**
 * Get the raw variance after a distribution has been fit
 * @return double
 */
double MinCovDetGaussian::raw_mean() const { return _impl->raw_mean; }

/**
 * Get the raw variance after a distribution has been fit
 * @return double
 */
double MinCovDetGaussian::raw_variance() const { return _impl->raw_variance; }

/**
 * Get the support after a distribution has been fit
 * @return vector of doubles
 */
std::vector<double> MinCovDetGaussian::support() const { return _impl->support; }

/**
 * Reset the raw support fraction
 * @param raw support_frac double between 0.5 and 1
 */
void MinCovDetGaussian::raw_support_fraction(double support_frac)
{
    assert(support_frac > 0.5 && support_frac <= 1);
    _impl->raw_support_frac = support_frac;
}

/**
 * Get the raw support fraction
 * N.B. The raw support fraction does not take into account the
 * length of the vector being fit
 * @return double between 0.5 and 1
 */
double MinCovDetGaussian::raw_support_fraction() const { return _impl->raw_support_frac; }

/**
 * Fit the observed data using the 1D Minimum covariance determinant method
 * @param nums a vector of doubles
 */
void MinCovDetGaussian::fit(const std::vector<double>& nums)
{
    _impl->reset();
    double support_frac = std::max(min_support_fraction(nums), _impl->raw_support_frac);
    unsigned n_support = get_n_support(nums, support_frac);

    std::vector<double> centered_nums(nums);

    if (n_support < nums.size())
    {
        std::vector<double> sorted_nums(nums);
        std::sort(sorted_nums.begin(), sorted_nums.end());

        std::vector<double> diff;

        auto it1 = sorted_nums.cbegin() + n_support;
        auto it2 = sorted_nums.cbegin();
        const auto it2_end = sorted_nums.begin() + (nums.size() - n_support);

        for (; it1 != sorted_nums.cend() && it2 != it2_end; ++it1, ++it2)
        {
            diff.emplace_back(*(it1) - *(it2));
        }

        std::vector<unsigned> shortest_half_indices = basics::min_element_indices(diff);
        std::vector<double> mean_intermediate;

        for (auto i : shortest_half_indices)
        {
            mean_intermediate.emplace_back(*(sorted_nums.cbegin() + i) + *(sorted_nums.cbegin() + i + n_support));
        }

        _impl->raw_mean = 0.5 * basics::mean(mean_intermediate);

        std::transform(
            centered_nums.begin(), centered_nums.end(), centered_nums.begin(),
            std::bind(std::minus<double>(), std::placeholders::_1, _impl->raw_mean));

        std::vector<size_t> indices(nums.size());
        std::size_t idx = 0;
        std::generate(indices.begin(), indices.end(), [&idx] { return idx++; });
        std::sort(indices.begin(), indices.end(), [&centered_nums](std::size_t i1, std::size_t i2) {
            return std::abs(centered_nums[i1]) < std::abs(centered_nums[i2]);
        });
        indices.resize(n_support);

        std::vector<double> var_nums;
        for (auto i : indices)
        {
            var_nums.emplace_back(nums[i]);
        }

        _impl->raw_variance = basics::var(var_nums);

        // Correct variance for N
        _impl->raw_variance *= (n_support - 1);
        _impl->raw_variance /= n_support;
    }
    else
    {
        std::pair<double, double> mean_var = basics::one_pass_mean_var(nums);
        mean(mean_var.first);
        _impl->raw_variance = mean_var.second;

        std::for_each(
            centered_nums.begin(), centered_nums.end(), std::bind(std::minus<double>(), std::placeholders::_1, mean()));
    }

    // Calculate squared Z-scores for observations (simplified from squared Mahalanobis distance in SciKitLearn)
    std::vector<double> zscores = basics::zscore(nums, _impl->raw_mean, _impl->raw_variance);
    std::transform(zscores.begin(), zscores.end(), zscores.begin(), [](double x) { return pow(x, 2); });

    static double correction_chi2_factor
        = boost::math::quantile(boost::math::complement(boost::math::chi_squared_distribution<double>(1.), 0.5));

    double correction = basics::median(zscores) / correction_chi2_factor;

    //    double corrected_covariance = _impl->raw_variance * correction;
    std::transform(
        zscores.begin(), zscores.end(), zscores.begin(),
        std::bind(std::divides<double>(), std::placeholders::_1, correction));

    static double mask_chi2_filter
        = boost::math::quantile(boost::math::complement(boost::math::chi_squared_distribution<double>(1.), 0.025));

    for (std::size_t i = 0; i < zscores.size(); ++i)
    {
        if (std::abs(zscores[i]) < mask_chi2_filter)
        {
            _impl->support.emplace_back(nums[i]);
        }
    }

    std::pair<double, double> mean_var = basics::one_pass_mean_var(_impl->support);
    mean(mean_var.first);
    variance(mean_var.second * (_impl->support.size() - 1) / _impl->support.size());
}
}
