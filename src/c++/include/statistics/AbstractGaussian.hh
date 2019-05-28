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
 * \brief Abstract base class to fit univariate Gaussian distributions
 *
 * \file AbstractGaussian.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace statistics
{

class AbstractGaussian
{
public:
    AbstractGaussian()
    {
        mean_ = std::numeric_limits<double>::quiet_NaN();
        variance_ = std::numeric_limits<double>::quiet_NaN();
    }

    AbstractGaussian(AbstractGaussian const&) noexcept = default;
    AbstractGaussian(AbstractGaussian&&) noexcept = default;

    AbstractGaussian& operator=(AbstractGaussian const&) noexcept = default;
    AbstractGaussian& operator=(AbstractGaussian&&) noexcept = default;

    virtual ~AbstractGaussian() = default;

    /**
     * Fit the observed data to a Gaussian (method specified in
     * derived class)
     */
    virtual void fit(const std::vector<double>&) = 0;

    /**
     * Get the mean of the fitted distribution
     * @return double
     */
    double mean() const { return mean_; };

    /**
     * Get the variance of the fitted distribution
     * @return double
     */
    double variance() const { return variance_; }

    /**
     * Get the standard deviation of the fitted distribution
     * @return double
     */
    double standard_deviation() const { return sqrt(variance_); };

protected:
    void mean(double mean) { mean_ = mean; }
    void variance(double variance) { variance_ = variance; }

private:
    double mean_, variance_;
};
}