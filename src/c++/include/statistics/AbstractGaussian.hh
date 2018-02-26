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