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
 * \brief Fits a univariate Gaussian distribution using the Minimum Covariance Determinant method
 * Described here: Rousseeuw, P. J. and Leroy, A. M. (2005) References, in Robust
 *                 Regression and Outlier Detection, John Wiley & Sons, chapter 4
 * Initial multivariate Gaussian method described here:
 *                 A Fast Algorithm for the Minimum Covariance Determinant Estimator,
 *                 1999, American Statistical Association and the American Society
 *                 for Quality, TECHNOMETRICS
 * Implementation borrowed from here:
 * https://github.com/scikit-learn/scikit-learn/blob/a24c8b46/sklearn/covariance/robust_covariance.py#L302
 *
 * \file MinCovDetGaussian.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#pragma once

#include "statistics/AbstractGaussian.hh"

namespace statistics
{
class MinCovDetGaussian : public AbstractGaussian
{
public:
    explicit MinCovDetGaussian();
    explicit MinCovDetGaussian(double support_fraction);

    MinCovDetGaussian(MinCovDetGaussian const&) = delete;
    MinCovDetGaussian(MinCovDetGaussian&&) noexcept;

    MinCovDetGaussian& operator=(MinCovDetGaussian const&) = delete;
    MinCovDetGaussian& operator=(MinCovDetGaussian&&) noexcept;

    ~MinCovDetGaussian() override;

    /**
     * Fit the observed data using the Minimum Covariance Determinant method
     */
    void fit(const std::vector<double>&) override;

    /**
     * Reset the support fraction for the distribution fit
     * @param support_fraction the percentage of observations to consider when estimating distribution parameters
     */
    void raw_support_fraction(double);

    /**
     * Get the current raw support fraction
     * @return a double in the range (0.5, 1]
     */
    double raw_support_fraction() const;

    /**
     * Get the raw variance of the MCD estimate
     * @return double
     */
    double raw_mean() const;

    /**
     * Get the raw variance of the MCD estimate
     * @return double
     */
    double raw_variance() const;

    /**
     * Get the support for the final MCD estimate
     * @return vector of doubles
     */
    std::vector<double> support() const;

private:
    struct MinCovDetGaussianImpl;
    std::unique_ptr<MinCovDetGaussianImpl> _impl;
};
}
