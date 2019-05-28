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
