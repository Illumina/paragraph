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
 * \brief A base class that fits a univariate Gaussian distribution using sample
 * mean / variance
 *
 * \file NaiveGaussian.hh
 * \author Mitch Bekritsky
 * \email mbekritsky@illumina.com
 *
 */

#pragma once

#include "statistics/AbstractGaussian.hh"
#include <vector>

namespace statistics
{

class NaiveGaussian : public AbstractGaussian
{
public:
    explicit NaiveGaussian()
        : AbstractGaussian()
    {
    }

    NaiveGaussian(NaiveGaussian&& rhs) noexcept
        : AbstractGaussian(std::move(rhs))
    {
    }
    NaiveGaussian(NaiveGaussian const&) = delete;

    NaiveGaussian& operator=(NaiveGaussian&& rhs) noexcept
    {
        AbstractGaussian::operator=(std::move(rhs));
        return *this;
    }

    NaiveGaussian& operator=(NaiveGaussian const&) = delete;

    /**
     * Fit the observed data with a naive Gaussian (i.e. simple mean and variance)
     */
    void fit(const std::vector<double>& nums) override;
};
}