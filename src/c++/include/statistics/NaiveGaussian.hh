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