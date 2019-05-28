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
 *  \brief Phred conversion
 *
 * \file Phred.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cmath>

namespace common
{

namespace phred
{
    static inline double pow10(double d) { return pow(10, d); }

    static inline double logErrorProbToPhred(double lp) { return -10 * lp; }

    static inline double errorProbToPhred(double p) { return -10 * log10(p); }

    static inline double phredToErrorProb(double ph) { return pow10(ph / -10); }

    static inline double phredToLogErrorProb(double ph) { return ph / -10; }
}
}
