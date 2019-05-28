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
 * \brief Helper to get CPU time
 *
 * \file Timing.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include <chrono>

#define TIMEIT(result, N, ts, ...)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;              \
        using std::chrono::duration_cast;                                                                              \
                                                                                                                       \
        auto t0 = std::chrono::high_resolution_clock::now();                                                           \
        for (int n = 0; n < N; ++n)                                                                                    \
        {                                                                                                              \
            ts timeme(__VA_ARGS__);                                                                                    \
        }                                                                                                              \
        auto t1 = std::chrono::high_resolution_clock::now();                                                           \
        auto t2 = std::chrono::high_resolution_clock::now();                                                           \
        for (int n = 0; n < N; ++n)                                                                                    \
        {                                                                                                              \
            ts timeme(__VA_ARGS__);                                                                                    \
            timeme();                                                                                                  \
        }                                                                                                              \
        auto t3 = std::chrono::high_resolution_clock::now();                                                           \
        auto ticks_per_iter = Cycle((t3 - t2) - (t1 - t0)) / N;                                                        \
        result = ticks_per_iter.count();                                                                               \
    } while (0)
