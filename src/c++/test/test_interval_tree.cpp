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
 *  \brief Interval tree test cases
 *
 *
 * \file test_haplotypes.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "IntervalTree.h"
#include "common/Error.hh"
#include "common/Fasta.hh"
#include "gtest/gtest.h"

#include "common.hh"

#include <assert.h>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>
#include <time.h>

using namespace std;

typedef Interval<bool> interval;
typedef vector<interval> intervalVector;
typedef IntervalTree<bool> intervalTree;
typedef vector<std::size_t> countsVector;

template <typename K> K randKey(K floor, K ceiling)
{
    K range = ceiling - floor;
    return floor + range * ((double)rand() / (double)(RAND_MAX + 1.0));
}

template <class T, typename K> Interval<T, K> randomInterval(K maxStart, K maxLength, K maxStop, const T& value)
{
    K start = randKey<K>(0, maxStart);
    K stop = min<K>(randKey<K>(start, start + maxLength), maxStop);
    return Interval<T, K>(start, stop, value);
}

TEST(IntervalTree, testIntervalTreeSanitySimple)
{
    // a simple sanity check
    intervalVector sanityIntervals;
    sanityIntervals.push_back(interval(60, 80, true));
    sanityIntervals.push_back(interval(20, 40, true));
    intervalTree sanityTree(sanityIntervals);

    intervalVector sanityResults;
    sanityTree.findOverlapping(30, 50, sanityResults);
    ASSERT_EQ(sanityResults.size(), (size_t)1);

    sanityResults.clear();
    sanityTree.findContained(15, 45, sanityResults);
    ASSERT_EQ(sanityResults.size(), (size_t)1);
}

TEST(IntervalTree, testIntervalTreeSpeed)
{
    srand((unsigned)time(NULL));

    intervalVector intervals;
    intervalVector queries;

    // generate a test set of target intervals
    for (int i = 0; i < 10000; ++i)
    {
        intervals.push_back(randomInterval<bool>(100000, 1000, 100000 + 1, true));
    }

    // and queries
    for (int i = 0; i < 5000; ++i)
    {
        queries.push_back(randomInterval<bool>(100000, 1000, 100000 + 1, true));
    }

    typedef chrono::high_resolution_clock Clock;
    typedef chrono::milliseconds milliseconds;

    // using brute-force search
    countsVector bruteforcecounts;
    Clock::time_point t0 = Clock::now();
    for (intervalVector::iterator q = queries.begin(); q != queries.end(); ++q)
    {
        intervalVector results;
        for (intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i)
        {
            if (i->start >= q->start && i->stop <= q->stop)
            {
                results.push_back(*i);
            }
        }
        bruteforcecounts.push_back(results.size());
    }

    Clock::time_point t1 = Clock::now();
    milliseconds ms = chrono::duration_cast<milliseconds>(t1 - t0);
    cout << "brute force:\t" << ms.count() << "ms" << endl;

    // using the interval tree
    intervalTree tree = intervalTree(intervals);
    countsVector treecounts;
    t0 = Clock::now();
    for (intervalVector::iterator q = queries.begin(); q != queries.end(); ++q)
    {
        intervalVector results;
        tree.findContained(q->start, q->stop, results);
        treecounts.push_back(results.size());
    }
    t1 = Clock::now();
    ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
    cout << "interval tree:\t" << ms.count() << "ms" << endl;

    // check that the same number of results are returned
    countsVector::iterator b = bruteforcecounts.begin();
    for (countsVector::iterator t = treecounts.begin(); t != treecounts.end(); ++t, ++b)
    {
        ASSERT_EQ(*t, *b);
    }
}
