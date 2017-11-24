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
 * \brief Genetics helper functions
 *
 * \file Genetics.hh
 * \author Richard Shaw
 * \email rshaw@illumina.com
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cctype>
#include <string>

namespace common
{
/**
 * @brief Is base a valid uppercase non-N base?
 *
 * @param base the base
 *
 */
static inline bool isRealBase(char base) { return ((base == 'A') || (base == 'C') || (base == 'G') || (base == 'T')); }

/**
 * @brief Return transition partner of base?
 *
 * @param base the base
 *
 */
static inline char transitionBase(char base)
{
    switch (base)
    {
    case 'A':
        return 'G';
    case 'C':
        return 'T';
    case 'G':
        return 'A';
    case 'T':
        return 'C';
    }

    return 'N';
}

/**
 * @brief Report SNV type
 *
 * @param refBase the reference base
 * @param altBase the variant base
 *
 * Return true if SNV is transversion, false if transition
 *
 */
static inline bool snvIsTransversion(char refBase, char altBase, bool& isValidSnv)
{
    assert(altBase != refBase);
    const char refBaseUpper(toupper(refBase));
    const char altBaseUpper(toupper(altBase));
    isValidSnv = true;

    if (!(isRealBase(refBaseUpper) && isRealBase(altBaseUpper)))
    {
        isValidSnv = false;
        return false;
    }

    return (altBaseUpper != transitionBase(refBaseUpper));
}

/**
 * @brief return the reverse complement of a string
 * @param in input DNA string
 * @return reverse-complement of string
 */
static inline std::string reverseComplement(std::string const& in)
{
    std::string result;
    result.resize(in.size());
    std::transform(in.crbegin(), in.crend(), result.begin(), [](char v) {
        const bool is_uppercase = isupper(v);
        switch (toupper(v))
        {
        case 'A':
            return is_uppercase ? 'T' : 't';
        case 'T':
            return is_uppercase ? 'A' : 'a';
        case 'U':
            return is_uppercase ? 'A' : 'a';
        case 'G':
            return is_uppercase ? 'C' : 'c';
        case 'C':
            return is_uppercase ? 'G' : 'g';
        case 'Y':
            return is_uppercase ? 'R' : 'r';
        case 'R':
            return is_uppercase ? 'Y' : 'y';
        case 'S':
            return is_uppercase ? 'S' : 's';
        case 'W':
            return is_uppercase ? 'W' : 'w';
        case 'K':
            return is_uppercase ? 'M' : 'm';
        case 'M':
            return is_uppercase ? 'K' : 'k';
        case 'B':
            return is_uppercase ? 'V' : 'v';
        case 'D':
            return is_uppercase ? 'H' : 'h';
        case 'H':
            return is_uppercase ? 'D' : 'd';
        case 'V':
            return is_uppercase ? 'B' : 'b';
        case 'N':
            return is_uppercase ? 'N' : 'n';
        default:
            return v;
        }
    });
    return result;
}
}
