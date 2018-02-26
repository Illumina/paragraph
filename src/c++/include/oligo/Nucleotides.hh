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
 * \brief Various bitwise operations for nucleotide encoding/decoding
 *
 * \file Nucleotides.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <vector>

#include "common/Error.hh"

namespace oligo
{

static const unsigned BITS_PER_BASE = 2;
static const unsigned BITS_PER_BASE_MASK = 3;

// valid ones are 0(A) 1(C) 2(G) and 3(T).
// for data, invalidOligo represents N
// for reference, invalidOligo indicates any non-ACGT base value
static const unsigned int INVALID_OLIGO = 4;

// It used to be that n in sequence would not match N in reference. This has caused problems generating NM bam tag
// changed to fix SAAC-697
static const char SEQUENCE_OLIGO_N = 'N';
static const char REFERENCE_OLIGO_N = 'N';

static const unsigned char BCL_QUALITY_MASK = 0xfc;
static const unsigned char BCL_BASE_MASK = 0x03;

/**
 * \brief Translator only ensures that no access is made outside the table 256 character space. No validation is done.
 */
template <bool withN = false, unsigned dfltVal = INVALID_OLIGO> struct Translator
{
    unsigned char operator[](const char& base) const { return getTranslatorTable()[static_cast<unsigned char>(base)]; }

private:
    static const unsigned char* getTranslatorTable()
    {
        static constexpr unsigned char table[256] = {
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            0 /*A*/,
            dfltVal,
            1 /*C*/,
            dfltVal,
            dfltVal,
            dfltVal,
            2 /*G*/,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            withN ? INVALID_OLIGO : dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            3 /*T*/,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            0 /*a*/,
            dfltVal,
            1 /*c*/,
            dfltVal,
            dfltVal,
            dfltVal,
            2 /*g*/,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            withN ? INVALID_OLIGO : dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            3 /*t*/,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
            dfltVal,
        };

        return table;

        // clang-format off
        BOOST_STATIC_ASSERT(table[(unsigned)'a'] == 0);
        BOOST_STATIC_ASSERT(table[(unsigned)'A'] == 0);
        BOOST_STATIC_ASSERT(table[(unsigned)'c'] == 1);
        BOOST_STATIC_ASSERT(table[(unsigned)'C'] == 1);
        BOOST_STATIC_ASSERT(table[(unsigned)'g'] == 2);
        BOOST_STATIC_ASSERT(table[(unsigned)'G'] == 2);
        BOOST_STATIC_ASSERT(table[(unsigned)'t'] == 3);
        BOOST_STATIC_ASSERT(table[(unsigned)'T'] == 3);
        BOOST_STATIC_ASSERT(table[(unsigned)'n'] == (withN ? INVALID_OLIGO : dfltVal));
        BOOST_STATIC_ASSERT(table[(unsigned)'N'] == (withN ? INVALID_OLIGO : dfltVal));
        // clang-format on
    }
};

inline unsigned int getValue(const char& base)
{
    static const oligo::Translator<> translator = {};
    return translator[base];
}

static const boost::array<char, 5> allBases = boost::assign::list_of('A')('C')('G')('T')('N');
inline char getBase(unsigned int base, const bool upperCase = true)
{
    base = std::min(base, static_cast<unsigned int>(allBases.size()));
    if (upperCase)
    {
        return allBases[base];
    }
    else
    {
        return allBases[base] + ('a' - 'A');
    }
}

/**
 * \return upercase base. Note that this one will not return N for bcl 0!
 */
inline char getUppercaseBase(unsigned int base) { return getBase(base, true); }

inline bool isBclN(const char bclByte) { return !(bclByte & BCL_QUALITY_MASK); }

inline unsigned char getQuality(char bclByte) { return static_cast<unsigned char>(bclByte) >> BITS_PER_BASE; }

/**
 * \return uppercase base or N for bcl N
 */
inline char getReferenceBaseFromBcl(unsigned char bcl)
{
    return oligo::isBclN(bcl) ? REFERENCE_OLIGO_N : getUppercaseBase(bcl & BCL_BASE_MASK);
}

/**
 * \return uppercase base or N for bcl N
 */
inline char getSequenceBaseFromBcl(unsigned char bcl)
{
    return oligo::isBclN(bcl) ? SEQUENCE_OLIGO_N : getUppercaseBase(bcl & BCL_BASE_MASK);
}

// Take a packed (2 bit per base) kmer and output to a string buffer
template <typename OutItr> void unpackKmer(const uint64_t kmer, const uint64_t kmerLength, OutItr itr)
{
    for (unsigned int i = 0; i < kmerLength * BITS_PER_BASE; i += BITS_PER_BASE)
    {
        *itr++ = getBase(((kmer >> i) & BCL_BASE_MASK));
    }
}

static const boost::array<char, 5> allReverseBases = boost::assign::list_of('T')('G')('C')('A')('N');
inline char getReverseBase(unsigned int base, const bool upperCase = true)
{
    //    return base > 3 ? getBase(base, upperCase) : getBase((~base) & 3, upperCase);
    base = std::min(base, static_cast<unsigned int>(allReverseBases.size()));
    if (upperCase)
    {
        return allReverseBases[base];
    }
    else
    {
        return allReverseBases[base] + ('a' - 'A');
    }
}

inline char reverseBase(const char base)
{
    switch (base)
    {
    case 'a':
        return 't';
    case 'A':
        return 'T';
    case 'c':
        return 'g';
    case 'C':
        return 'G';
    case 'g':
        return 'c';
    case 'G':
        return 'C';
    case 't':
        return 'a';
    case 'T':
        return 'A';
    case 'n':
        return 'n';
    default:
        return 'N';
    }
}

/**
 * \brief reverse-complements the base bits of a bcl byte
 *
 * \returns 0 for 0, reverse-complemented lower bits with quality bits unchanged
 */
inline unsigned char getReverseBcl(const unsigned char bcl)
{
    return !isBclN(bcl) ? (bcl & BCL_QUALITY_MASK) | (BCL_BASE_MASK - (bcl & BCL_BASE_MASK)) : 0;
}

template <int bitsPerBase, typename KmerT> std::string bases(KmerT kmer, unsigned kmerLength)
{
    static const unsigned kmerMask = ~(~0U << bitsPerBase);

    std::string s;
    while (kmerLength--)
    {
        s.push_back(getBase((kmer >> bitsPerBase * kmerLength) & kmerMask));
    }
    return s;
}

template <unsigned bitsPerBase, typename KmerT> struct Bases
{
    const KmerT kmer_;
    const unsigned kmerLength_;
    static const unsigned bitsPerBase_ = bitsPerBase;
    static const unsigned kmerMask_ = ~(~0U << bitsPerBase);
    Bases(const KmerT kmer, unsigned kmerLength)
        : kmer_(kmer)
        , kmerLength_(kmerLength)
    {
    }
};

template <unsigned bitsPerBase, typename KmerT> struct ReverseBases
{
    const KmerT kmer_;
    const unsigned kmerLength_;
    static const unsigned bitsPerBase_ = bitsPerBase;
    static const unsigned kmerMask_ = ~(~0U << bitsPerBase);
    ReverseBases(const KmerT kmer, unsigned kmerLength)
        : kmer_(kmer)
        , kmerLength_(kmerLength)
    {
    }
};

template <typename BasesT> std::ostream& printBases(std::ostream& os, const BasesT& b)
{
    unsigned pos = b.kmerLength_;
    while (pos--)
    {
        static const unsigned kmerMask_ = BasesT::kmerMask_;
        os << getBase((b.kmer_ >> BasesT::bitsPerBase_ * pos) & kmerMask_);
    }
    return os;
}

template <typename BasesT> std::ostream& printReverseBases(std::ostream& os, const BasesT& b)
{
    for (unsigned pos = 0; b.kmerLength_ > pos; ++pos)
    {
        static const unsigned kmerMask_ = BasesT::kmerMask_;
        os << getReverseBase(unsigned((b.kmer_ >> BasesT::bitsPerBase_ * pos) & kmerMask_));
    }
    return os;
}

template <typename BasesIteratorT> inline std::string bclToString(BasesIteratorT basesIterator, unsigned length)
{
    std::string ret;
    while (length--)
    {
        ret += isBclN(*basesIterator) ? SEQUENCE_OLIGO_N : getBase(BCL_BASE_MASK & *basesIterator);
        ++basesIterator;
    }
    return ret;
}

inline std::string bclToRString(const unsigned char* basesIterator, unsigned length)
{
    std::string ret;
    while (length--)
    {
        ret += oligo::isBclN(*(basesIterator + length))
            ? SEQUENCE_OLIGO_N
            : oligo::reverseBase(oligo::getBase(BCL_BASE_MASK & *(basesIterator + length)));
    }
    return ret;
}

template <typename FwdIteratorT> uint64_t pack32BclBases(FwdIteratorT bcl)
{
    uint64_t ret = 0;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++);
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 2;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 4;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 6;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 8;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 10;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 12;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 14;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 16;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 18;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 20;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 22;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 24;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 26;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 28;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 30;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 32;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 34;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 36;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 38;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 40;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 42;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 44;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 46;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 48;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 50;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 52;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 54;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 56;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 58;
    ret |= uint64_t(BCL_BASE_MASK & *bcl++) << 60;
    ret |= uint64_t(BCL_BASE_MASK & *bcl) << 62;
    return ret;
}

template <typename FwdIteratorT> uint64_t packBclBases(FwdIteratorT bclBegin, FwdIteratorT bclEnd)
{
    uint64_t ret = 0;
    const int length = bclEnd - bclBegin;
    assert(length <= 32); // && "Cannot pack more than 64 bases");

    for (int i = 0; i < length; ++i)
    {
        ret |= (BCL_BASE_MASK & *bclBegin++) << (i * BITS_PER_BASE);
    }
    // Nothing left to do with remaining: 32 - length, 0 init will be fine
    return ret;
}

} // namespace oligo
