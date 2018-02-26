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
 * \brief Kmer manipulation helpers
 *
 * \file Kmer.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <iomanip>
#include <string>

#include <boost/foreach.hpp>
#include <boost/integer.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/list/list10_c.hpp>
#include <boost/mpl/modulus.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/vector/vector50_c.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/type_traits/is_same.hpp>

#include "oligo/Nucleotides.hh"

namespace oligo
{

template <unsigned K> struct KmerBitsType
{
    typedef typename boost::uint_t<K * BITS_PER_BASE>::least BitsType;
};

typedef boost::mpl::vector1_c<size_t, 10> SUPPORTED_KMERS_10;
typedef boost::mpl::push_back<SUPPORTED_KMERS_10, boost::mpl::size_t<11>>::type SUPPORTED_KMERS_11;
typedef boost::mpl::push_back<SUPPORTED_KMERS_11, boost::mpl::size_t<12>>::type SUPPORTED_KMERS_12;
typedef boost::mpl::push_back<SUPPORTED_KMERS_12, boost::mpl::size_t<13>>::type SUPPORTED_KMERS_13;
typedef boost::mpl::push_back<SUPPORTED_KMERS_13, boost::mpl::size_t<14>>::type SUPPORTED_KMERS_14;
typedef boost::mpl::push_back<SUPPORTED_KMERS_14, boost::mpl::size_t<15>>::type SUPPORTED_KMERS_15;
typedef boost::mpl::push_back<SUPPORTED_KMERS_15, boost::mpl::size_t<16>>::type SUPPORTED_KMERS_16;
typedef boost::mpl::push_back<SUPPORTED_KMERS_16, boost::mpl::size_t<17>>::type SUPPORTED_KMERS_17;
typedef boost::mpl::push_back<SUPPORTED_KMERS_17, boost::mpl::size_t<18>>::type SUPPORTED_KMERS_18;
typedef boost::mpl::push_back<SUPPORTED_KMERS_18, boost::mpl::size_t<19>>::type SUPPORTED_KMERS_19;
typedef boost::mpl::push_back<SUPPORTED_KMERS_19, boost::mpl::size_t<20>>::type SUPPORTED_KMERS_20;

typedef SUPPORTED_KMERS_20 SUPPORTED_KMERS;

static const unsigned FIRST_SUPPORTED_KMER = boost::mpl::front<SUPPORTED_KMERS>::type::value;
static const unsigned LAST_SUPPORTED_KMER = boost::mpl::back<SUPPORTED_KMERS>::type::value;

template <class It, class End> bool isSupportedKmerLength(unsigned seedLength, boost::mpl::true_ endofvec)
{
    return false;
}

template <class It, class End> bool isSupportedKmerLength(unsigned seedLength, boost::mpl::false_)
{
    if (seedLength == boost::mpl::deref<It>::type::value)
    {
        return true;
    }
    else
    {
        typedef typename boost::mpl::next<It>::type Next;
        return isSupportedKmerLength<Next, End>(seedLength, typename boost::is_same<Next, End>::type());
    }
}

inline bool isSupportedKmerLength(unsigned length)
{
    typedef boost::mpl::begin<SUPPORTED_KMERS>::type begin;
    typedef boost::mpl::end<SUPPORTED_KMERS>::type end;

    return isSupportedKmerLength<begin, end>(length, boost::is_same<begin, end>::type());
}

inline std::string supportedKmersString()
{
    std::string ret;
    for (unsigned i = FIRST_SUPPORTED_KMER; i <= LAST_SUPPORTED_KMER; ++i)
    {
        if (isSupportedKmerLength(i))
        {
            ret += std::to_string(i) + ' ';
        }
    }

    return ret;
}

template <typename BitsT, typename DerivedT> class ArithmeticOperations
{
    DerivedT& derived() { return *static_cast<DerivedT*>(this); }
    const DerivedT& derived() const { return *static_cast<const DerivedT*>(this); }

public:
    typedef BitsT BitsType;

    template <typename ShiftT> DerivedT& operator<<=(const ShiftT& lshift)
    {
        derived().bits_ <<= lshift;
        derived().bits_ &= derived().BITS_MASK();
        return derived();
    }
    template <typename ShiftT> DerivedT& operator>>=(const ShiftT& rshift)
    {
        derived().bits_ >>= rshift;
        return derived();
    }
    DerivedT& operator|=(const DerivedT& that)
    {
        derived().bits_ |= that.bits_;
        return derived();
    }
    DerivedT& operator&=(const DerivedT& that)
    {
        derived().bits_ &= that.bits_;
        return derived();
    }
    DerivedT& operator^=(const DerivedT& that)
    {
        derived().bits_ ^= that.bits_;
        return derived();
    }

    template <typename ShiftT> DerivedT operator<<(const ShiftT& lshift) const
    {
        DerivedT ret(derived());
        ret <<= lshift;
        return ret;
    }
    template <typename ShiftT> DerivedT operator>>(const ShiftT& rshift) const
    {
        DerivedT ret(derived());
        ret >>= rshift;
        return ret;
    }
    DerivedT operator|(const DerivedT& right) const
    {
        DerivedT ret(derived());
        ret |= right;
        return ret;
    }
    DerivedT operator&(const DerivedT& right) const
    {
        DerivedT ret(derived());
        ret &= right;
        return ret;
    }

    DerivedT operator^(const DerivedT& right) const
    {
        DerivedT ret(derived());
        ret ^= right;
        return ret;
    }

    unsigned operator&(const unsigned& right) const { return derived().bits_ & right; }

    bool operator!() const { return !derived().bits_; }
    bool operator<(const DerivedT& right) const { return derived().bits_ < right.bits_; }
    bool operator>(const DerivedT& right) const { return right.bits_ < derived().bits_; }
    bool operator!=(const DerivedT& right) const { return derived().bits_ != right.bits_; }
    bool operator==(const DerivedT& right) const { return !(derived().bits_ != right.bits_); }
    DerivedT operator~() const { return DerivedT(~derived().bits_ & DerivedT::BITS_MASK()); }
};

#pragma pack(push, 1)
template <unsigned K>
struct BasicKmerType : public ArithmeticOperations<typename KmerBitsType<K>::BitsType, BasicKmerType<K>>
{
    typedef ArithmeticOperations<typename KmerBitsType<K>::BitsType, BasicKmerType<K>> BaseT;
    typedef typename BaseT::BitsType BitsType;
    BitsType bits_;

    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = K;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;

    static BitsType BITS_MASK() { return ~BitsType(0U) >> ((sizeof(BitsType) * BITS_PER_BYTE) - KMER_BITS); }

    explicit BasicKmerType(BitsType u)
        : bits_(u)
    {
        assert((u & BITS_MASK()) == u); // && "Invalid initialization value supplied to k-mer");
    }
};
#pragma pack(pop)

typedef BasicKmerType<32> KmerType;
BOOST_STATIC_ASSERT_MSG(8 == sizeof(KmerType), "Unexpected object type size");

typedef BasicKmerType<16> ShortKmerType;
BOOST_STATIC_ASSERT_MSG(4 == sizeof(ShortKmerType), "Unexpected object type size");

typedef BasicKmerType<18> SeedKmerType;
BOOST_STATIC_ASSERT_MSG(8 == sizeof(SeedKmerType), "Unexpected object type size");

typedef BasicKmerType<8> VeryShortKmerType;
BOOST_STATIC_ASSERT_MSG(2 == sizeof(VeryShortKmerType), "Unexpected object type size");

template <typename KmerT> inline KmerT getMaxKmer(const unsigned kmerLength)
{
    return ~(~KmerT(0) << BITS_PER_BASE * kmerLength);
}

template <unsigned kmerLength, typename KmerT> struct MaxKmer
{
    // gcc 6.1 thinks ~KmerT(0) is signed -1 for KmerT being unsigned short.
    //    static const KmerT value = ~(KmerT(-1UL) << 2 * kmerLength);
    //    static const KmerT value = ~(~KmerT(0) << 2 * kmerLength);
    //    static const KmerT value = ~((KmerT(0)-KmerT(1)) << 2 * kmerLength);
    static const KmerT value = (KmerT(1) << BITS_PER_BASE * kmerLength) - 1;
};

template <typename KmerT, bool withSuffix> struct KmerTraitsImpl
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = KmerT::KMER_BASES;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;
    static const unsigned KMER_BYTES = (KMER_BITS + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
    static const std::size_t UNUSED_BITS = KMER_BYTES * BITS_PER_BYTE - KMER_BITS;

    BOOST_STATIC_ASSERT(!(KMER_BASES % 2));
    typedef BasicKmerType<KMER_BASES / 2> SuffixType;

    static const unsigned SUFFIX_BITS = KMER_BITS / 2;
};

template <typename KmerT> struct KmerTraitsImpl<KmerT, false>
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = KmerT::KMER_BASES;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;
    static const unsigned KMER_BYTES = (KMER_BITS + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
    static const std::size_t UNUSED_BITS = KMER_BYTES * BITS_PER_BYTE - KMER_BITS;
};

template <typename KmerT>
struct KmerTraits : public KmerTraitsImpl<
                        KmerT,
                        boost::mpl::equal_to<
                            boost::mpl::modulus<boost::mpl::int_<KmerT::KMER_BASES>, boost::mpl::int_<2>>,
                            boost::mpl::int_<0>>::value>
{
};

template <unsigned bytes> struct IntegralKmerTraits
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BITS = bytes * BITS_PER_BYTE;
    static const unsigned KMER_BYTES = bytes;
    static const unsigned KMER_BASES = KMER_BITS / BITS_PER_BASE;
    static const std::size_t UNUSED_BITS = 0;

    typedef typename boost::uint_t<KMER_BITS / 2>::exact SuffixType;
    static const unsigned SUFFIX_BITS = sizeof(SuffixType) * BITS_PER_BYTE;
};

template <> struct KmerTraits<unsigned short> : public IntegralKmerTraits<sizeof(unsigned short)>
{
};

template <> struct KmerTraits<unsigned int> : public IntegralKmerTraits<sizeof(unsigned int)>
{
};

/*
 * \brief shift left by specified amount of BASES. Also ensures that the result is not undefined in
 *        case the number of bases equals to KmerTraits<KmerT>::KMER_BASES. See below.
 *
 * "The count operand can be an immediate value or register CL.
 * The count is masked to five bits, which limits the count range to 0 to 31."
 * See http://www.intel.com/design/intarch/manuals/243191.htm
 */
template <typename KmerT> KmerT shlBases(KmerT kmer, const std::size_t bases)
{
    assert(bases <= KmerTraits<KmerT>::KMER_BASES);
    //&& "Shifting more than max number of bases indicates an error in the code");
    if (KmerTraits<KmerT>::KMER_BASES == bases)
    {
        return KmerT(0);
    }
    return kmer << (oligo::BITS_PER_BASE * bases);
}

/*
 * \brief shift left by specified amount of BASES. Also ensures that the result is not undefined in
 *        case the number of bases equals to KmerTraits<KmerT>::KMER_BASES. See below.
 *
 * "The count operand can be an immediate value or register CL.
 * The count is masked to five bits, which limits the count range to 0 to 31."
 * See http://www.intel.com/design/intarch/manuals/243191.htm
 */
template <typename KmerT> KmerT safeShl(KmerT kmer, const std::size_t bits)
{
    assert(bits <= KmerTraits<KmerT>::KMER_BITS);
    // && "Shifting more than max number of bits indicates an error in the code");
    if (KmerTraits<KmerT>::KMER_BITS == bits)
    {
        return KmerT(0);
    }
    return kmer << bits;
}

template <typename KmerT> std::string bases(KmerT kmer)
{
    return bases<BITS_PER_BASE>(kmer, KmerTraits<KmerT>::KMER_BASES);
}

template <typename KmerT> inline std::string reverseBases(KmerT kmer)
{
    std::string s = bases(~kmer);
    std::reverse(s.begin(), s.end());
    return s;
}

inline unsigned char rc(const unsigned char& forward)
{
    typedef unsigned char KmerT;
    KmerT kmer(forward);
    KmerT reversed(0);
    kmer = ~kmer; // complement all the bases
    for (unsigned i = 0; sizeof(KmerT) * 8 / BITS_PER_BASE > i; ++i)
    {
        reversed <<= BITS_PER_BASE;
        reversed |= (kmer & KmerT(BITS_PER_BASE_MASK));
        kmer >>= BITS_PER_BASE;
    }
    return reversed;
}

static const unsigned char BYTE_REVERSE_COMPLEMENTS[] = {
    rc(0x00), rc(0x01), rc(0x02), rc(0x03), rc(0x04), rc(0x05), rc(0x06), rc(0x07), rc(0x08), rc(0x09), rc(0x0a),
    rc(0x0b), rc(0x0c), rc(0x0d), rc(0x0e), rc(0x0f), rc(0x10), rc(0x11), rc(0x12), rc(0x13), rc(0x14), rc(0x15),
    rc(0x16), rc(0x17), rc(0x18), rc(0x19), rc(0x1a), rc(0x1b), rc(0x1c), rc(0x1d), rc(0x1e), rc(0x1f), rc(0x20),
    rc(0x21), rc(0x22), rc(0x23), rc(0x24), rc(0x25), rc(0x26), rc(0x27), rc(0x28), rc(0x29), rc(0x2a), rc(0x2b),
    rc(0x2c), rc(0x2d), rc(0x2e), rc(0x2f), rc(0x30), rc(0x31), rc(0x32), rc(0x33), rc(0x34), rc(0x35), rc(0x36),
    rc(0x37), rc(0x38), rc(0x39), rc(0x3a), rc(0x3b), rc(0x3c), rc(0x3d), rc(0x3e), rc(0x3f), rc(0x40), rc(0x41),
    rc(0x42), rc(0x43), rc(0x44), rc(0x45), rc(0x46), rc(0x47), rc(0x48), rc(0x49), rc(0x4a), rc(0x4b), rc(0x4c),
    rc(0x4d), rc(0x4e), rc(0x4f), rc(0x50), rc(0x51), rc(0x52), rc(0x53), rc(0x54), rc(0x55), rc(0x56), rc(0x57),
    rc(0x58), rc(0x59), rc(0x5a), rc(0x5b), rc(0x5c), rc(0x5d), rc(0x5e), rc(0x5f), rc(0x60), rc(0x61), rc(0x62),
    rc(0x63), rc(0x64), rc(0x65), rc(0x66), rc(0x67), rc(0x68), rc(0x69), rc(0x6a), rc(0x6b), rc(0x6c), rc(0x6d),
    rc(0x6e), rc(0x6f), rc(0x70), rc(0x71), rc(0x72), rc(0x73), rc(0x74), rc(0x75), rc(0x76), rc(0x77), rc(0x78),
    rc(0x79), rc(0x7a), rc(0x7b), rc(0x7c), rc(0x7d), rc(0x7e), rc(0x7f), rc(0x80), rc(0x81), rc(0x82), rc(0x83),
    rc(0x84), rc(0x85), rc(0x86), rc(0x87), rc(0x88), rc(0x89), rc(0x8a), rc(0x8b), rc(0x8c), rc(0x8d), rc(0x8e),
    rc(0x8f), rc(0x90), rc(0x91), rc(0x92), rc(0x93), rc(0x94), rc(0x95), rc(0x96), rc(0x97), rc(0x98), rc(0x99),
    rc(0x9a), rc(0x9b), rc(0x9c), rc(0x9d), rc(0x9e), rc(0x9f), rc(0xa0), rc(0xa1), rc(0xa2), rc(0xa3), rc(0xa4),
    rc(0xa5), rc(0xa6), rc(0xa7), rc(0xa8), rc(0xa9), rc(0xaa), rc(0xab), rc(0xac), rc(0xad), rc(0xae), rc(0xaf),
    rc(0xb0), rc(0xb1), rc(0xb2), rc(0xb3), rc(0xb4), rc(0xb5), rc(0xb6), rc(0xb7), rc(0xb8), rc(0xb9), rc(0xba),
    rc(0xbb), rc(0xbc), rc(0xbd), rc(0xbe), rc(0xbf), rc(0xc0), rc(0xc1), rc(0xc2), rc(0xc3), rc(0xc4), rc(0xc5),
    rc(0xc6), rc(0xc7), rc(0xc8), rc(0xc9), rc(0xca), rc(0xcb), rc(0xcc), rc(0xcd), rc(0xce), rc(0xcf), rc(0xd0),
    rc(0xd1), rc(0xd2), rc(0xd3), rc(0xd4), rc(0xd5), rc(0xd6), rc(0xd7), rc(0xd8), rc(0xd9), rc(0xda), rc(0xdb),
    rc(0xdc), rc(0xdd), rc(0xde), rc(0xdf), rc(0xe0), rc(0xe1), rc(0xe2), rc(0xe3), rc(0xe4), rc(0xe5), rc(0xe6),
    rc(0xe7), rc(0xe8), rc(0xe9), rc(0xea), rc(0xeb), rc(0xec), rc(0xed), rc(0xee), rc(0xef), rc(0xf0), rc(0xf1),
    rc(0xf2), rc(0xf3), rc(0xf4), rc(0xf5), rc(0xf6), rc(0xf7), rc(0xf8), rc(0xf9), rc(0xfa), rc(0xfb), rc(0xfc),
    rc(0xfd), rc(0xfe), rc(0xff),
};
template <typename KmerT> KmerT reverseComplement(KmerT kmer)
{
    unsigned char* const begin = reinterpret_cast<unsigned char*>(&kmer);
    unsigned char* const end = begin + KmerTraits<KmerT>::KMER_BYTES;
    std::reverse(begin, end);

    BOOST_FOREACH (unsigned char& b, std::make_pair(begin, end))
    {
        b = BYTE_REVERSE_COMPLEMENTS[b];
    }

    static const std::size_t UNUSED_BITS = KmerTraits<KmerT>::UNUSED_BITS;
    if (UNUSED_BITS)
    {
        kmer >>= UNUSED_BITS;
    }

    return kmer;
}

template <typename KmerT>
inline std::ostream& operator<<(std::ostream& os, const oligo::Bases<BITS_PER_BASE, KmerT>& bases)
{
    return printBases(os, bases);
}

template <typename KmerT>
inline std::ostream& operator<<(std::ostream& os, const oligo::ReverseBases<BITS_PER_BASE, KmerT>& bases)
{
    return printReverseBases(os, bases);
}

template <typename T> std::ostream& traceHexValue(std::ostream& os, const T& kmer)
{
    const std::ios_base::fmtflags ff = os.flags();
    os << std::hex << std::setfill('0') << std::setw(sizeof(kmer) * 2) << (uint64_t)kmer;
    os.flags(ff);
    return os;
}

template <unsigned K> std::ostream& operator<<(std::ostream& os, const BasicKmerType<K>& kmer)
{
    return traceHexValue(os, kmer.bits_);
}

} // namespace oligo
