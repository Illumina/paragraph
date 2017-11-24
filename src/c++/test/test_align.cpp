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
 * \brief Alignment testing
 *
 * \file test_align.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "gtest/gtest.h"

#include "common/Alignment.hh"
#include "common/Timing.hh"

#include <cstdlib>

using namespace common;

TEST(Alignment, KlibBasic)
{
    Alignment* aln = makeAlignment("klib");

    std::cerr << "Testing klib alignment."
              << "\n";

    aln->setRef("AAATGACGGATTG");
    aln->setQuery("TGGGA");

    // AAATGACGGATTG
    //      **|||
    // -----TGGGA---

    aln->dump();

    uint32_t* icigar = NULL;
    int r0 = -1, r1 = -1, a0 = -1, a1 = -1;
    int softclipped = -1, mismatches = -1, matches = -1, ins = -1, del = -1;
    int ncigar = -1;
    aln->getCigar(r0, r1, a0, a1, ncigar, icigar);
    getCigarStats("AAATGACGGATTG", "TGGGA", r0, a0, icigar, ncigar, softclipped, matches, mismatches, ins, del);

    ASSERT_EQ(softclipped, 10);
    ASSERT_EQ(matches, 3);
    ASSERT_EQ(mismatches, 0);
    ASSERT_EQ(del, 0);
    ASSERT_EQ(ins, 0);

    delete aln;
}

TEST(Align, KlibIndel)
{
    Alignment* aln = makeAlignment("klib");
    int r0, r1, a0, a1;
    std::string cig;

    std::cerr << "Testing klib indel alignment."
              << "\n";

    // ref del / alt ins

    // AAATGAC-----GGATTG
    //
    // AAATGACCACCAGGATTG

    aln->setRef("AAATGACGGATTG");
    aln->setQuery("AAATGACCACCAGGATTG");
    aln->dump();

    aln->getCigar(r0, r1, a0, a1, cig);

    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 12);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 17);
    ASSERT_EQ(cig, "7M5D6M");

    // ref ins / alt del

    // AAATGACCACCAGGATTG
    //
    // AAATGAC-----GGATTG

    aln->setRef("AAATGACCACCAGGATTG");
    aln->setQuery("AAATGACGGATTG");
    aln->dump();

    aln->getCigar(r0, r1, a0, a1, cig);

    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 17);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 12);
    ASSERT_EQ(cig, "7M5I6M");

    // complex

    aln->setRef("AAATGACGGGGCATTGCCA");
    aln->setQuery("AAATGACCACCAGGATTGCCA");

    // AAATGAC---GGGGCATTG
    //
    // AAATGACCACCAGG-ATTG

    aln->dump();

    aln->getCigar(r0, r1, a0, a1, cig);

    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 18);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 20);
    ASSERT_EQ(cig, "9M3D2M1I7M");

    delete aln;

    aln = makeAlignment("klibg");
    // complex different scores

    AlignmentParameters ap;
    const int8_t bwa_subs_mat[] = {
        1,  -4, -4, -4, 0, // a
        -4, 1,  -4, -4, 0, // c
        -4, -4, 1,  -4, 0, // g
        -4, -4, -4, 1,  0, // t
        0,  0,  0,  0,  0, // n
    };

    memcpy(ap.subs_mat, bwa_subs_mat, 25 * sizeof(int8_t));
    ap.gapo = 6;
    ap.gape = 1;

    aln->setParameters(ap);

    aln->setRef("CCCCAAATGACGGGGCATTGCCACACA");
    aln->setQuery("CCCCAAATGACCACCAGGATTGCCACACA");

    // CCCCAAATGAC-----GGGGCATTGCCACACA
    //                   *
    // CCCCAAATGACCACCAGGA---TTGCCACACA

    aln->dump();

    aln->getCigar(r0, r1, a0, a1, cig);
    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 26);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 28);
    ASSERT_EQ(cig, "11M5D2M3I11M");

    // REF   CTGTG------TGTGTGTGTGAAAA
    // ALT1  CTGTGTGTGTGTGTGTGTGTGAAAA
    // ALT2  CTGTGTGTGAGTGTGTGTGTGAAAA

    // => Left-shifting ALT1

    // REF   C------TGTGTGTGTGTGTGAAAA
    // ALT1  CTGTGTGTGTGTGTGTGTGTGAAAA

    // C -> CTGTGTG at pos. 1

    // => Left-shifting ALT2

    // REF   CTGT------GTGTGTGTGTGAAAA
    // ALT2  CTGTGTGTGAGTGTGTGTGTGAAAA

    // T -> TGTGTGA at pos. 4

    // Ref vs. Alt1:  9M6D10M (Score 7)
    // --------------------------------
    // CTGTGTGTG------TGTGTGAAAA
    // CTGTGTGTGAGTGTGTGTGTGAAAA

    // Ref vs. Alt1:   15M6D4M (Score 7)
    // ---------------------------------
    // CTGTGTGTGTGTGT------GAAAA
    // CTGTGTGTGTGTGTGTGTGTGAAAA

    // Alt1 vs. Alt2:  25M (Score 20)
    // ------------------------------
    // CTGTGTGTGTGTGTGTGTGTGAAAA
    //          *
    // CTGTGTGTGAGTGTGTGTGTGAAAA

    aln->setRef("CTGTGTGTGTGTGTGAAAA");
    aln->setQuery("CTGTGTGTGTGTGTGTGTGTGAAAA");
    aln->dump();
    // Score: 7
    // ref: 0-18
    // alt: 0-24
    // Cigar: 15M6D4M
    aln->getCigar(r0, r1, a0, a1, cig);
    ASSERT_EQ(aln->getScore(), 7);
    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 18);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 24);
    ASSERT_EQ(cig, "15M6D4M");

    aln->setRef("CTGTGTGTGTGTGTGAAAA");
    aln->setQuery("CTGTGTGTGAGTGTGTGTGTGAAAA");
    aln->dump();
    // Score: 7
    // ref: 0-18
    // alt: 0-24
    // Cigar: 9M6D10M
    aln->getCigar(r0, r1, a0, a1, cig);
    ASSERT_EQ(aln->getScore(), 7);
    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 18);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 24);
    ASSERT_EQ(cig, "9M6D10M");

    aln->setRef("CTGTGTGTGTGTGTGTGTGTGAAAA");
    aln->setQuery("CTGTGTGTGAGTGTGTGTGTGAAAA");

    aln->dump();
    // Score: 20
    // ref: 0-24
    // alt: 0-24
    // Cigar: 25M
    aln->getCigar(r0, r1, a0, a1, cig);
    ASSERT_EQ(aln->getScore(), 20);
    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 24);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 24);
    ASSERT_EQ(cig, "25M");

    // DEBUG REF:
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG
    // DEBUG ALT: GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAGCCCACGGTGTTTTGCG
    aln->setRef("GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCG"
                "TCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG");
    aln->setQuery("GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAGCCCACGGTGTTTTGCG");
    // Score: -55
    // ref: 0-149
    // alt: 0-68
    // Cigar: 53M6I9M66I6M9I1M

    aln->dump();
    aln->getCigar(r0, r1, a0, a1, cig);
    ASSERT_EQ(aln->getScore(), -55);
    ASSERT_EQ(r0, 0);
    ASSERT_EQ(r1, 149);
    ASSERT_EQ(a0, 0);
    ASSERT_EQ(a1, 68);
    ASSERT_EQ(cig, "53M6I9M66I6M9I1M");

    uint32_t* icigar = NULL;
    int softclipped = -1, mismatches = -1, matches = -1, ins = -1, del = -1;
    int ncigar = -1;
    aln->getCigar(r0, r1, a0, a1, ncigar, icigar);

    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG
    //                                                    *        *  * *   *
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG------------------------------------------------------------------TTTTGC---------G

    getCigarStats(
        "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAG"
        "GCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG",
        "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAGCCCACGGTGTTTTGCG", r0, a0, icigar, ncigar, softclipped,
        matches, mismatches, ins, del);

    ASSERT_EQ(softclipped, 0);
    ASSERT_EQ(matches, 52 + 7 + 4 + 1);
    ASSERT_EQ(mismatches, 5);
    ASSERT_EQ(del, 66 + 6 + 9);
    ASSERT_EQ(ins, 0);

    delete aln;
}

struct AlignmentTimer
{
    AlignmentTimer(const char* type, size_t len1, size_t len2)
    {
        aln = makeAlignment(type);
        for (size_t i = 0; i < len1; ++i)
        {
            ref += chars[rand() & 3];
        }
        for (size_t i = 0; i < len2; ++i)
        {
            alt += chars[rand() & 3];
        }
    }

    virtual ~AlignmentTimer() { delete aln; }

    void operator()()
    {
        aln->setRef(ref.c_str());
        aln->setQuery(alt.c_str());
        int j = aln->getScore();
        srand(j);
    }

    Alignment* aln;
    std::string ref;
    std::string alt;

    const char chars[4] = { 'A', 'C', 'G', 'T' };
};

TEST(Align, PerformanceReads)
{
    std::cerr << "Testing Ksw performance for read realignment."
              << "\n";

    std::cerr << "\n";
    std::cerr << "size\tklib\tklibg\n";
    int count = 0;

    const int SIZE1 = 160;
    const int SIZE2 = 100;
#ifdef _DEBUG
    const int N = 1;
    const int REPS = 16;
#else
    const int N = 100;
    const int REPS = 32;
#endif
    while (count++ < REPS)
    {
        double resultk1, resultk2;
        TIMEIT(resultk1, N, AlignmentTimer, "klib", SIZE1, SIZE2);
        TIMEIT(resultk2, N, AlignmentTimer, "klibg", SIZE1, SIZE2);

        std::cerr << SIZE1 << "\t" << SIZE2 << "\t" << resultk1 << "\t" << resultk2 << "\t"
                  << "\n";
    }
    std::cerr << "\n";
}

TEST(Align, PerformanceSquares)
{
    std::cerr << "Testing Ksw performance for Haplotype comparison."
              << "\n";

    std::cerr << "\n";
    std::cerr << "size\tklib\tklibg\n";
    int size = 128;

#ifdef _DEBUG
    const int N = 2;
    const int MAX = 260;
    const int I = 128;
#else
    const int N = 100;
    const int MAX = 515;
    const int I = 64;
#endif
    while (size < MAX)
    {
        double resultk1, resultk2;
        TIMEIT(resultk1, N, AlignmentTimer, "klib", size, size);
        TIMEIT(resultk2, N, AlignmentTimer, "klibg", size, size);

        std::cerr << size << "\t" << resultk1 << "\t" << resultk2 << "\n";

        size += I;
    }
    std::cerr << "\n";
}
