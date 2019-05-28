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
 *  \brief Test cases for haplotype comparison
 *
 *
 * \file test_refvar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common.hh"
#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "variant/RefVar.hh"

// for testing realignRefVar
#include "common/Alignment.hh"

using namespace common;
using namespace variant;

TEST(RefVar, LeftRightShiftStringSimple)
{
    const std::string ref = "AAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTT";
    RefVar r;
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    leftShift(ref, r);
    ASSERT_EQ(r.start, 17);
    ASSERT_EQ(r.end, 17);
    ASSERT_EQ(r.alt, "CGGGTTT");

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    // limit shift distance
    leftShift(ref, r, 19);
    ASSERT_EQ(r.start, 19);
    ASSERT_EQ(r.end, 19);
    ASSERT_EQ(r.alt, "GGTTTGG");

    r.start = 5;
    r.end = 6;
    r.alt = "C";

    rightShift(ref, r);
    ASSERT_EQ(r.start, 8);
    ASSERT_EQ(r.end, 9);
    ASSERT_EQ(r.alt, "C");

    // make sure we don't shift over the end
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    rightShift(ref, r);

    ASSERT_EQ(r.start, 35);
    ASSERT_EQ(r.end, 35);
    ASSERT_EQ(r.alt, "TGGGTTT");
}

TEST(RefVar, LeftShiftSimple)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";
    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    leftShift(f, "chrQ", r);
    ASSERT_EQ(r.start, 17);
    ASSERT_EQ(r.end, 17);
    ASSERT_EQ(r.alt, "CGGGTTT");

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";

    // limit shift distance
    leftShift(f, "chrQ", r, 19);
    ASSERT_EQ(r.start, 19);
    ASSERT_EQ(r.end, 19);
    ASSERT_EQ(r.alt, "GGTTTGG");
}

TEST(RefVar, ShiftNs)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 22;
    r.end = 23;
    r.alt = "C";

    leftShift(f, "chrU", r);
    ASSERT_EQ(r.start, 20);
    ASSERT_EQ(r.end, 21);
    ASSERT_EQ(r.alt, "A");

    r.start = 52;
    r.end = 52;
    r.alt = "GT";

    // limit shift distance
    rightShift(f, "chrU", r);
    ASSERT_EQ(r.start, 55);
    ASSERT_EQ(r.end, 55);
    ASSERT_EQ(r.alt, "TT");
}

TEST(RefVar, Shift_HAP_64)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 23;
    r.end = 23;
    r.alt = "C";

    leftShift(f, "chrU", r);
    ASSERT_EQ(r.start, 23);
    ASSERT_EQ(r.end, 23);
    ASSERT_EQ(r.alt, "C");

    r.start = 23;
    r.end = 23;
    r.alt = "C";

    // limit shift distance
    rightShift(f, "chrU", r);
    ASSERT_EQ(r.start, 23);
    ASSERT_EQ(r.end, 23);
    ASSERT_EQ(r.alt, "C");
}

TEST(RefVar, LeftShiftList)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    std::list<RefVar> l;

    RefVar r;

    r.start = 5;
    r.end = 6;
    r.alt = "C";
    l.push_back(r);

    r.start = 8;
    r.end = 9;
    r.alt = "C";
    l.push_back(r);

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    l.push_back(r);

    leftShift(f, "chrQ", l);

    int count = 0;
    for (RefVar const& r : l)
    {
        switch (count++)
        {
        case 0:
            ASSERT_EQ(r.start, 5);
            ASSERT_EQ(r.end, 6);
            ASSERT_EQ(r.alt, "C");
            break;
        case 1:
            ASSERT_EQ(r.start, 7);
            ASSERT_EQ(r.end, 8);
            ASSERT_EQ(r.alt, "A");
            break;
        case 2:
            ASSERT_EQ(r.start, 17);
            ASSERT_EQ(r.end, 17);
            ASSERT_EQ(r.alt, "CGGGTTT");
            break;
        default:
            break;
        }
    }
}

TEST(RefVar, RightShiftSimple)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    RefVar r;
    r.start = 5;
    r.end = 6;
    r.alt = "C";

    rightShift(f, "chrQ", r);
    ASSERT_EQ(r.start, 8);
    ASSERT_EQ(r.end, 9);
    ASSERT_EQ(r.alt, "C");

    // make sure we don't shift over the end
    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    rightShift(f, "chrQ", r);

    ASSERT_EQ(r.start, 35);
    ASSERT_EQ(r.end, 35);
    ASSERT_EQ(r.alt, "TGGGTTT");
}

TEST(RefVar, RightShiftList)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    std::list<RefVar> l;

    RefVar r;

    r.start = 5;
    r.end = 6;
    r.alt = "C";
    l.push_back(r);

    r.start = 7;
    r.end = 8;
    r.alt = "A";
    l.push_back(r);

    r.start = 30;
    r.end = 29;
    r.alt = "GGGTTT";
    l.push_back(r);

    rightShift(f, "chrQ", l);

    int count = 0;
    for (RefVar const& r : l)
    {
        switch (count++)
        {
        case 0:
            ASSERT_EQ(r.start, 6);
            ASSERT_EQ(r.end, 7);
            ASSERT_EQ(r.alt, "A");
            break;
        case 1:
            ASSERT_EQ(r.start, 8);
            ASSERT_EQ(r.end, 9);
            ASSERT_EQ(r.alt, "C");
            break;
        case 2:
            ASSERT_EQ(r.start, 35);
            ASSERT_EQ(r.end, 35);
            ASSERT_EQ(r.alt, "TGGGTTT");
            break;
        default:
            break;
        }
    }
}

TEST(RefVar, RightShiftMulti)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    // DEBUG *RV: 50-49:ATG
    // DEBUG *RV: 51-52:TG
    // DEBUG *RV: 56-57:GC
    // DEBUG *RV: 59-59:A
    // DEBUG *RV: 61-63:ACA
    // DEBUG *RV: 64-63:CCACA
    // DEBUG *RV: 67-67:C
    // DEBUG *RV: 69-70:
    // DEBUG *RV: 76-75:CAGCA
    // DEBUG *RV: 78-78:A
    // DEBUG *RV: 80-79:ACGAC
    // DEBUG *RV: 81-81:T
    // DEBUG *RV: 84-86:GTC

    RefVar r;
    std::list<RefVar> rl;
    r.start = 50;
    r.end = 49;
    r.alt = "ATG";
    rl.push_back(r);

    r.start = 51;
    r.end = 52;
    r.alt = "TG";
    rl.push_back(r);

    r.start = 56;
    r.end = 57;
    r.alt = "GC";
    rl.push_back(r);

    r.start = 59;
    r.end = 59;
    r.alt = "A";
    rl.push_back(r);

    r.start = 61;
    r.end = 63;
    r.alt = "ACA";
    rl.push_back(r);

    r.start = 64;
    r.end = 63;
    r.alt = "CCACA";
    rl.push_back(r);

    r.start = 67;
    r.end = 67;
    r.alt = "C";
    rl.push_back(r);

    r.start = 69;
    r.end = 70;
    r.alt = "";
    rl.push_back(r);

    r.start = 76;
    r.end = 75;
    r.alt = "CAGCA";
    rl.push_back(r);

    r.start = 78;
    r.end = 78;
    r.alt = "A";
    rl.push_back(r);

    r.start = 80;
    r.end = 79;
    r.alt = "ACGAC";
    rl.push_back(r);

    r.start = 81;
    r.end = 81;
    r.alt = "T";
    rl.push_back(r);

    r.start = 84;
    r.end = 86;
    r.alt = "GTC";
    rl.push_back(r);

    leftShift(f, "chrT", rl);

    const char* expected[14]
        = { "49-49:TATG", "51-52:TG",     "56-57:GC", "59-59:A",      "61-63:ACA", "64-63:CCACA", "67-67:C",
            "68-70:A",    "75-75:CCAGCA", "78-78:A",  "79-79:CACGAC", "81-81:T",   "84-86:GTC" };

    int count = 0;
    for (RefVar const& rv : rl)
    {
        ASSERT_TRUE(count < 14);
        std::ostringstream s;
        s << rv;
        ASSERT_EQ(s.str(), expected[count]);
        ++count;
    }
}

TEST(RefVar, Alleles)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    std::vector<RefVar> rl;

    // chrQ: AAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTT

    RefVar r;

    // long del
    r.start = 6;
    r.end = 11;
    r.alt = "A";
    rl.push_back(r);

    // short del
    r.start = 6;
    r.end = 9;
    r.alt = "A";
    rl.push_back(r);

    std::vector<std::string> alleles;
    toAlleles(f, "chrQ", rl, alleles);

    ASSERT_EQ(alleles.size(), (size_t)3);
    ASSERT_EQ("AAACCC", alleles[0]);
    ASSERT_EQ("A", alleles[1]);
    ASSERT_EQ("ACC", alleles[2]);
}

TEST(RefVar, AppendToVarList)
{
    std::list<RefVar> rvlist1;
    std::list<RefVar> rvlist2;
    std::list<RefVar> rvlist3;

    // Alignment matrix:
    //   0     .    :
    // REF AACACAC-AT
    //     ||||  | ||
    // PG  AATGTGCAAT
    //     ||||||  ||
    // QRY AACATGA-AT

    appendToVarList(1, 'A', 'A', rvlist1);
    appendToVarList(2, 'A', 'A', rvlist1);

    appendToVarList(3, 'C', 'T', rvlist2);
    appendToVarList(3, 'C', 'C', rvlist3);

    appendToVarList(4, 'A', 'G', rvlist2);
    appendToVarList(4, 'A', 'A', rvlist3);

    appendToVarList(5, 'C', 'T', rvlist1);
    appendToVarList(6, 'A', 'G', rvlist1);

    appendToVarList(7, 'C', 'C', rvlist2);
    appendToVarList(7, 'C', 'A', rvlist3);

    appendToVarList(8, '-', 'A', rvlist2);
    appendToVarList(8, '-', '-', rvlist3);

    appendToVarList(9, 'A', 'A', rvlist1);
    appendToVarList(10, 'T', 'T', rvlist1);

    // for(auto & x : rvlist1)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    // for(auto & x : rvlist2)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    // for(auto & x : rvlist3)
    // {
    //     std::cout << x << "; ";
    // }
    // std::cout << "\n";

    ASSERT_EQ(rvlist1.size(), (size_t)1);
    ASSERT_EQ(rvlist2.size(), (size_t)2);
    ASSERT_EQ(rvlist3.size(), (size_t)1);

    ASSERT_EQ(rvlist1.front().repr(), "5-6:TG");
    ASSERT_EQ(rvlist2.front().repr(), "3-4:TG");
    ASSERT_EQ(rvlist2.back().repr(), "8-7:A");
    ASSERT_EQ(rvlist3.front().repr(), "7-7:A");
}

TEST(RefVar, AppendToVarList2)
{
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG
    //                                                    *        *  * *   *
    // GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG------------------------------------------------------------------TTTTGC---------G

    std::string s1 = "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCCTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGC"
                     "ACTCGTCAATCAGGCTCCATTGGGAATTCCCCGAGATTCTTGTCACAGGACGG";
    std::string s2 = "GAAGTACAGAGTCGATTTGGACGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG-----------------------------"
                     "-------------------------------------TTTTGC---------G";

    std::list<RefVar> rvlist;

    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(i + 1, s1[i], s2[i], rvlist);
    }

    std::ostringstream oss;
    for (auto& x : rvlist)
    {
        oss << x << "; ";
    }
    ASSERT_EQ(oss.str(), "52-52:A; 54-59:; 61-61:C; 64-64:C; 69-134:; 136-136:T; 140-149:C; ");
}

TEST(RefVar, AppendToVarList3)
{
    std::string s1 = "GAAGTACAGAGTCGATTTG--"
                     "GACGCTTTTCGATGAGCATCGTTCACGTACCGCGAGGACCTTCAGGGTGCCCAAGAGAAGCCCTACCGCTTGGATAGCACTCGTCAATCAGGCTCCA"
                     "TTGGGAATTCCCCGAGATTCTTGTCACAGGACGG";
    std::string s2 = "GAAGTACAGAGTCGATTTGGA--CGCTTTTCGATGAGCATCGTTCACGTACCGAG------CCCACGGTG---------------------------"
                     "---------------------------------------TTTTGC---------G";

    std::list<RefVar> rvlist;

    size_t refpos = 0;
    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(refpos + 1, s1[i], s2[i], rvlist);
        if (s1[i] != '-')
        {
            ++refpos;
        }
    }

    std::ostringstream oss;
    for (auto& x : rvlist)
    {
        oss << x << "; ";
    }
    ASSERT_EQ(oss.str(), "20-21:GA; 52-52:A; 54-61:CC; 64-64:C; 69-134:; 136-136:T; 140-149:C; ");
}

TEST(RefVar, AppendToVarList4)
{
    std::string s1 = "TG--T";
    std::string s2 = "TGGAT";

    std::list<RefVar> rvlist;

    size_t refpos = 0;
    for (size_t i = 0; i < s1.size(); ++i)
    {
        appendToVarList(refpos + 1, s1[i], s2[i], rvlist);
        if (s1[i] != '-')
        {
            ++refpos;
        }
    }

    std::ostringstream oss;
    for (auto& x : rvlist)
    {
        oss << x << "; ";
    }
    ASSERT_EQ(oss.str(), "3-2:GA; ");
}

TEST(RefVar, Apply)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    RefVar rv;
    int64_t s, e;

    rv.start = 8;
    rv.end = 7;
    rv.alt = "AAACCC";
    s = 0;
    e = 15;

    std::string res = rv.apply(f, "chrT", s, e);

    // ref: ATTCTGAC------ATACAGGT
    // alt: ATTCTGACaaacccATACAGGT

    ASSERT_EQ(res, "ATTCTGACAAACCCATACAGGT");

    rv.start = 8;
    rv.end = 10;
    rv.alt = "AAACCC";
    s = 0;
    e = 15;

    // ref: ATTCTGAC------ACAGGT
    // alt: ATTCTGACaaacccACAGGT

    res = rv.apply(f, "chrT", s, e);
    ASSERT_EQ(res, "ATTCTGACAAACCCCAGGT");
}

TEST(RefVar, Primitives)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());

    //>chrS:10-15
    // CGACT

    RefVar rv;
    rv.start = 9;
    rv.end = 14;
    // complex insertion
    rv.alt = "TGCCTTT";

    std::list<RefVar> rvl;
    toPrimitives(f, "chrS", rv, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        ASSERT_EQ(oss.str(), "9-9:T; 11-11:C; 15-14:T; ");
    }

    rv.start = 9;
    rv.end = 14;
    // complex deletion
    rv.alt = "TGCC";

    rvl.clear();
    toPrimitives(f, "chrS", rv, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        ASSERT_EQ(oss.str(), "9-9:T; 11-11:C; 13-14:; ");
    }
}

TEST(RefVar, PrimitiveAlign)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());
    Alignment* aln = makeAlignment("klibg");

    //>chrS
    // CGACT

    RefVar rv;
    rv.start = 9;
    rv.end = 14;
    // complex insertion
    rv.alt = "TGCCTTT";

    std::list<RefVar> rvl;
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        ASSERT_EQ(oss.str(), "9-9:T; 11-11:C; 15-14:T; ");
    }

    rv.start = 9;
    rv.end = 14;
    // complex deletion
    rv.alt = "TGCC";

    rvl.clear();
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        ASSERT_EQ(oss.str(), "9-9:T; 11-11:C; 13-14:; ");
    }
    delete aln;
}

TEST(RefVar, PrimitiveAlign2)
{
    std::string tp = g_testenv->getBasePath() + "/../share/test-data/misc/chrQ.fa";

    FastaFile f(tp.c_str());
    Alignment* aln = makeAlignment("klibg");

    //>chrS:10
    // CGACTTGAGACATACACCTGCGCCTAATCACTTCAGAG...

    RefVar rv;
    rv.start = 9;
    rv.end = 40;
    // complex insertion
    rv.alt = "CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCA";

    std::list<RefVar> rvl;
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        // REF: CGACTTGAGACATACACCTGCGCCTAATCACTTCAGAG--------------
        //           *   *     *  *                   iiiiiiiiiiiiii
        // ALT: CGACTAGAGTCATACCCCCGCGCCTAATCACTTCACAGCTAATCACTAATCA
        ASSERT_EQ(oss.str(), "14-14:A; 18-18:T; 24-24:C; 27-27:C; 41-40:TCACAGCTAATCACTAATCA; ");
    }

    rv.start = 9;
    rv.end = 47;
    // complex deletion
    rv.alt = "CGAACCGAGACATACAGCCTACTTCACAT";

    rvl.clear();
    realignRefVar(f, "chrS", rv, aln, rvl);

    {
        std::ostringstream oss;
        for (auto& x : rvl)
        {
            oss << x << "; ";
        }
        // REF: CGACTTGAGACATACA CCTGC GCCTA ATCA CTTCAGAGG
        //         ***           ddddd       dddd      * *
        // ALT: CGAACCGAGACATACA ----- GCCTA ---- CTTCACAT-
        ASSERT_EQ(oss.str(), "12-12:A; 13-13:C; 14-14:C; 25-29:; 35-38:; 44-44:C; 46-46:T; 47-47:; ");
    }
    delete aln;
}

TEST(RefVar, Cigar2RefVar)
{
    {
        const std::string ref = "XXCYY";
        const std::string alt = "YYTZZ";
        const std::string cigar = "2S1X2S";
        int ref_left = -1, alt_left = -1;
        const std::list<RefVar> rvl = cigarToRefVar(ref, alt, cigar, ref_left, alt_left);
        ASSERT_EQ(4, ref_left);
        ASSERT_EQ(0, alt_left);
        {
            std::ostringstream oss;
            for (auto& x : rvl)
            {
                oss << x << "; ";
            }
            ASSERT_EQ(oss.str(), "0-0:T; ");
        }
    }
    {
        const std::string ref = "GGCTT";
        const std::string alt = "XXGGTTTXX";
        const std::string cigar = "2S5M2S";
        int ref_left = -1, alt_left = -1;
        const std::list<RefVar> rvl = cigarToRefVar(ref, alt, cigar, ref_left, alt_left, true);
        ASSERT_EQ(0, ref_left);
        ASSERT_EQ(0, alt_left);
        {
            std::ostringstream oss;
            for (auto& x : rvl)
            {
                oss << x << "; ";
            }
            ASSERT_EQ(oss.str(), "0-1:.; 2-2:T; 3-4:.; ");
        }
    }
    {
        const std::string ref = "CTC";
        const std::string alt = "XXXXCTCCCYYYYY";
        const std::string cigar = "4S3M2I5S";
        int ref_left = -1, alt_left = -1;
        const std::list<RefVar> rvl = cigarToRefVar(ref, alt, cigar, ref_left, alt_left);
        ASSERT_EQ(0, ref_left);
        ASSERT_EQ(0, alt_left);
        {
            std::ostringstream oss;
            for (auto& x : rvl)
            {
                oss << x << "; ";
            }
            ASSERT_EQ(oss.str(), "3-2:CC; ");
        }
    }
    {
        const std::string ref = "CCCTC";
        const std::string alt = "CTCCC";
        const std::string cigar = "2D3M2I";
        int ref_left = -1, alt_left = -1;
        const std::list<RefVar> rvl = cigarToRefVar(ref, alt, cigar, ref_left, alt_left, true);
        ASSERT_EQ(0, ref_left);
        ASSERT_EQ(0, alt_left);
        {
            std::ostringstream oss;
            for (auto& x : rvl)
            {
                oss << x << "; ";
            }
            ASSERT_EQ(oss.str(), "0-1:; 2-4:.; 5-4:CC; ");
        }
    }
    {
        const std::string ref = "CGACTTGAGACATACACCTGCGCCTAATCACTTCAGAGG";
        const std::string alt = "CGAACCGAGACATACAGCCTACTTCACAT";
        const std::string cigar = "16M5D5M4D8M1D";
        int ref_left = -1, alt_left = -1;
        const std::list<RefVar> rvl = cigarToRefVar(ref, alt, cigar, ref_left, alt_left);
        ASSERT_EQ(0, ref_left);
        ASSERT_EQ(0, alt_left);
        {
            std::ostringstream oss;
            for (auto& x : rvl)
            {
                oss << x << "; ";
            }
            // REF: CGACTTGAGACATACA CCTGC GCCTA ATCA CTTCAGAGG
            //         ***           ddddd       dddd      * *
            // ALT: CGAACCGAGACATACA ----- GCCTA ---- CTTCACAT-
            ASSERT_EQ(oss.str(), "3-3:A; 4-4:C; 5-5:C; 16-20:; 26-29:; 35-35:C; 37-37:T; 38-38:; ");
        }
    }
}
