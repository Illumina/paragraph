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
 * \brief Alignment interface
 *
 * \file Alignment.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "variant/RefVar.hh"

namespace common
{

struct AlignmentParameters
{
    AlignmentParameters(
        const int8_t match = 2, const int8_t mismatch = -2, const int8_t gap_open = 3, const int8_t gap_extension = 1)
    {
        //  a   c   g   t   n
        const int8_t default_substitution_scores[] = {
            match,    mismatch, mismatch, mismatch, 0, // a
            mismatch, match,    mismatch, mismatch, 0, // c
            mismatch, mismatch, match,    mismatch, 0, // g
            mismatch, mismatch, mismatch, match,    0, // t
            0,        0,        0,        0,        0, // n
        };

        memcpy(subs_mat, default_substitution_scores, 25 * sizeof(int8_t));
        gapo = gap_open;
        gape = gap_extension;
    }

    int8_t subs_mat[25];

    uint8_t gapo;
    uint8_t gape;

    int8_t maxScore() const
    {
        int8_t max_score = std::numeric_limits<int8_t>::min();

        for (int i = 0; i < 25; ++i)
        {
            max_score = std::max(max_score, subs_mat[i]);
        }
        return max_score;
    }
};

/**
 * @brief Sequence alignment class interface
 */
class Alignment
{
public:
    virtual ~Alignment();

    virtual void setParameters(AlignmentParameters const& ap) = 0;
    virtual void getParameters(AlignmentParameters& ap) = 0;

    /* get the best possible score for comparing two sequences of length len */
    virtual int bestScore(int len);

    /**
     * @brief set target sequence
     */
    virtual void setRef(const char* seq) = 0;
    /**
     * @brief set query sequence
     */
    virtual void setQuery(const char* seq) = 0;
    /**
     * @brief Get the alignment score
     */
    virtual int getScore() = 0;
    /**
     * @brief Get a human-readable cigar string + start + end
     */
    virtual void getCigar(int& r0, int& r1, int& a0, int& a1, std::string& cig) = 0;

    /**
     * @brief Get int* cigar string + start + end
     */
    virtual void getCigar(int& r0, int& r1, int& a0, int& a1, int& n_cigar, uint32_t*& cigar) = 0;

    /**
     * @brief Debug dump, optional
     */
    virtual void dump()
    {
        int r0, r1, a0, a1;
        std::string cig;
        getCigar(r0, r1, a0, a1, cig);
        std::cerr << "Score: " << getScore() << "\n"
                  << "ref: " << r0 << "-" << r1 << "\n"
                  << "alt: " << a0 << "-" << a1 << "\n"
                  << "Cigar: " << cig << "\n";
    }
};

/** Alignment result storage */
struct AlignmentResult
{
    explicit AlignmentResult(Alignment* aln)
    {
        score = aln->getScore();
        uint32_t* tmp;
        aln->getCigar(s1, e1, s2, e2, n_cigar, tmp);
        cigar = new uint32_t[n_cigar];
        memcpy(cigar, tmp, sizeof(uint32_t) * n_cigar);
    }

    AlignmentResult()
        : cigar(NULL)
    {
    }
    ~AlignmentResult()
    {
        if (cigar)
        {
            delete[] cigar;
        }
    }
    AlignmentResult(AlignmentResult const& rhs)
        : score(rhs.score)
        , hap1(rhs.hap1)
        , hap2(rhs.hap2)
        , s1(rhs.s1)
        , e1(rhs.e1)
        , s2(rhs.s2)
        , e2(rhs.e2)
        , n_cigar(rhs.n_cigar)
    {
        if (rhs.n_cigar > 0)
        {
            cigar = new uint32_t[rhs.n_cigar];
            memcpy(cigar, rhs.cigar, sizeof(uint32_t) * rhs.n_cigar);
        }
        else
        {
            cigar = NULL;
        }
    }

    AlignmentResult& operator=(AlignmentResult const& rhs)
    {
        if (&rhs == this)
        {
            return *this;
        }
        if (cigar)
        {
            delete[] cigar;
        }
        score = rhs.score;
        hap1 = rhs.hap1;
        hap2 = rhs.hap2;
        s1 = rhs.s1;
        e1 = rhs.e1;
        s2 = rhs.s2;
        e2 = rhs.e2;
        n_cigar = rhs.n_cigar;
        if (rhs.n_cigar > 0)
        {
            cigar = new uint32_t[rhs.n_cigar];
            memcpy(cigar, rhs.cigar, sizeof(uint32_t) * rhs.n_cigar);
        }
        else
        {
            cigar = NULL;
        }
        return *this;
    }

    int score = 0;
    std::string hap1;
    std::string hap2;
    int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
    uint32_t* cigar = 0;
    int n_cigar = 0;
};

enum CigarCode
{
    ALIGN = 0, // 'M'
    INSERT = 1, // 'I'
    DELETE = 2, // 'D'
    SKIP = 3, // 'N' Essentially same as 'D' but not treated as a deletion.
    // Can be used for intron when aligning RNA sample against whole genome reference
    SOFT_CLIP = 4, // 'S'
    HARD_CLIP = 5, // 'H'
    PAD = 6, // 'P'
    MATCH = 7, // '='
    MISMATCH = 8, // 'X'
    UNKNOWN // '?'
};

class CigarCoder
{
    union {
        struct
        {
            CigarCode opCode_ : 4;
            uint32_t length_ : 28;
        };
        uint32_t value_;
    };

public:
    static const uint32_t LENGTH_MAX = ~(~uint32_t(0) << 28);
    CigarCoder(uint64_t length, CigarCode opCode)
        : opCode_(opCode)
        , length_(LENGTH_MAX < length ? LENGTH_MAX : length)
    {
        assert(ALIGN <= opCode && UNKNOWN >= opCode); //, "Invalid CIGAR code " << opCode);
        assert(length == length_); //, "Supplied length value does not fit the length_ data field: " << length << "
                                   // length_:" << length_);
    }

    explicit CigarCoder(uint32_t value)
        : value_(value)
    {
    }

    uint32_t getValue() const { return value_; }
    uint32_t getLength() const { return length_; }
    void dec(uint32_t by) { length_ -= by; }
    CigarCode getCode() const { return opCode_; }
};

typedef std::vector<uint32_t> Cigar;

template <typename CigarIT> std::string cigarToString(CigarIT pathCigarIt, const CigarIT pathCigarEnd)
{
    static const std::vector<char> CIGAR_CHARS = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?' };
    std::string ret;
    while (pathCigarEnd != pathCigarIt)
    {
        const CigarCoder c(*pathCigarIt++);
        ret += std::to_string(c.getLength()) + CIGAR_CHARS[c.getCode()];
    }

    return ret;
}

std::string
makeCigarBit(std::string::const_iterator itRef, std::string::const_iterator& itSeq, std::size_t length, int& matches);

/**
 * @brief Factory interface. Caller has to delete the returned
 * pointer
 */
extern Alignment* makeAlignment(const char* type);

/**
 * @brief Format int encoded Cigar string
 *
 * @param tb for padding with "S" : begin
 * @param te for padding with "S" : end
 * @param altlen for padding with "S" : length of alternate sequence
 * @param n_cigar length of cigar
 * @param cigar int* to cigar entries
 *
 */
std::string makeCigar(int tb, int te, int altlen, int n_cigar, uint32_t* cigar);

/** make variants from a cigar string */
extern void getVariantsFromCigar(
    std::string const& ref, std::string const& alt, int r0, int a0, uint32_t* cigar, int n_cigar,
    std::list<variant::RefVar>& target);

/** get stats from a cigar string */
extern void getCigarStats(
    std::string const& ref, std::string const& alt, int r0, int a0, uint32_t* cigar, int n_cigar, int& softclipped,
    int& matches, int& mismatches, int& ins, int& del);

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param aln the alignment interface to use
 * @param vars the primitive records
 */
void realignRefVar(
    FastaFile const& f, const char* chr, variant::RefVar const& rv, Alignment* aln, std::list<variant::RefVar>& vars);

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param aln the alignment interface to use
 * @param snps the number of snps
 * @param ins the number of insertions
 * @param dels the number of deletions
 * @param homref the number of calls with no variation
 */
void realignRefVar(
    FastaFile const& f, const char* chr, variant::RefVar const& rv, Alignment* aln, size_t& snps, size_t& ins,
    size_t& dels, size_t& homref, size_t& transitions, size_t& transversions);
}
