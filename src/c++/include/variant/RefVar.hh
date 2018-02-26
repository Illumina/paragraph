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
 *  \brief Helper for storing and normalizing reference variation records
 *
 *
 * \file RefVar.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cassert>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "common/BCFHelpers.hh"
#include "common/Error.hh"
#include "common/Fasta.hh"

namespace variant
{

/** short variation record */
typedef struct _RefVar
{
    _RefVar()
        : flags(0)
    {
    }

    _RefVar(int64_t _start, int64_t _end, std::string _alt, int64_t _flags = 0)
        : start(_start)
        , end(_end)
        , alt(_alt)
        , flags(_flags)
    {
        if (alt == "<DEL>")
        {
            alt = "";
        }
        else
        {
            if (alt.find("<") != std::string::npos || alt.find(">") != std::string::npos)
            {
                error("Invalid allele at %i: %s", _start, alt.c_str());
            }
        }
    }

    int64_t start, end;
    std::string alt;

    int64_t flags;

    inline std::string repr() const { return std::to_string(start) + "-" + std::to_string(end) + ":" + alt; }

    /** apply and retrieve modified reference sequence between start and end
     *  Note that _start and _end will be extended if they aren't contained within
     *  [start, end]
     */
    inline std::string apply(common::FastaFile const& ref, std::string const& chr, int64_t& _start, int64_t& _end) const
    {
        if (_start > start)
        {
            _start = start;
        }
        if (_end < end)
        {
            _end = end;
        }

        std::string result = ref.query(chr.c_str(), _start, _end);

        // we have start >= _start
        int64_t pos = start - _start;
        int64_t reflen = end - start + 1;

        assert(reflen >= 0);

        // support insertions with no reference character
        if (reflen == 0)
        {
            result.insert(pos, alt);
        }
        else
        {
            result.replace(pos, reflen, alt);
        }

        return result;
    }

    bool operator<(const _RefVar& rhs) const { return start < rhs.start; }

    bool operator==(const _RefVar& rhs) const { return start == rhs.start && end == rhs.end && alt == rhs.alt; }

} RefVar;

static inline std::ostream& operator<<(std::ostream& o, RefVar const& r)
{
    o << r.repr();
    return o;
}

/**
 * @brief Trim reference bases at left/right end of variant
 *
 * @param f fasta file
 * @param chr chromosome in fasta file
 * @param rv RefVar record
 * @param refpadding enable single-base reference padding
 */
void trimLeft(std::string const& ref, RefVar& rv, bool refpadding = true, int64_t rel_start = 0);

void trimRight(std::string const& ref, RefVar& rv, bool refpadding = true, int64_t rel_start = 0);

void trimLeft(common::FastaFile const& f, const char* chr, RefVar& rv, bool refpadding = true);

void trimRight(common::FastaFile const& f, const char* chr, RefVar& rv, bool refpadding = true);

/**
 * @brief Left/right shifting w.r.t reference fasta
 *
 * Left/right boundary position can be given to prevent overlap with other variation
 *
 */
extern void leftShift(std::string const& ref, RefVar& rv, int64_t pos_min = -1);

extern void rightShift(std::string const& ref, RefVar& rv, int64_t pos_max = std::numeric_limits<int64_t>::max());

extern void leftShift(common::FastaFile const& f, const char* chr, RefVar& rv, int64_t pos_min = -1);

extern void rightShift(
    common::FastaFile const& f, const char* chr, RefVar& rv, int64_t pos_max = std::numeric_limits<int64_t>::max());

/**
 * @brief List functions making sure variants aren't pushed past each other.
 *
 * variants in the input list must be sorted by starting position
 *
 * Functions templated so they work on things derived from RefVar
 *
 */
static inline void leftShift(common::FastaFile const& f, const char* chr, std::list<RefVar>& rv, int64_t pos_min = -1)
{
    int64_t last_start = -1;
    for (RefVar& r : rv)
    {
        if (last_start > 0 && r.start < last_start)
        {
            error("Variants out of order at %s:%i", chr, last_start);
        }
        leftShift(f, chr, (RefVar&)r, pos_min);
        pos_min = r.end + 1;
        last_start = r.start;
    }
}

static inline void rightShift(
    common::FastaFile const& f, const char* chr, std::list<RefVar>& rv,
    int64_t pos_max = std::numeric_limits<int64_t>::max())
{
    for (auto it = rv.rbegin(); it != rv.rend(); ++it)
    {
        RefVar const& r = *it;
        rightShift(f, chr, (RefVar&)r, pos_max);
        pos_max = r.start - 1;
    }
}

/**
 * @brief Quickly decompose a RefVar into primitive variants (subst / ins / del)
 *
 * This function will output matches / mismatches first, and then
 * finish with the insertion / deletion part. The alignment is assumed to be
 * like this for a complex insertion:
 *
 * (reflen < altlen)
 * REF:  XXXYY-----YYZZZ
 *       |||||iiiii
 * ALT:  AAAAABBBBB
 *
 * Or like this for a complex deletion
 *
 * (reflen > altlen)
 * REF:  XXXYY--ZZZ
 *       |||||dd
 * ALT:  AAAAA--
 *
 * To realign variants that do not follow this scheme, see realignRefvar in Alignment.hh
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param vars the primitive records
 */
extern int
toPrimitives(common::FastaFile const& f, const char* chr, RefVar const& rv, std::list<variant::RefVar>& vars);

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param snps the number of snps
 * @param ins the number of insertions
 * @param dels the number of deletions
 * @param homref the number of calls with no variation
 */
void countRefVarPrimitives(
    common::FastaFile const& f, const char* chr, variant::RefVar const& rv, size_t& snps, size_t& ins, size_t& dels,
    size_t& homref, size_t& transitions, size_t& transversions);

/**
 * Convert a list of RefVar records to allele strings
 */
extern void
toAlleles(common::FastaFile const& f, const char* chr, std::vector<RefVar> const& in, std::vector<std::string>& out);

/**
 * Convert a CIGAR string to a list of variants
 *
 * @param ref reference sequence string
 * @param alt alt sequence string
 * @param cigar cigar for aligning ref to alt
 * @param ref_left return the number of reference characters left after processing cigar ops
 * @param alt_left return the number of alt characters left after processing cigar ops
 * @param ref_matches when set to true, this will return ref-matching records with "." as the alt
 */
extern std::list<RefVar> cigarToRefVar(
    std::string const& ref, std::string const& alt, std::string const& cigar, int& ref_left, int& alt_left,
    bool ref_matches = false);

/**
 * Helper to aggregate reference variants base by base
 *
 * Refpos must be > vars.back().end ; variants with different flags will not be combined
 *
 */
void appendToVarList(int64_t refpos, char refchr, char altchr, std::list<variant::RefVar>& vars, int flags = 0);

/**
 * @brief returns true if rv is representing the reference
 *
 * TODO change signature to move chr before rv
 *
 * @param ref refererence sequence fasta file
 * @param rv  an allele
 * @param chr the chromosome
 */
static inline bool isReference(const common::FastaFile& ref, const variant::RefVar& rv, const char* chr)
{
    return (rv.alt == ref.query(chr, rv.start, rv.end));
}

/**
 * Get the length of the shortest ALT allele for a vector of refVars
 * @param rv
 * @return
 */
static inline int64_t getMinAltLen(std::vector<variant::RefVar>& rv)
{
    if (rv.empty())
    {
        return 0;
    }
    size_t min_altlen = (int64_t)rv[0].alt.size();
    for (size_t i = 1; i < rv.size(); i++)
    {
        min_altlen = min_altlen < rv[i].alt.size() ? min_altlen : rv[i].alt.size();
    }
    return ((int64_t)min_altlen);
}

/**
 * A simple debugging function for printing a RefVar
 */
static inline void printRefVar(common::FastaFile const& f, const char* chr, std::vector<variant::RefVar> const& in)
{

    std::cerr << "[";
    for (size_t i = 0; i < in.size(); i++)
    {
        if (i > 0)
        {
            std::cerr << ",";
        }
        std::vector<variant::RefVar> v;
        v.push_back(in[i]);
        std::vector<std::string> tmp;
        variant::toAlleles(f, chr, v, tmp);
        std::cerr << in[i].start << ":" << tmp[0] << "/" << tmp[1];
    }
    std::cerr << "]" << std::endl;
}

/**
 * A simple debugging function for printing a RefVar
 */
static inline void printRefVar(common::FastaFile const& f, const char* chr, variant::RefVar const& in)
{
    std::vector<variant::RefVar> dummy;
    dummy.push_back(in);
    printRefVar(f, chr, dummy);
}

} // namespace variant
