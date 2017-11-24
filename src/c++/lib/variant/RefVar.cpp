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
 * \brief Reference Variation helper implementation
 *
 *
 * \file RefVar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "variant/RefVar.hh"
#include "common/Error.hh"
#include "common/Genetics.hh"

#include <algorithm>
#include <boost/concept_check.hpp>
#include <cassert>
#include <limits>
//#include <etip.h>

using namespace common;

// #define DEBUG_REFVAR

namespace variant
{

void trimLeft(std::string const& ref, RefVar& rv, bool refpadding, int64_t rel_start)
{
    size_t ref_min = refpadding ? 1 : 0;

    while (ref.size() - rel_start > ref_min && rv.alt.size() - rel_start > ref_min
           && ref[rel_start] == rv.alt[rel_start])
    {
        rel_start++;
        rv.start++;
    }
    if (rel_start > 0)
    {
        rv.alt = rv.alt.substr((unsigned long)rel_start);
    }
}

void trimRight(std::string const& ref, RefVar& rv, bool refpadding, int64_t rel_start)
{
    // trim right
    int64_t reflen = rv.end - rv.start + 1;
    int64_t altlen = (int64_t)rv.alt.size();
    int64_t min_len = refpadding ? 1 : 0;

    if (reflen <= min_len || altlen <= min_len)
    {
        return;
    }

#ifdef DEBUG_REFVAR
    std::cerr << rv << " -- ref sq = " << ref << "\n";
#endif

    while (reflen > min_len && altlen > min_len && ref[rel_start + reflen - 1] == rv.alt[altlen - 1])
    {
        altlen--;
        reflen--;
    }
    rv.end = rv.start + reflen - 1;
    if (altlen > 0)
    {
        rv.alt = rv.alt.substr(0, (unsigned long)altlen);
    }
    else
    {
        rv.alt = "";
    }
}

void trimLeft(FastaFile const& f, const char* chr, RefVar& rv, bool refpadding)
{
    // trim left
    std::string ref = f.query(chr, rv.start, rv.end);
    int64_t rel_start = 0;
    trimLeft(ref, rv, refpadding, rel_start);
}

void trimRight(FastaFile const& f, const char* chr, RefVar& rv, bool refpadding)
{
    std::string ref = f.query(chr, rv.start, rv.end);
    trimRight(ref, rv, refpadding, 0);
}

void leftShift(std::string const& ref, RefVar& rv, int64_t pos_min)
{
    pos_min = std::max(pos_min, (int64_t)0);

    trimLeft(ref, rv);
    trimRight(ref, rv);

    int64_t reflen = rv.end - rv.start + 1;

    if (reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    if (reflen >= 0 && reflen == (signed)rv.alt.size())
    {
        const std::string ref_allele = ref.substr((unsigned long)rv.start, (unsigned long)(rv.end - rv.start + 1));
        if (ref_allele == rv.alt)
        {
            return;
        }
    }

    bool done = false;
    while (!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if (rv.start <= pos_min)
        {
            break;
        }

        if (rv.start < 1 || ref.size() == 0 || ((signed)ref.size()) < rv.start + reflen
            || ref[rv.start - 1] == 'N') // don't shift over Ns
        {
            break;
        }

        // right trim.
        if (reflen > 0 && rv.alt.size() > 0 && ref[rv.start + reflen - 1] == rv.alt.back())
        {
            reflen--;
            rv.end--;
            rv.alt = rv.alt.substr(0, rv.alt.size() - 1);
            done = false;
        }

        if (reflen == 0 || rv.alt.size() == 0)
        {
            rv.start--;
            rv.alt = ref.substr((unsigned long)rv.start, 1) + rv.alt;
            done = false;
        }
    }
    trimLeft(ref, rv);
    trimRight(ref, rv);
}

void rightShift(std::string const& ref, RefVar& rv, int64_t pos_max)
{
    trimLeft(ref, rv);
    trimRight(ref, rv);

    int64_t reflen = rv.end - rv.start + 1;
    if (reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    if (reflen >= 0 && reflen == (signed)rv.alt.size())
    {
        const std::string ref_allele = ref.substr((unsigned long)rv.start, (unsigned long)(rv.end - rv.start + 1));
        if (ref_allele == rv.alt)
        {
            return;
        }
    }

    bool done = false;
    while (!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if (rv.end >= pos_max)
        {
            break;
        }

        if (ref.size() == 0 || ((signed)ref.size()) <= rv.start + reflen || ref[rv.start + reflen] == 'N')
        {
            break;
        }

        // left trim.
        if (reflen > 0 && rv.alt.size() > 0 && ref[rv.start] == rv.alt[0])
        {
            reflen--;
            rv.start++;
            rv.alt = rv.alt.substr(1, rv.alt.size());
            done = false;
        }

        // right-extend
        if (reflen == 0 || rv.alt.size() == 0)
        {
            int64_t refnext = rv.start + reflen;
            rv.end++;
            rv.alt = rv.alt + ref.substr((unsigned long)refnext, 1);
            done = false;
        }
    }

    trimLeft(ref, rv);
    trimRight(ref, rv);
}

void leftShift(FastaFile const& f, const char* chr, RefVar& rv, int64_t pos_min)
{
    int64_t rstart = -1, rend = -1, reflen;

    // vaguely like
    // http://genome.sph.umich.edu/wiki/File:Variant_normalization_algorithm.png
    bool done = false;
    std::string ref;

    pos_min = std::max(pos_min, (int64_t)0);

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);

    reflen = rv.end - rv.start + 1;
    // check for all ref match (HAP-64)
    if (reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    if (reflen >= 0 && reflen == (signed)rv.alt.size())
    {
        std::string ref_allele = f.query(chr, rv.start, rv.end);
        if (ref_allele == rv.alt)
        {
            return;
        }
    }

    while (!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if (rstart < 0 || rv.start <= rstart || rend < 0 || rv.end > rend)
        {
            rstart = rv.start - 20;
            if (reflen <= 0)
            {
                rend = rv.start;
            }
            else
            {
                rend = rv.end;
            }

            if (rstart < 0)
            {
                rstart = 0;
            }
            ref = f.query(chr, rstart, rend);
        }
        if (rv.start <= pos_min)
        {
            break;
        }

        int64_t rel_start = rv.start - rstart;

        if (rel_start < 1 || ref.size() == 0 || ((signed)ref.size()) < rel_start + reflen
            || ref[rel_start - 1] == 'N') // don't shift over Ns
        {
            done = true;
            break;
        }

        // right trim.
        if (reflen > 0 && rv.alt.size() > 0 && ref[rel_start + reflen - 1] == rv.alt.back())
        {
            reflen--;
            rv.end--;
            rv.alt = rv.alt.substr(0, rv.alt.size() - 1);
            done = false;
        }

        if (reflen == 0 || rv.alt.size() == 0)
        {
            rv.start--;
            rv.alt = ref.substr((unsigned long)(rel_start - 1), 1) + rv.alt;
            done = false;
        }
    }
    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);
}

void rightShift(FastaFile const& f, const char* chr, RefVar& rv, int64_t pos_max)
{
    int64_t rstart = -1, rend = -1, reflen;

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);

    reflen = rv.end - rv.start + 1;
    // check for all ref match (HAP-64)
    if (reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    if (reflen >= 0 && reflen == (signed)rv.alt.size())
    {
        std::string ref = f.query(chr, rv.start, rv.end);
        if (ref == rv.alt)
        {
            return;
        }
    }

    // adapted from
    // http://genome.sph.umich.edu/wiki/File:Variant_normalization_algorithm.png
    bool done = false;
    std::string ref;
    while (!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if (rstart < 0 || rv.start < rstart || rend < 0 || rv.end >= rend)
        {
            rstart = rv.start;
            if (reflen <= 0)
            {
                rend = rv.start + 20;
            }
            else
            {
                rend = rv.end + 20;
            }

            if (rstart < 0)
            {
                rstart = 0;
            }
            ref = f.query(chr, rstart, rend);
        }
        if (rv.end >= pos_max)
        {
            break;
        }

        int64_t rel_start = rv.start - rstart;

        if (ref.size() == 0 || ((signed)ref.size()) <= rel_start + reflen || ref[rel_start + reflen] == 'N')
        {
            done = true;
            break;
        }

        // left trim.
        if (reflen > 0 && rv.alt.size() > 0 && ref[rel_start] == rv.alt[0])
        {
            reflen--;
            rel_start++;
            rv.start++;
            rv.alt = rv.alt.substr(1, rv.alt.size());
            done = false;
        }

        // right-extend
        if (reflen == 0 || rv.alt.size() == 0)
        {
            int64_t refnext = rel_start + reflen;
            rv.end++;
            rv.alt = rv.alt + ref.substr(refnext, 1);
            done = false;
        }
    }

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);
}

/**
 * Convert a list of RefVar records to allele strings
 */
extern void toAlleles(FastaFile const& f, const char* chr, std::vector<RefVar> const& in, std::vector<std::string>& out)
{
    int64_t minpos = std::numeric_limits<int64_t>::max(), maxpos = 0;

    if (in.size() == 0)
    {
        return;
    }

    for (size_t s = 0; s < in.size(); ++s)
    {
        RefVar const& rv = in[s];
        minpos = std::min(minpos, rv.start);
        minpos = std::min(minpos, rv.end);
        // insertions
        if (rv.start > rv.end)
        {
            maxpos = std::max(maxpos, rv.end);
        }
        else
        {
            maxpos = std::max(maxpos, rv.start);
            maxpos = std::max(maxpos, rv.end);
        }
    }
    std::string ref = f.query(chr, minpos, maxpos);
    if ((signed)ref.size() != maxpos - minpos + 1)
    {
        error("Cannot query reference sequence at %s:%i-%i", ref.c_str(), minpos, maxpos);
    }
    out.resize(in.size() + 1);
    out[0] = ref;
    for (size_t i = 0; i < in.size(); ++i)
    {
        // ref allele for this one:
        // ref.substr(in[i].start - minpos, in[i].end - in[i].start + 1);
        int64_t refstart = in[i].start - minpos;
        int64_t reflen = in[i].end - in[i].start + 1;
        if (reflen == (signed)ref.size() && refstart == 0)
        {
            out[i + 1] = in[i].alt;
        }
        else
        {
            // create alt allele by replacing
            out[i + 1] = ref;
            out[i + 1].replace(refstart, reflen, in[i].alt);
        }
        if (out[i + 1].size() == 0)
        {
            out[i + 1] = "<DEL>";
        }
    }
}

/**
 * Make RefVar
 */
void mkRefVar(int64_t refpos, char refchr, char altchr, RefVar& rv)
{
    rv.start = refpos;
    rv.end = refpos;
    if (refchr == '-')
    {
        // alt insertion, refchr != altchr so altchr != '-'
        rv.end--;
        rv.alt = altchr;
    }
    else if (altchr == '-')
    {
        // ref deletion
        rv.alt = "";
    }
    else
    {
        // snp
        rv.alt = altchr;
    }
}

/**
 * Append 1 char to RefVar.
 *
 * Return false if this is not possible because the variants are incompatible
 */
bool appendToRefVar(int64_t refpos, char refchr, char altchr, RefVar& rv)
{
#ifdef DEBUG_REFVAR
    std::cerr << "\t" << refpos << " " << refchr << " " << altchr << " " << rv << "\n";
#endif
    if (refchr == altchr)
    {
        // already checked in function just below actually
        // we test anyway in case we ever want to expose this
        // function
        return false;
    }
    // see what we have currently
    if ((refpos == rv.start - 1 || refpos == rv.start) && rv.end == rv.start - 1
        && rv.alt.size() > 0) // rv is an insertion just before refpos, or a complex insertion ending at refpos
    {
        // we've been passed an insertion and we can extend?
        if (refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            rv.alt += altchr;
        }
        // we've been passed a substitution?
        else if (refchr != '-' && altchr != '-')
        {
            if (rv.start > refpos)
            {
                // first subst after insertion -- fix start and end
                rv.start = refpos;
                rv.end = refpos;
            }
            else
            {
                rv.end++;
            }
            rv.alt += altchr;
        }
        else if (altchr == '-') //  -> && refchr != '-'; would have caught this above otherwise
        {
            // we've been passed a deletion
            if (rv.start > refpos)
            {
                // first deletion after insertion -- fix start and end
                rv.start = refpos;
                rv.end = refpos;
            }
            else
            {
                rv.end++;
            }
        }
        else
        {
            return false;
        }
        return true;
    }

    // here, we append only
    if (refpos != rv.end + 1)
    {
        return false;
    }

    if (rv.end >= rv.start && rv.alt.size() > 0) // block substitution
    {
        // we've been passed a substitution?
        if (refchr != '-' && altchr != '-')
        {
            rv.end++;
            rv.alt += altchr;
        }
        else if (refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            // we've been passed an insertion
            rv.alt += altchr;
        }
        else if (altchr == '-') //  -> && refchr != '-' since refchr != altchr
        {
            // we've been passed a deletion
            rv.end++;
        }
        else
        {
            return false;
        }
    }
    else if (rv.end >= rv.start && rv.alt.size() == 0) // deletion
    {
        // we've been passed a substitution?
        if (refchr != '-' && altchr != '-')
        {
            rv.end++;
            rv.alt += altchr;
        }
        else if (refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            // we've been passed an insertion
            rv.alt += altchr;
        }
        else if (altchr == '-') //  -> && refchr != '-' since refchr != altchr
        {
            // we've been passed a deletion
            rv.end++;
        }
        else
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del)
 *
 * @param rv the RefVar record
 * @param vars the primitive records
 */
int toPrimitives(FastaFile const& f, const char* chr, RefVar const& rv, std::list<variant::RefVar>& vars)
{
    int64_t rstart = rv.start, rend = rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t)rv.alt.size();

    std::string refseq;
    std::string altseq(rv.alt);

    if (reflen > 0)
    {
        refseq = f.query(chr, rstart, rend);
    }

    // from the left, split off SNPs / matches
    size_t pos = 0;
    int count = 0;
    while (reflen > 0 && altlen > 0)
    {
        char r = refseq[pos];
        char a = altseq[pos];
        if (r != a)
        {
            vars.push_back(RefVar{ rstart, rstart, std::string(1, a), rv.flags });
            count++;
        }
        ++rstart;
        ++pos;
        --reflen;
        --altlen;
    }
    // now either reflen == 0 or altlen == 0
    if (reflen > 0)
    {
        // del
        vars.push_back(RefVar{ rstart, rend, "", rv.flags });
        count++;
    }
    else if (altlen > 0)
    {
        // ins
        vars.push_back(RefVar{ rstart, rstart - 1, altseq.substr(pos), rv.flags });
        count++;
    }
    // else nothing left
    return (count);
}

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
    FastaFile const& f, const char* chr, variant::RefVar const& rv, size_t& snps, size_t& ins, size_t& dels,
    size_t& homref, size_t& transitions, size_t& transversions)
{
    int64_t rstart = rv.start, rend = rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t)rv.alt.size();

    std::string refseq;
    std::string altseq(rv.alt);

    if (reflen <= 0)
    {
        if (altlen > 0)
        {
            ins += altlen;
        }
        else
        {
            ++homref;
        }
        return;
    }
    // reflen > 0
    refseq = f.query(chr, rstart, rend);

    // from the left, split off SNPs / matches
    size_t pos = 0;
    bool isValidSnv(false);

    while (reflen > 0 && altlen > 0)
    {
        char r = refseq[pos];
        char a = altseq[pos];
        if (r != a)
        {
            ++snps;
            const bool isTransversion(snvIsTransversion(r, a, isValidSnv));

            if (isValidSnv)
            {
                if (isTransversion)
                {
                    ++transversions;
                }
                else
                {
                    ++transitions;
                }
            }
        }
        else
        {
            ++homref;
        }
        ++rstart;
        ++pos;
        --reflen;
        --altlen;
    }
    // now either reflen == 0 or altlen == 0
    if (reflen > 0)
    {
        // del
        dels += reflen;
    }
    if (altlen > 0)
    {
        // ins
        ins += altlen;
    }
    // else nothing left
}

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
std::list<RefVar> cigarToRefVar(
    std::string const& refseq, std::string const& altseq, std::string const& cigar, int& ref_left, int& alt_left,
    bool ref_matches)
{
    std::list<RefVar> result;

    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";
    rv.flags = 0;

    size_t refpos = 0;
    size_t altpos = 0;
    size_t cigarpos = 0;

    while (cigarpos < cigar.size())
    {
        std::string count_str;
        while (cigarpos < cigar.size() && cigar[cigarpos] >= '0' && cigar[cigarpos] <= '9')
        {
            count_str.push_back(cigar[cigarpos]);
            ++cigarpos;
        }
        auto count = (uint32_t)atoll(count_str.c_str());
        if (cigarpos >= cigar.size())
        {
            break;
        }
        auto op = (uint8_t)cigar[cigarpos++];
        if (count == 0)
        {
            error("Invalid CIGAR string with zero-length operation: %s", cigar.c_str());
        }
        switch (op)
        {
        case 'S':
            refpos += count;
            break;
        case 'M': // 'M'
        case '=': // '='
        case 'X': // 'X'
        {
            int ref_match_count = 0;
            for (uint32_t j = 0; j < count; ++j)
            {
                if (refpos >= refseq.size() || altpos > altseq.size())
                {
                    break;
                }
                if (refseq[refpos] != altseq[altpos])
                {
                    if (ref_match_count != 0)
                    {
                        rv.start = (int64_t)refpos - ref_match_count;
                        rv.end = (int64_t)refpos - 1;
                        rv.flags = (int64_t)altpos - ref_match_count;
                        rv.alt = ".";
                        result.push_back(rv);
                        ref_match_count = 0;
                    }
                    rv.start = (int64_t)refpos;
                    rv.end = (int64_t)refpos;
                    rv.alt = altseq[altpos];
                    rv.flags = (int64_t)altpos; // return alt pos for extracting quality later
                    result.push_back(rv);
                }
                else if (ref_matches)
                {
                    ++ref_match_count;
                }
                ++refpos;
                ++altpos;
            }
            if (ref_match_count != 0)
            {
                rv.start = (int64_t)refpos - ref_match_count;
                rv.end = (int64_t)refpos - 1;
                rv.alt = ".";
                rv.flags = (int64_t)altpos - ref_match_count;
                result.push_back(rv);
            }
            break;
        }
        case 'I': // 'I' -> REF insertion in ALT = ALT deletion
            rv.start = (int64_t)refpos;
            rv.end = (int64_t)(refpos + count - 1);
            rv.alt = "";
            rv.flags = (int64_t)altpos;
            result.push_back(rv);
            // shift the reference position
            refpos += count;
            break;
        case 'D': // 'D' -> REF deletion = ALT insertion;
            // insert before reference pos

            // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
            // this is interpreted as an insertion before pos.
            rv.start = (int64_t)refpos;
            rv.end = (int64_t)(refpos - 1);
            rv.alt = altseq.substr(altpos, count);
            rv.flags = (int64_t)altpos;
            result.push_back(rv);
            altpos += count;
            break;
        default:
            error("Unknown CIGAR operation: %c", op);
        }
        if (refpos >= refseq.size() || altpos > altseq.size())
        {
            break;
        }
    }
    ref_left = (int)(refseq.size() - refpos);
    alt_left = (int)(altseq.size() - altpos);
    return result;
}

/**
 * Helper to aggregate reference variants base by base
 *
 * Refpos must be > vars.back().end
 *
 */
void appendToVarList(int64_t refpos, char refchr, char altchr, std::list<variant::RefVar>& vars, int flags)
{
#ifdef DEBUG_REFVAR
    std::cerr << refpos << " " << refchr << " " << altchr << " nv: " << vars.size() << "\n";
#endif
    if (refchr == altchr)
    {
        return;
    }
    if (vars.empty())
    {
        RefVar rv;
        mkRefVar(refpos, refchr, altchr, rv);
        rv.flags = flags;
        vars.push_back(rv);
    }
    else
    {
        RefVar& xrv = vars.back();
        if (flags != xrv.flags || !appendToRefVar(refpos, refchr, altchr, xrv))
        {
            RefVar rv;
            mkRefVar(refpos, refchr, altchr, rv);
            rv.flags = flags;
            vars.push_back(rv);
        }
    }
}

} // namespace variant
