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
 *  \brief Alignment and Factory implementation
 *
 * \file Alignment.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Alignment.hh"

#include <cstdlib>
#include <cstring>
#include <sstream>

#include "KlibGlobal.hh"
#include "common/Klib.hh"

#include "common/Error.hh"

using namespace variant;

namespace common
{

Alignment::~Alignment() {}

/* get the best possible score for comparing two sequences of length len */
int Alignment::bestScore(int len)
{
    AlignmentParameters ap;
    getParameters(ap);
    return len * ap.maxScore();
}

Alignment* makeAlignment(const char* type)
{
    if (strstr(type, "klibg") == type)
    {
        return new KlibGlobalAlignment();
    }
    else if (strstr(type, "klib") == type)
    {
        return new KlibAlignment();
    }
    else
    {
        error("Unknown alignment type '%s'", type);
    }
    return NULL;
}

inline char getCigarOp(const char s, const char r) { return s == r ? 'M' : s == 'N' ? 'N' : r == 'N' ? 'N' : 'X'; }

/**
 * \param  matches will be incremented by the number of matches found
 * \return CIGAR with mismatches clipped on left if first is set
 *         or on right if last is set. If the entire sequence is
 *         clipped, an empty CIGAR is returned
 */
std::string
makeCigarBit(std::string::const_iterator itRef, std::string::const_iterator& itSeq, std::size_t length, int& matches)
{
    std::string ret;

    char lastOp = 0;
    std::size_t lastOpLen = 0;
    for (; length; ++itSeq, ++itRef, --length, ++lastOpLen)
    {
        const char op = getCigarOp(*itRef, *itSeq);
        if (op != lastOp)
        {
            if (lastOpLen)
            {
                ret += std::to_string(lastOpLen) + lastOp;
                if ('M' == lastOp)
                {
                    matches += lastOpLen;
                }
            }
            lastOp = op;
            lastOpLen = 0;
        }
    }
    if (lastOpLen)
    {
        ret += std::to_string(lastOpLen) + lastOp;
        if ('M' == lastOp)
        {
            matches += lastOpLen;
        }
    }

    return ret;
}

/**
 * @brief Format int encoded Cigar string
 *
 * @param qb for padding with "S" : begin
 * @param qe for padding with "S" : end
 * @param altlen for padding with "S" : length of alternate sequence
 * @param ncigar length of cigar
 * @param cigar int* to cigar entries
 *
 */
std::string makeCigar(int qb, int qe, int altlen, int ncigar, uint32_t* cigar)
{
    std::string cig;

    if (ncigar > 0)
    {
        std::ostringstream cigar_string;
        if (qb > 0)
        {
            cigar_string << qb << 'S';
        }

        for (int i = 0; i < ncigar; ++i)
        {
            cigar_string << (cigar[i] >> 4);
            uint8_t op = static_cast<uint8_t>(cigar[i] & 0x000f);
            switch (op)
            {
            case 0:
                cigar_string << 'M';
                break;
            case 1:
                cigar_string << 'I';
                break;
            case 2:
                cigar_string << 'D';
                break;
            }
        }

        int end = altlen - qe - 1;
        if (end > 0)
        {
            cigar_string << end << 'S';
        }

        cig = cigar_string.str();
    }

    return cig;
}

/** make variants from a cigar string */
void getVariantsFromCigar(
    std::string const& ref, std::string const& alt, int r0, int a0, uint32_t* cigar, int ncigar,
    std::list<variant::RefVar>& target)
{

    bool have_rv = false;
    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";

    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = cigar[i] >> 4;
        uint8_t op = cigar[i] & 0x000f;
        switch (op)
        {
        case 0: // 'M'
        case 7: // '='
        case 8: // 'X'
            for (uint32_t j = 0; j < count; ++j)
            {
                // check match / mismatch because ksw doesn't give this to us
                if (ref[refpos] != alt[altpos])
                {
                    // push a previous block which we can't append to
                    if (have_rv && rv.end < refpos - 1)
                    {
                        target.push_back(rv);
                        have_rv = false;
                        rv.start = -1;
                        rv.end = -1;
                        rv.alt = "";
                    }

                    if (rv.start < 0)
                    {
                        rv.start = refpos;
                    }
                    rv.end = refpos;
                    rv.alt += alt[altpos];
                    have_rv = true;
                }
                else if (have_rv)
                {
                    target.push_back(rv);
                    have_rv = false;
                    rv.start = -1;
                    rv.end = -1;
                    rv.alt = "";
                }
                ++refpos;
                ++altpos;
            }
            break;
        case 2: // 'D' -> REF deletion = ALT insertion;
            if (have_rv)
            {
                target.push_back(rv);
                have_rv = false;
                rv.start = -1;
                rv.end = -1;
                rv.alt = "";
            }
            rv.start = refpos;
            rv.end = refpos + count - 1;
            rv.alt = "";
            have_rv = true;

            // shift the reference position
            refpos += count;
            break;
        case 1: // 'I' -> REF insertion in ALT = ALT deletion
            if (have_rv)
            {
                target.push_back(rv);
                have_rv = false;
                rv.start = -1;
                rv.end = -1;
                rv.alt = "";
            }
            // insert before reference pos

            // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
            // this is interpreted as an insertion before pos.
            rv.start = refpos;
            rv.end = refpos - 1;
            rv.alt = alt.substr(altpos, count);

            have_rv = true;
            // shift the reference position up by one
            altpos += count;
            break;
        }
    }

    if (have_rv)
    {
        target.push_back(rv);
    }
}

/** get stats from a cigar string */
void getCigarStats(
    std::string const& ref, std::string const& alt, int r0, int a0, uint32_t* cigar, int ncigar, int& softclipped,
    int& matches, int& mismatches, int& ins, int& del)
{
    softclipped = a0;
    matches = 0;
    mismatches = 0;
    ins = 0;
    del = 0;
    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = cigar[i] >> 4;
        uint8_t op = cigar[i] & 0x000f;
        switch (op)
        {
        case 0: // 'M'
        case 7: // '='
        case 8: // 'X'
            for (uint32_t j = 0; j < count; ++j)
            {
                // check match / mismatch because ksw doesn't give this to us
                if (ref[refpos] != alt[altpos])
                {
                    ++mismatches;
                }
                else
                {
                    ++matches;
                }
                ++refpos;
                ++altpos;
            }
            break;
        case 1: // 'I' -> REF insertion in ALT = ALT deletion
            ins += count;
            // shift the reference position
            altpos += count;
            break;
        case 2: // 'D' -> REF deletion = ALT insertion;
            del += count;
            // shift the reference position up by one
            refpos += count;
            break;
        }
    }

    // TODO we might want to check if we're at the end of alt here.
    softclipped += alt.size() - altpos;
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param in_rv the RefVar record
 * @param aln the alignment interface to use
 * @param vars the primitive records
 */
void realignRefVar(
    FastaFile const& f, const char* chr, RefVar const& in_rv, Alignment* aln, std::list<variant::RefVar>& vars)
{
    int64_t rstart = in_rv.start, rend = in_rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t)in_rv.alt.size();

    if (reflen < 2 || altlen < 2)
    {
        // no complex ref / alt => use fast and simple function
        toPrimitives(f, chr, in_rv, vars);
        return;
    }

    std::string refseq = f.query(chr, rstart, rend);
    std::string altseq = in_rv.alt;

    aln->setRef(refseq.c_str());
    aln->setQuery(altseq.c_str());

    uint32_t* icigar;
    int r0, r1, a0, a1;
    int ncigar = 0;
    aln->getCigar(r0, r1, a0, a1, ncigar, icigar);

    RefVar rv;
    rv.start = -1;
    rv.end = -1;
    rv.alt = "";

    int refpos = r0;
    int altpos = a0;

    for (int i = 0; i < ncigar; ++i)
    {
        uint32_t count = icigar[i] >> 4;
        uint8_t op = icigar[i] & 0x000f;
        switch (op)
        {
        case 0: // 'M'
        case 7: // '='
        case 8: // 'X'
            for (uint32_t j = 0; j < count; ++j)
            {
                // check match / mismatch because ksw doesn't give this to us
                if (refseq[refpos] != altseq[altpos])
                {
                    rv.start = rstart + refpos;
                    rv.end = rstart + refpos;
                    rv.alt = altseq[altpos];
                    vars.push_back(rv);
                }
                ++refpos;
                ++altpos;
            }
            break;
        case 2: // 'D' -> REF deletion = ALT insertion;
            rv.start = rstart + refpos;
            rv.end = rstart + refpos + count - 1;
            rv.alt = "";
            vars.push_back(rv);
            // shift the reference position
            refpos += count;
            break;
        case 1: // 'I' -> REF insertion in ALT = ALT deletion
            // insert before reference pos

            // reference length = end - start + 1 == refpos-1 - refpos + 1 == 0
            // this is interpreted as an insertion before pos.
            rv.start = rstart + refpos;
            rv.end = rstart + refpos - 1;
            rv.alt = altseq.substr(altpos, count);
            vars.push_back(rv);

            altpos += count;
            break;
        }
    }
}
}
