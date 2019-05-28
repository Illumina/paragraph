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
 * \brief klib Smith Waterman wrapped. See external/klib.
 *
 * \file Klib.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Klib.hh"
#include "KlibImpl.hh"

namespace common
{

KlibAlignment::KlibAlignment()
    : _impl(new KlibAlignmentImpl())
{
    setParameters(AlignmentParameters());
}

KlibAlignment::KlibAlignment(KlibAlignment&& that)
    : _impl(that._impl)
{
    that._impl = nullptr;
}

KlibAlignment::~KlibAlignment() { delete _impl; }

void KlibAlignment::setParameters(AlignmentParameters const& ap)
{
    for (int i = 0; i < 25; ++i)
    {
        _impl->mat[i] = ap.subs_mat[i];
    }
    _impl->gapo = ap.gapo;
    _impl->gape = ap.gape;
}

void KlibAlignment::getParameters(AlignmentParameters& ap)
{
    for (int i = 0; i < 25; ++i)
    {
        ap.subs_mat[i] = _impl->mat[i];
    }
    ap.gapo = static_cast<uint8_t>(_impl->gapo);
    ap.gape = static_cast<uint8_t>(_impl->gape);
}

/**
 * @brief set target sequence
 */
void KlibAlignment::setRef(const char* seq)
{
    _impl->valid_result = false;
    _impl->reflen = static_cast<int>(strlen(seq));
    _impl->ref = std::shared_ptr<uint8_t>(new uint8_t[_impl->reflen], [](uint8_t* p) { delete[] p; });
    translate(seq, _impl->ref.get(), _impl->reflen);
}

/*
 * @brief set query sequence
 */
void KlibAlignment::setQuery(const char* seq)
{
    if (_impl->qprofile)
    {
        free(_impl->qprofile);
        _impl->qprofile = nullptr;
    }
    _impl->valid_result = false;
    _impl->altlen = static_cast<int>(strlen(seq));
    _impl->alt = std::shared_ptr<uint8_t>(new uint8_t[_impl->altlen], [](uint8_t* p) { delete[] p; });
    translate(seq, _impl->alt.get(), _impl->altlen);
}

/**
 * @brief Get the alignment score
 */
int KlibAlignment::getScore()
{
    if (!_impl->valid_result)
    {
        this->update();
    }
    return _impl->result.score;
}

/**
 * @brief Get a cigar string
 */
void KlibAlignment::getCigar(int& r0, int& r1, int& a0, int& a1, std::string& cig)
{
    if (!_impl->valid_result)
    {
        this->update();
    }
    r0 = _impl->result.tb;
    r1 = _impl->result.te;
    a0 = _impl->result.qb;
    a1 = _impl->result.qe;
    cig = makeCigar(_impl->result.qb, _impl->result.qe, _impl->altlen, _impl->cigar_len, _impl->cigar);
}

/**
 * @brief Get int* cigar string + start + end
 */
void KlibAlignment::getCigar(int& r0, int& r1, int& a0, int& a1, int& n_cigar, uint32_t*& cigar)
{
    if (!_impl->valid_result)
    {
        this->update();
    }
    r0 = _impl->result.tb;
    r1 = _impl->result.te;
    a0 = _impl->result.qb;
    a1 = _impl->result.qe;
    n_cigar = _impl->cigar_len;
    cigar = _impl->cigar;
}

/**
 * @brief Debug dump, optional
 */
void KlibAlignment::update()
{
    _impl->result = ksw_align(
        _impl->altlen, _impl->alt.get(), _impl->reflen, _impl->ref.get(), 5, _impl->mat, _impl->gapo, _impl->gape,
        KSW_XSTART, // add flags here
        &(_impl->qprofile));

    if (_impl->cigar)
    {
        free(_impl->cigar);
        _impl->cigar = nullptr;
        _impl->cigar_len = 0;
    }

    ksw_global(
        _impl->result.qe - _impl->result.qb + 1, _impl->alt.get() + _impl->result.qb,
        _impl->result.te - _impl->result.tb + 1, _impl->ref.get() + _impl->result.tb, 5, _impl->mat, _impl->gapo,
        _impl->gape, _impl->reflen, &_impl->cigar_len, &_impl->cigar);

    _impl->valid_result = true;
}
}
