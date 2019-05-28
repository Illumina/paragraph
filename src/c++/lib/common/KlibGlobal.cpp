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
 *  \brief Implementation: global alignment via klib
 *
 *
 * \file KlibGlobal.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "KlibGlobal.hh"
#include "KlibImpl.hh"

#include "common/Error.hh"

// see ksw.c
#define MINUS_INF -0x40000000

namespace common
{

void KlibGlobalAlignment::update()
{
    _impl->result.qb = 0;
    _impl->result.qe = _impl->altlen - 1;
    _impl->result.tb = 0;
    _impl->result.te = _impl->reflen - 1;

    if (_impl->cigar)
    {
        free(_impl->cigar);
        _impl->cigar = NULL;
        _impl->cigar_len = 0;
    }

    _impl->result.score = ksw_global(
        _impl->result.qe - _impl->result.qb + 1, _impl->alt.get() + _impl->result.qb,
        _impl->result.te - _impl->result.tb + 1, _impl->ref.get() + _impl->result.tb, 5, _impl->mat, _impl->gapo,
        _impl->gape, std::max(_impl->reflen, _impl->altlen), &_impl->cigar_len, &_impl->cigar);

    if (_impl->result.score <= MINUS_INF)
    {
        error("Failed to globally align.");
    }

    _impl->valid_result = true;
}
}
