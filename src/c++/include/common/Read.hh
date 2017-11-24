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
 *  \brief BAM alignment data structure
 *
 * \file BamAlignment.hh
 * \author Egor Dolzhenko
 * \email edolzhenko@illumina.com
 *
 */

#pragma once

#include <google/protobuf/util/json_util.h>
#include <string>

#include "json/json.h"

#include "read.pb.h"

namespace common
{

// Holds a BAM alignment and decodes them from HTSLib structs.
class Read : public ::reads::Read
{
public:
    Read()
    {
        set_chrom_id(-1);
        set_pos(-1);
        set_mapq(0);
        set_is_mapped(false);
        set_is_first_mate(true);
        set_is_mate_mapped(false);
        set_mate_chrom_id(-1);
        set_mate_pos(-1);
    }

    Read(Read const& rhs) = default;

    Read(std::string const& fragment_id, std::string const& bases, std::string const& quals)
        : Read()
    {
        set_fragment_id(fragment_id);
        set_bases(bases);
        set_quals(quals);
    }

    void setCoreInfo(std::string const& fragment_id, std::string const& bases, std::string const& quals)
    {
        set_fragment_id(fragment_id);
        set_bases(bases);
        set_quals(quals);
    }

    bool is_initialized() const { return !bases().empty(); }

    bool operator==(const Read& other) const
    {
        return fragment_id() == other.fragment_id() && bases() == other.bases() && quals() == other.quals()
            && chrom_id() == other.chrom_id() && pos() == other.pos() && mapq() == other.mapq()
            && is_mapped() == other.is_mapped() && is_first_mate() == other.is_first_mate()
            && is_mate_mapped() == other.is_mate_mapped() && mate_chrom_id() == other.mate_chrom_id()
            && mate_pos() == other.mate_pos();
    }

    Json::Value asJson() const
    {
        Json::Value val;
        std::string str;
        google::protobuf::util::MessageToJsonString(*((Message*)this), &str);
        Json::Reader reader;
        reader.parse(str, val);
        return val;
    }
};

typedef std::unique_ptr<Read> p_Read;
typedef std::vector<p_Read> ReadBuffer;

template <typename read_container> static inline ReadBuffer toReadBuffer(read_container const& reads)
{
    ReadBuffer rb;
    for (Read const& read : reads)
    {
        rb.emplace_back(new Read(read));
    }
    return rb;
}
}
