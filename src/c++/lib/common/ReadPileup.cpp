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

#include "common/ReadPileup.hh"

#include <limits>
#include <list>
#include <vector>

#include "common/Error.hh"

namespace common
{

struct ReadPileupOrderByPos
{
    bool operator()(ReadPileupInfo const* lhs, ReadPileupInfo const* rhs) const { return lhs->pos() < rhs->pos(); }
};

struct SinglePosition : public ReadPileupInfo
{
    explicit SinglePosition(int64_t pos__)
        : pos_(pos__)
    {
    }
    int64_t pos_ = 0;

    int64_t pos() const override { return pos_; }
    int64_t end() const override { return pos_; }
};

struct ReadPileup::ReadPileupImpl
{
    std::list<std::unique_ptr<ReadPileupInfo>> read_list;
    std::vector<ReadPileupInfo*> reads;
    int64_t max_distance = 0l;
};

ReadPileup::ReadPileup()
    : _impl(new ReadPileupImpl())
{
}

ReadPileup::ReadPileup(ReadPileup&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

ReadPileup::~ReadPileup() = default;

ReadPileup& ReadPileup::operator=(ReadPileup&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

/**
 * Add a read to the pileup
 * @param readinfo Read info object
 */
void ReadPileup::addRead(std::unique_ptr<ReadPileupInfo> readinfo)
{
    if (!_impl->reads.empty())
    {
        assert(readinfo->pos() >= _impl->reads.back()->pos());
    }
    _impl->max_distance = std::max(_impl->max_distance, readinfo->end() - readinfo->pos() + 2);
    _impl->reads.emplace_back(readinfo.get());
    _impl->read_list.emplace_back(std::move(readinfo));
}

/**
 * Pile up reads at pos
 * @param pos position to calculuate pileup for
 * @param processor function that is called successively for every read that overlaps pos
 */
void ReadPileup::pileup(int64_t pos, std::function<void(ReadPileupInfo*)>&& processor) const
{
    SinglePosition location{ pos };
    ReadPileupOrderByPos order;
    auto last
        = std::upper_bound(_impl->reads.begin(), _impl->reads.end(), static_cast<ReadPileupInfo*>(&location), order);
    if (last == _impl->reads.end() && !_impl->reads.empty())
    {
        last = std::prev(last);
    }
    while (last != _impl->reads.end())
    {
        if ((*last)->pos() <= pos && (*last)->end() >= pos)
        {
            processor((*last));
        }
        else if ((*last)->end() < pos - _impl->max_distance)
        {
            break;
        }
        if (last != _impl->reads.begin())
        {
            --last;
        }
        else
        {
            break;
        }
    }
}

/**
 * Clear out reads that end before pos
 */
void ReadPileup::flush(int64_t pos)
{
    size_t deleted = 0;
    while (!_impl->read_list.empty() && _impl->read_list.front()->end() < pos)
    {
        _impl->read_list.pop_front();
        ++deleted;
    }
    if (deleted > 0)
    {
        _impl->reads.resize(_impl->read_list.size());
        size_t index = 0;
        for (auto& read : _impl->read_list)
        {
            _impl->reads[index++] = read.get();
        }
    }
}
}
