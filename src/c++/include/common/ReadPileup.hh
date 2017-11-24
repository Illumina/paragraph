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

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>

namespace common
{

class ReadPileupInfo
{
public:
    virtual ~ReadPileupInfo() = default;
    virtual int64_t pos() const = 0;
    virtual int64_t end() const = 0;
};

struct SimpleRead : public ReadPileupInfo
{
    explicit SimpleRead(int64_t pos__, int64_t len = 150)
        : pos_(pos__)
        , end_(pos__ + len - 1)
    {
    }
    int64_t pos_;
    int64_t end_;
    int64_t pos() const override { return pos_; }
    int64_t end() const override { return end_; }

    static std::unique_ptr<ReadPileupInfo> make(int64_t pos, int64_t len = 150)
    {
        std::unique_ptr<ReadPileupInfo> result{ new SimpleRead(pos, len) };
        return result;
    }
};

/**
 * Class for storing read segments
 */
class ReadPileup
{
public:
    ReadPileup();
    ReadPileup(ReadPileup const&) = delete;
    ReadPileup(ReadPileup&&) noexcept;
    virtual ~ReadPileup();
    ReadPileup& operator=(ReadPileup const&) = delete;
    ReadPileup& operator=(ReadPileup&&) noexcept;

    /**
     * Add a read to the pileup
     * @param readinfo Read info object
     */
    void addRead(std::unique_ptr<ReadPileupInfo> readinfo);

    /**
     * Pile up reads at pos
     * @param pos position to calculuate pileup for
     * @param processor function that is called successively for every read that overlaps pos
     */
    void pileup(int64_t pos, std::function<void(ReadPileupInfo*)>&& processor) const;

    /**
     * Clear out reads that end before pos
     */
    void flush(int64_t pos);

private:
    struct ReadPileupImpl;
    std::unique_ptr<ReadPileupImpl> _impl;
};
}
