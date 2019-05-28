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
