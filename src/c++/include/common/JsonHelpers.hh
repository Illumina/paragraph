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
 * @brief data structure
 *
 * @author Sai Chen & Peter Krusche & Egor Dolzhenko
 * @email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#pragma once

#include <fstream>
#include <functional>
#include <json/json.h>
#include <memory>
#include <sstream>

#include <boost/filesystem.hpp>
#include <htslib/hts.h>

#include "Error.hh"
#include "StringUtil.hh"

namespace common
{

/**
 * Helper function to return JSON object from a file or from a string
 * @param file_or_value filename or string value
 * @return JSON value
 */
static inline Json::Value getJSON(std::string const& file_or_value)
{
    Json::Value result;
    Json::Reader reader;
    boost::filesystem::path p(file_or_value);
    if (!boost::filesystem::is_regular_file(file_or_value))
    {
        std::istringstream input(file_or_value);
        reader.parse(input, result);
    }
    else
    {
        if (stringutil::endsWith(file_or_value, "gz"))
        {
            auto jsonfile
                = std::unique_ptr<htsFile, std::function<void(htsFile*)>>{ hts_open(file_or_value.c_str(), "rz"),
                                                                           hts_close };
            if (!jsonfile)
            {
                error("Cannot open %s", file_or_value.c_str());
            }

            std::string value;

            int res = 1;
            while (res > 0)
            {
                kstring_t str{ 0, 0, nullptr };
                res = hts_getline(jsonfile.get(), '\n', &str);
                if (res > 0)
                {
                    value += str.s;
                }
                if (str.s != nullptr)
                {
                    free(str.s);
                }
            }
            std::istringstream input(value);
            reader.parse(input, result);
        }
        else
        {
            std::ifstream input(file_or_value);
            reader.parse(input, result);
        }
    }
    return result;
}
}
