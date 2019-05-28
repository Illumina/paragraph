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

#include "StringUtil.hh"

#include "Error.hh"

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
    boost::filesystem::path p(file_or_value);
    if (!boost::filesystem::is_regular_file(file_or_value))
    {
        std::istringstream input(file_or_value);
        input >> result;
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
            input >> result;
        }
        else
        {
            std::ifstream input(file_or_value);
            input >> result;
        }
    }
    return result;
}

/**
 * Write a Json::Value object to string
 * indent = false also suppresses line breaks
 */
static inline std::string writeJson(Json::Value const& val, bool indent = true, int float_precision = 5)
{
    Json::StreamWriterBuilder wBuilder;
    wBuilder["precision"] = float_precision;
    wBuilder["useSpecialFloats"] = false;
    if (!indent)
    {
        wBuilder["indentation"] = "";
    }
    return Json::writeString(wBuilder, val);
}
}