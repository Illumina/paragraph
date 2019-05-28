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
 *  \brief String helper functions
 *
 * \file StringUtil.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace common
{
namespace stringutil
{
    /**
     * @brief Split a string along separators
     */
    static inline void
    split(std::string str, std::vector<std::string>& result, const std::string& seps = " ,", bool return_empty = false)
    {
        bool has_last = false;
        while (str.size() > 0)
        {
            size_t pos = str.find_first_of(seps);
            if (pos != std::string::npos)
            {
                if (return_empty || pos > 0)
                {
                    result.push_back(str.substr(0, pos));
                }
                str = str.substr(pos + 1);
                has_last = true;
            }
            else
            {
                result.push_back(str);
                str = "";
                has_last = false;
            }
        }
        if (has_last && return_empty)
        {
            result.push_back("");
        }
    }

    /**
     * @brief Test if string has given suffix
     *
     * @return true if str ends with suffix
     */
    static inline bool endsWith(std::string const& str, std::string const& suffix)
    {
        if (suffix.size() > str.size())
        {
            return false;
        }
        return str.substr(str.size() - suffix.size()) == suffix;
    }

    /**
     * @brief Replace all instances of find with replace in str
     */
    static inline std::string replaceAll(std::string str, std::string const& find, std::string const& replace)
    {
        size_t start_pos = 0;
        while ((start_pos = str.find(find, start_pos)) != std::string::npos)
        {
            str.replace(start_pos, find.length(), replace);
            start_pos += replace.length();
        }
        return str;
    }

    /**
     * @brief Replace all instances of find with replace in str
     */
    static inline void replaceAllInplace(std::string& str, std::string const& find, std::string const& replace)
    {
        size_t start_pos = 0;
        while ((start_pos = str.find(find, start_pos)) != std::string::npos)
        {
            str.replace(start_pos, find.length(), replace);
            start_pos += replace.length();
        }
    }

    /**
     * Format genomic coordinates
     */
    static inline std::string formatPos(const char* chr, int64_t pos = -1, int64_t end = -1)
    {
        std::stringstream ss;

        if (chr)
        {
            ss << chr;
        }
        if (pos >= 0)
        {
            ss << ":";
            ss << (pos + 1);
            if (end >= 0)
            {
                ss << "-" << (end + 1);
            }
        }

        return ss.str();
    }

    /**
     * Format genomic coordinates
     */
    static inline std::string formatPos(std::string const& chr, int64_t pos = -1, int64_t end = -1)
    {
        return formatPos(chr.c_str(), pos, end);
    }

    /**
     * @brief Parse coordinates
     *
     * @param input input string, e.g. "chr1:1,000-2000"
     * @param chr the contig name, e.g. "chr1"
     * @param start the start, e.g. 999
     * @param end the end, e.g. 1999
     *
     * Returned coordinates are 0-based, input coordinates are 1-based.
     *
     */
    static inline void parsePos(std::string input, std::string& chr, int64_t& start, int64_t& end)
    {
        std::vector<std::string> spl;
        split(input, spl, " :-");

        if (spl.size() >= 1)
        {
            chr = spl[0];
        }
        if (spl.size() >= 2)
        {
            start = std::stoll(replaceAll(spl[1], ",", "")) - 1;
        }
        if (spl.size() >= 3)
        {
            end = std::stoll(replaceAll(spl[2], ",", "")) - 1;
        }
    }

    /**
     * upper-case a string
     */
    static inline void toUpper(std::string& str) { std::transform(str.begin(), str.end(), str.begin(), ::toupper); }

    /**
     * Returns true if query is a substring of str.
     */
    static inline bool occurs(const std::string& query, const std::string& str)
    {
        return str.find(query) != std::string::npos;
    }

    /**
     * Convert double to string with a fixed precision
     */
    static inline std::string doubleToStr(double x, int precision, bool fixed = false)
    {
        std::stringstream s;
        s.precision(5);
        if (fixed)
        {
            s << std::fixed;
        }
        s << x;
        return s.str();
    }
}
}
