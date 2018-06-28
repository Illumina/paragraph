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
 * \brief Error handling and logging helper
 *
 * \file Error.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Error.hh"

#include <cstdlib>

std::string DEFAULT_LOGGER_NAME = "grmpy";

/**
 * Initialize a logger.
 * @param name   logger name
 * @param async  log asynchronously (use true when there are many messages)
 * @param level  minimum level to log at
 */
void initLogging(const char* name, const char* filename, bool async, const char* level)
{
    if (async)
    {
        auto flush_interval = std::chrono::seconds(10);
        spdlog::set_async_mode(8192, spdlog::async_overflow_policy::block_retry, nullptr, flush_interval);
    }
    else
    {
        spdlog::set_sync_mode();
    }

    if (std::string(level) == "trace")
    {
        spdlog::set_level(spdlog::level::trace);
    }
    else if (std::string(level) == "debug")
    {
        spdlog::set_level(spdlog::level::debug);
    }
    else if (std::string(level) == "info")
    {
        spdlog::set_level(spdlog::level::info);
    }
    else if (std::string(level) == "warning")
    {
        spdlog::set_level(spdlog::level::warn);
    }
    else if (std::string(level) == "error")
    {
        spdlog::set_level(spdlog::level::err);
    }
    else
    {
        error("Unknown log level: %s", level);
    }

    // LOG creates + registers the logger, volatile so
    // it doesn't get optimized out
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        if (filename == nullptr || strlen(filename) == 0)
        {
            logger = spdlog::stderr_color_mt(name);
        }
        else
        {
            logger = spdlog::basic_logger_mt(name, filename);
        }
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        error("Log initialization failed: %s", ex.what());
    }
    DEFAULT_LOGGER_NAME = name;

    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%t] [%l] %v");
    // necessary for MSVC, also seems to avoid a memory leak
    std::atexit(spdlog::drop_all);
}
