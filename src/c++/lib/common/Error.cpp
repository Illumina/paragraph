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
