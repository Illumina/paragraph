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
 * \file Error.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

// include this here before we redefine error below
#include "spdlog/spdlog.h"

#ifndef ERRMSG_MAX
#define ERRMSG_MAX 2048
#endif

#ifndef CYGWIN

#ifdef __GNUC__

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>

#endif

inline void error_stop(const char* f, int l, const char* format, ...)
{

#ifdef _DEBUG
    fprintf(stderr, "\n[DEBUG] Error at %s:%i\n", f, l);
    { // extra Debug print when errors occur
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }

// in GNUC, do backtrace
#ifdef __GNUC__
    void* array[20];

    size_t size;
    size = (size_t)backtrace(array, 20);
    backtrace_symbols_fd(array, size, 2);
#endif
    fprintf(stderr, "\n[DEBUG END]\n");
#endif

    char errmsg[ERRMSG_MAX];
    va_list ap;
    va_start(ap, format);
    vsnprintf(errmsg, ERRMSG_MAX, format, ap);
    va_end(ap);
    throw std::runtime_error(errmsg);
}

#else // CYGWIN / windows version

#include <cstdio>
#include <cstdlib>
#include <imagehlp.h>
#include <inttypes.h>
#include <stdexcept>
#include <windows.h>

static inline int addr2line(void const* const addr)
{
    char addr2line_cmd[512] = { 0 };
    HMODULE hModule = GetModuleHandleW(NULL);
    WCHAR program_name[MAX_PATH];
    GetModuleFileNameW(hModule, program_name, MAX_PATH);

    /* have addr2line map the address to the relent line in the code */
    sprintf(addr2line_cmd, "addr2line -f -p -e `cygpath -u '%.256ls'` %p", program_name, addr);

    /* This will print a nicely formatted string specifying the
       function and source line of the address */
    return system(addr2line_cmd);
}

inline void error_stop(const char* f, int l, const char* format, ...)
{
#ifdef _DEBUG
    fprintf(stderr, "\n[DEBUG] Error at %s:%i\n", f, l);
    { // extra Debug print when errors occur
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }

    const int kmaxcallers = 62;

    void* callers[kmaxcallers];
    int count = CaptureStackBackTrace(0, kmaxcallers, callers, NULL);

    for (int i = 0; i < count; i++)
    {
        addr2line(callers[i]);
    }

    fprintf(stderr, "Call Stack:\n");
    for (int i = 0; i < count; i++)
    {
        printf("frame: %2d | 0x%" PRIX64 "\n", i, (uint64_t)callers[i]);
    }

    fprintf(stderr, "\n[DEBUG END]\n");
#endif

    char errmsg[ERRMSG_MAX];
    va_list ap;
    va_start(ap, format);
    vsnprintf(errmsg, ERRMSG_MAX, format, ap);
    va_end(ap);
    throw std::runtime_error(errmsg);
}

#define error(...)                                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        error_stop(__FILE__, __LINE__, __VA_ARGS__);                                                                   \
    } while (0)

#endif

// DEBUG error message
#ifdef error
#undef error
#endif

#define error(...)                                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        error_stop(__FILE__, __LINE__, __VA_ARGS__);                                                                   \
    } while (false)

#include <boost/filesystem.hpp>

/**
 * Error if given file doesn't exist
 * @param path_str
 */
static inline void assertFileExists(const std::string& path_str)
{
    // allow S3 paths.
    if (path_str.substr(0, 5) == "s3://")
    {
        return;
    }
    boost::filesystem::path p(path_str);
    if (!boost::filesystem::is_regular_file(p))
    {
        error("ERROR: %s does not exist", path_str.c_str());
    }
}

/**
 * Error if one of the files doesn't exist
 * @param path_str
 */
template <typename It> void assertFilesExist(It begin, It end) { std::for_each(begin, end, &assertFileExists); }

/**
 * Error if one or more paths have matching file names
 * @param path_str
 */
template <typename It> void assertFileNamesUnique(It begin, It end)
{
    std::vector<boost::filesystem::path> filePaths(begin, end);
    std::sort(
        filePaths.begin(), filePaths.end(),
        [](const boost::filesystem::path& left, const boost::filesystem::path& right) {
            return left.filename() < right.filename();
        });

    const auto dupe = std::adjacent_find(
        filePaths.begin(), filePaths.end(),
        [](const boost::filesystem::path& left, const boost::filesystem::path& right) {
            return left.filename() == right.filename();
        });
    if (filePaths.end() != dupe)
    {
        error(
            "ERROR: file names conflict between these paths: '%s' and '%s'", dupe->string().c_str(),
            (dupe + 1)->string().c_str());
    }
}

/**
 * Initialize a logger.
 * @param name   logger name
 * @param filename output file name (empty string -> use stderr)
 * @param async  log asynchronously (use true when there are many messages)
 * @param level  minimum level to log at
 */
void initLogging(const char* name = "grmpy", const char* filename = "", bool async = true, const char* level = "info");

extern std::string DEFAULT_LOGGER_NAME;
/**
 * Get / initialize the logger
 * @param name   name of the logger
 * @return the spdlog logger
 */
static inline std::shared_ptr<spdlog::logger> LOG(const char* name = nullptr)
{
    if (name == nullptr)
    {
        name = DEFAULT_LOGGER_NAME.c_str();
    }
    auto log = spdlog::get(name);
    if (!log)
    {
        log = spdlog::stderr_color_mt(name);
    }
    return log;
}

    /* this needs to come at the end after all #include statements s.t. it doesn't get
       reset by any header that includes cassert */

#ifdef assert
#undef assert
#endif

/** redefine assert to use error and make backtraces */
#ifdef assert
#undef assert
#endif
#define assert(x)                                                                                                      \
    do                                                                                                                 \
    {                                                                                                                  \
        if (!(x))                                                                                                      \
            error("Assertion failed: %s", #x);                                                                         \
    } while (0)
