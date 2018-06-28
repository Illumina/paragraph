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
 * Generic program command line parsing and startup behavior
 *
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>

#include "common/Error.hh"

namespace common
{

namespace bpo = boost::program_options;

/**
 ** Encapsulation of the pocessing of the command line options.
 **
 ** TODO: add config file and environment options
 **/
class Options : boost::noncopyable
{
public:
    enum Action
    {
        RUN,
        HELP,
        VERSION,
        ABORT
    };
    Options();
    virtual ~Options() {}
    Action parse(const char* moduleName, int argc, const char* const argv[]);
    std::string usage() const;
    const std::string& logLevel() const { return logLevel_; }
    const std::string& logFile() const { return logFile_; }
    bool logAsync() const { return logAsync_; }

protected:
    bpo::options_description namedOptions_;
    bpo::options_description unnamedOptions_;
    bpo::positional_options_description positionalOptions_;

    typedef boost::shared_ptr<boost::program_options::option_description> OptionDescriptionPtr;
    typedef std::vector<OptionDescriptionPtr> OptionDescriptionPtrs;
    std::string helpDefaults(const OptionDescriptionPtrs& options) const;
    std::string help(const OptionDescriptionPtrs& options, const bool markdown) const;

private:
    virtual std::string usagePrefix() const = 0;
    virtual std::string usageSuffix() const { return ""; }
    virtual void postProcess(bpo::variables_map&) {}

    static const unsigned MARKDOWN_LINE_LENGTH = 120;
    // "parse" will store the state in vm_ so that "usage" can access the details of parsed command line
    bpo::variables_map vm_;

    std::string logLevel_;
    std::string logFile_;
    bool logAsync_;
};

/**
 ** Unified behavior of all programs.
 **/
template <class O> void run(void (*callback)(const O&), const char* moduleName, int argc, const char* argv[])
{
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        O options;
        const typename O::Action action = options.parse(moduleName, argc, argv);

        logger = LOG();

        if (O::RUN == action)
        {
            callback(options);
        }
        else if (O::HELP == action)
        {
            std::cout << options.usage() << std::endl;
        }
        else if (O::VERSION == action)
        {
            std::cout << "TODO: have compile-time version constant" << std::endl;
        }
        else
        {
            //            std::clog << options.usage() << std::endl;
            exit(1);
        }
    }
    catch (const std::exception& e)
    {
        if (logger)
        {
            logger->critical(e.what());
        }
        else
        {
            std::cerr << e.what() << std::endl;
        }
        exit(1);
    }
}

} // namespace common
