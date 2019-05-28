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
