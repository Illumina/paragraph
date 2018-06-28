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
 * \brief Graph caller parameters
 *
 * \file Parameters.cpp
 * \author Peter Krusche & Egor Dolzhenko
 * \email pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "paragraph/Parameters.hh"
#include "common/Error.hh"
#include <fstream>

#include "common/StringUtil.hh"
#include "json/json.h"

namespace paragraph
{

void Parameters::load(
    const std::string& graph_path, const std::string& reference_path, const std::string& override_target_regions)
{
    reference_path_ = reference_path;

    Json::Value root;
    std::ifstream graph_desc(graph_path);
    graph_desc >> root;

    // compatibility with graph key
    if (root.isMember("graph"))
    {
        for (auto& key_name : root["graph"].getMemberNames())
        {
            root[key_name] = root["graph"][key_name];
        }
        root.removeMember("graph");
    }
    description_ = root;

    if (!override_target_regions.empty())
    {
        std::vector<std::string> regions;
        common::stringutil::split(override_target_regions, regions);
        target_regions_.insert(target_regions_.end(), regions.begin(), regions.end());
    }
    else
    {
        // add target regions
        if (!root.isMember("target_regions") || root["target_regions"].type() != Json::ValueType::arrayValue)
        {
            error("Graph description is missing \"target_regions\" key.");
        }
        for (auto& r : root["target_regions"])
        {
            auto s = r.asString();
            target_regions_.emplace_back(s);
        }
    }
    if (root.isMember("max_reads"))
    {
        max_reads_ = root["max_reads"].asUInt64();
    }

    for (auto& node : description_["nodes"])
    {
        if (node.isMember("sequence") && node["sequence"].asString().size() > longest_alt_insertion_)
        {
            longest_alt_insertion_ = node["sequence"].asString().size();
        }
    }
}
}
