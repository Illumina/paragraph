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
