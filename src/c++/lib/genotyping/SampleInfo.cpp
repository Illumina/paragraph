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

#include "genotyping/SampleInfo.hh"
#include "common/StringUtil.hh"

#include <fstream>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

#include "common/Error.hh"
#include "common/JsonHelpers.hh"

using std::list;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace genotyping
{

void SampleInfo::set_autosome_depth(double autosome_depth)
{
    autosome_depth_ = autosome_depth;
    if (depth_sd_ == 0)
    {
        depth_sd_ = sqrt(autosome_depth_ * 5);
    }
}

void SampleInfo::set_sex(std::string sex_string)
{
    std::transform(sex_string.begin(), sex_string.end(), sex_string.begin(), ::tolower);
    if (sex_string.size() > 0 && tolower(sex_string[0]) == 'm')
    {
        sex_ = MALE;
    }
    else if (sex_string.size() > 0 && tolower(sex_string[0]) == 'f')
    {
        sex_ = FEMALE;
    }
    else if (sex_string.size() > 0 && tolower(sex_string[0]) == 'u')
    {
        sex_ = UNKNOWN;
    }
    else
    {
        error("illegal sex string: %s", sex_string.c_str());
    }
}

/**
 * Load manifest file. A manifest contains
 *
 *  * Sample names
 *  * Sample BAM/CRAM file locations
 *  * Depth estimates, or location of idxdepth output
 *  * (optionally) location of alignment JSON
 *
 * @param filename file name of manifest
 * @return list of sample info records
 */
Samples loadManifest(const std::string& filename)
{
    std::ifstream manifest_file(filename.c_str(), std::ifstream::in);
    if (!manifest_file.is_open())
    {
        error("Unable to open manifest: %s", filename.c_str());
    }
    Samples sampleinfo;

    string line;
    vector<string> header;
    map<string, size_t> header_map;
    while (std::getline(manifest_file, line))
    {
        line = common::stringutil::replaceAll(line, "\n", "");
        line = common::stringutil::replaceAll(line, "#", "");
        if (line.empty())
        {
            continue;
        }

        if (header.empty())
        {
            common::stringutil::split(line, header, "\t,");
            static const set<string> legal_header_columns
                = { "id",    "path",        "index_path", "paragraph",      "idxdepth",
                    "depth", "read length", "sex",        "depth variance", "depth sd" };
            size_t j = 0;
            for (auto& h : header)
            {
                std::transform(h.begin(), h.end(), h.begin(), ::tolower);
                if (legal_header_columns.count(h) == 0)
                {
                    error("Unknown column %s in manifest", h.c_str());
                }
                header_map[h] = j++;
            }
            static const std::list<string> required_headers = { "id", "path" };
            for (const auto& rh : required_headers)
            {
                if (header_map.find(rh) == header_map.end())
                {
                    error("Required header %s not present in manifest", rh.c_str());
                }
            }
            if (!((header_map.count("idxdepth") != 0u)
                  || ((header_map.count("depth") != 0u) && (header_map.count("read length") != 0u))))
            {
                error("Manifest header must either specify index depth locations or depth and read length.");
            }
            continue;
        }

        vector<string> tokens;
        common::stringutil::split(line, tokens, "\t,");
        tokens.resize(header.size(), "");

        genotyping::SampleInfo sid;
        sid.set_sample_name(tokens[header_map["id"]]);

        auto find_file = [&sid, &filename](const std::string& bam_filename) -> std::string {
            if (bam_filename.compare(0, 5, "s3://") == 0 || bam_filename.compare(0, 7, "http://") == 0
                || bam_filename.compare(0, 8, "https://") == 0)
            {
                return bam_filename;
            }

            boost::filesystem::path p(bam_filename);
            if (!boost::filesystem::is_regular_file(p))
            {
                p = boost::filesystem::path(filename).parent_path() / p;
            }
            if (!boost::filesystem::is_regular_file(p))
            {
                error("Sample %s: File not found: %s", sid.sample_name().c_str(), p.c_str());
            }
            return p.string();
        };

        sid.set_filename(find_file(tokens[header_map["path"]]));
        if (header_map.count("index_path") != 0u)
        {
            sid.set_index_filename(find_file(tokens[header_map["index_path"]]));
        }

        double depth = -1;
        int read_length = -1;
        if ((header_map.count("depth") != 0u) && (header_map.count("read length") != 0u))
        {
            try
            {
                depth = stod(tokens[header_map["depth"]]);
                read_length = stoi(tokens[header_map["read length"]]);
            }
            catch (const std::exception&)
            {
            }
        }
        if ((depth < 0 || read_length < 0) && (header_map.count("idxdepth") != 0u))
        {
            std::string idxdepth_filename = tokens[header_map["idxdepth"]];
            boost::filesystem::path p(idxdepth_filename);
            if (!boost::filesystem::is_regular_file(p))
            {
                p = boost::filesystem::path(filename).parent_path() / p;
            }
            if (boost::filesystem::is_regular_file(p))
            {
                idxdepth_filename = p.string();
            }
            try
            {
                Json::Value idxdepth_json = common::getJSON(idxdepth_filename);
                if (read_length < 0 && idxdepth_json.isMember("read_length"))
                {
                    read_length = idxdepth_json["read_length"].asInt();
                }
                if (depth < 0 && idxdepth_json.isMember("autosome") && idxdepth_json["autosome"].isMember("depth"))
                {
                    depth = idxdepth_json["autosome"]["depth"].asDouble();
                }
            }
            catch (std::exception const& e)
            {
                if (!idxdepth_filename.empty())
                {
                    LOG()->warn(
                        "Cannot read idxdepth for sample %s: %s -- %s", sid.sample_name().c_str(),
                        idxdepth_filename.c_str(), e.what());
                }
            }
        }

        if (depth <= 0 || read_length <= 0)
        {
            error("No depth / read length estimate for sample %s", sid.sample_name().c_str());
        }

        sid.set_autosome_depth(depth);
        sid.set_read_length(static_cast<unsigned int>(read_length));

        if (header_map.count("depth sd") != 0u)
        {
            double depth_sd = 0;
            try
            {
                depth_sd = stod(tokens[header_map["depth sd"]]);
            }
            catch (const std::exception&)
            {
            }
            if (depth_sd <= 0)
            {
                error("Depth sd is not positive in sample %s", sid.sample_name().c_str());
            }
            sid.set_depth_sd(depth_sd);
        }
        else
        {
            if (header_map.count("depth variance") != 0u)
            {
                double depth_variance = 0;
                try
                {
                    depth_variance = stod(tokens[header_map["depth variance"]]);
                }
                catch (const std::exception&)
                {
                }
                if (depth_variance <= 0)
                {
                    error("Depth variance is not positive in sample %s", sid.sample_name().c_str());
                }
                double depth_sd = sqrt(depth_variance);
                sid.set_depth_sd(depth_sd);
            }
        }

        if (header_map.count("sex") != 0u)
        {
            string sex_string = tokens[header_map["sex"]];
            sid.set_sex(sex_string);
        }

        if (header_map.count("paragraph") != 0u)
        {
            std::string paragraph_filename = tokens[header_map["paragraph"]];
            boost::filesystem::path p(paragraph_filename);
            if (!boost::filesystem::is_regular_file(p))
            {
                p = boost::filesystem::path(filename).parent_path() / p;
            }
            if (boost::filesystem::is_regular_file(p))
            {
                paragraph_filename = p.string();
            }
            if (!paragraph_filename.empty())
            {
                try
                {
                    const Json::Value paragraph_json = common::getJSON(paragraph_filename);
                    sid.set_alignment_data(paragraph_json);
                }
                catch (std::exception const& e)
                {
                    if (!paragraph_filename.empty())
                    {
                        LOG()->warn(
                            "Cannot read paragraph JSON for sample %s: %s -- %s", sid.sample_name().c_str(),
                            paragraph_filename.c_str(), e.what());
                    }
                }
            }
        }
        sampleinfo.push_back(sid);
    }
    return sampleinfo;
}
}
