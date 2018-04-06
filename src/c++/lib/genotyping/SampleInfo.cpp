//
// Copyright (c) 2016 Illumina, Inc.
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

#include "genotyping/SampleInfo.hh"
#include "common/StringUtil.hh"

#include <fstream>
#include <list>
#include <map>
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
            static const set<string> legal_header_columns = {
                "id", "path", "index_path", "paragraph", "idxdepth", "depth", "read length",
            };
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
            const auto idxdepth_filename = tokens[header_map["idxdepth"]];
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

        if (header_map.count("paragraph") != 0u)
        {
            const auto paragraph_filename = tokens[header_map["paragraph"]];
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
