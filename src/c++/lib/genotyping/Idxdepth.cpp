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

#include "genotyping/Idxdepth.hh"
#include <fstream>
#include <iterator>
#include <sstream>

using std::string;
using std::vector;

namespace genotyping
{
Idxdepth::Idxdepth(int read_length, int autosome_depth, int sample_size)
{
    if (sample_size > 0)
    {
        for (int i = 0; i < sample_size; i++)
        {
            string sample_name = "Dummy" + std::to_string(i + 1);
            sample_names_.push_back(sample_name);
        }
        SingleIdxdepth sid;
        sid.read_length = read_length;
        sid.autosome_depth = autosome_depth;
        idx_stats_.resize(sample_size, sid);
    }
}

void Idxdepth::load(const string& manifest_path)
{
    std::ifstream manifest_file;
    manifest_file.open(manifest_path, std::ifstream::in);
    if (!manifest_file.is_open())
    {
        string message = "Unable to open manifest:\n" + manifest_path;
        throw std::ios_base::failure(message);
    }
    string line;
    while (std::getline(manifest_file, line))
    {
        if (line[0] == '#')
        {
            try
            {
                is_good_manifest_header(line);
            }
            catch (const std::exception& e)
            {
                throw e.what();
            }
            continue;
        }
        std::istringstream iss(line);
        vector<string> tokens;
        copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(), back_inserter(tokens));
        if (tokens.size() != num_tokens_)
        {
            throw std::logic_error("Abnormal field count in manifest.");
        }
        struct SingleIdxdepth sid;
        try
        {
            sid.autosome_depth = stod(tokens[depth_token_index_]);
            sid.read_length = static_cast<int32_t>(stoi(tokens[read_len_token_index_]));
        }
        catch (const std::exception& e)
        {
            throw e.what();
        }
        sample_names_.push_back(tokens[0]);
        idx_stats_.push_back(sid);
    }
    manifest_file.close();
    if (idx_stats_.empty())
    {
        throw std::logic_error("Didn't load any sample info from manifest file.");
    }
}

void Idxdepth::is_good_manifest_header(string& line)
{
    std::stringstream iss(line);
    vector<string> tokens;
    copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(), back_inserter(tokens));
    if (tokens.size() != num_tokens_)
    {
        throw std::logic_error("Valid header should have 4 columns with: ID, Path, Depth, Read length.");
    }
}
}
