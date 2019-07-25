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

#include <chrono>
#include <htslib/sam.h>
#include <regex>
#include <set>

#include "idxdepth/DepthEstimation.hh"

#include "common/BamReader.hh"
#include "common/Fasta.hh"
#include "common/StringUtil.hh"
#include "common/Timing.hh"

extern "C" {
#include <htslib/cram.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
};

#include "common/Error.hh"

namespace idxdepth
{

/**
 * Estimate index and filtered read depth from a BAM / CRAM file
 * @param parameters parameters object
 * @return depth information structure
 */
Json::Value estimateDepths(Parameters const& parameters)
{
    Json::Value output;
    common::FastaFile reference(parameters.reference_path());

    std::set<std::string> reference_chromosomes;
    for (const auto& chr : reference.getContigNames())
    {
        reference_chromosomes.insert(chr);
    }

    std::unique_ptr<htsFile, std::function<void(htsFile*)>> reads{ hts_open(parameters.bam_path().c_str(), "r"),
                                                                   hts_close };
    assert(reads);

    std::unique_ptr<bam_hdr_t, std::function<void(bam_hdr_t*)>> header{ sam_hdr_read(reads.get()), bam_hdr_destroy };
    assert(header);

    std::set<std::string> bam_chromosomes;
    std::set<std::string> both_chromosomes;

    std::set<std::string> autosome;
    std::set<std::string> sex_chromosomes;

    std::regex include_regex(parameters.include_regex());
    std::regex auto_regex(parameters.autosome_regex());
    std::regex sex_chromosome_regex(parameters.sex_chromosome_regex());

    for (int i = 0; i < header->n_targets; ++i)
    {
        const std::string contig = std::string(header->target_name[i]);

        if (!parameters.include_alt_contig() && contig.length() > 5)
        {
            continue;
        }

        if (parameters.include_regex().empty() || std::regex_match(contig, include_regex))
        {
            bam_chromosomes.insert(contig);
        }
        if (reference_chromosomes.count(contig) == 0)
        {
            error(
                "BAM does not match reference, chromosome %s is not present in %s.", contig.c_str(),
                parameters.reference_path().c_str());
        }
        else
        {
            if (header->target_len[i] != reference.contigSize(contig))
            {
                error(
                    "Contig lengths don't match for BAM (%i) and FASTA (%i)", header->target_len[i],
                    reference.contigSize(contig));
            }
            both_chromosomes.emplace(contig);
            if (std::regex_match(contig, auto_regex))
            {
                autosome.insert(contig);
            }
            if (std::regex_match(contig, sex_chromosome_regex))
            {
                sex_chromosomes.insert(contig);
            }
        }
    }
    if (both_chromosomes.size() != reference_chromosomes.size())
    {
        LOG()->warn("BAM header only has a subset of the reference chromosomes -- please make sure they match!");
    }

    std::map<std::string, size_t> reads_per_chromosome;
    if (common::stringutil::endsWith(parameters.bam_path(), ".bam"))
    {
        std::unique_ptr<hts_idx_t, std::function<void(hts_idx_t*)>> index{
            sam_index_load(reads.get(), (parameters.bam_path() + ".bai").c_str()), hts_idx_destroy
        };
        assert(index);

        reads_per_chromosome["*"] = 0;
        for (int i = 0; i < header->n_targets; ++i)
        {
            uint64_t u, v;
            hts_idx_get_stat(index.get(), i, &u, &v);
            reads_per_chromosome[header->target_name[i]] = u;
            reads_per_chromosome["*"] += v;
        }
        // Dump information about unmapped reads
        reads_per_chromosome["*"] += hts_idx_get_n_no_coor(index.get());
    }
    else
    {
        reads_per_chromosome["*"] = 0;
        for (int i = 0; i < header->n_targets; ++i)
        {
            reads_per_chromosome[header->target_name[i]] = 0;
        }
    }
    output["unaligned_reads"] = (Json::UInt64)reads_per_chromosome["*"];
    output["reference"] = parameters.reference_path();
    output["bam_path"] = parameters.bam_path();
    output["contigs"] = Json::arrayValue;

    std::unordered_map<std::string, std::unique_ptr<common::DepthInfo>> per_chromosome_depths;

    size_t read_length = 0;
    bool has_read_length = false;
    bool read_length_unique = true;
    std::mutex input_mutex;
    std::mutex output_mutex;
    int next_contig = 0;

    auto worker = [&]() {
        common::BamReader depthEstimator(
            parameters.bam_path(), parameters.bam_index_path(), parameters.reference_path());
        while (next_contig < header->n_targets)
        {
            int this_contig = -1;
            {
                std::lock_guard<std::mutex> input_lock(input_mutex);
                this_contig = next_contig;
                next_contig++;
            }
            if (this_contig >= 0 && this_contig < header->n_targets)
            {
                typedef std::chrono::duration<double, typename std::chrono::milliseconds::period> Milliseconds;
                using std::chrono::duration_cast;
                auto t0 = std::chrono::high_resolution_clock::now();
                const int i = this_contig;
                const std::string contig = header->target_name[i];
                std::ostringstream tid;
                tid << std::this_thread::get_id();

                if (bam_chromosomes.count(contig) == 0)
                {
                    LOG()->info("Thread {} skipping {}", tid.str(), contig);
                    continue;
                }

                LOG()->info("Thread {} estimating depth for {}", tid.str(), contig);

                Json::Value contig_info;
                contig_info["name"] = contig;
                contig_info["length"] = header->target_len[i];
                contig_info["non_n_length"] = (Json::UInt64)reference.contigNonNSize(contig);

                auto dp_info = depthEstimator.estimateDepth(header->target_name[i]);

                contig_info["depth"] = dp_info->depth_median;
                contig_info["depth_variance"] = dp_info->depth_variance;
                contig_info["reads_for_estimation"] = (Json::UInt64)dp_info->read_count;

                if (reads_per_chromosome[header->target_name[i]] > 0)
                {
                    contig_info["reads"] = (Json::UInt64)reads_per_chromosome[contig];
                    contig_info["index_depth"]
                        = (Json::UInt64)(dp_info->read_length * reads_per_chromosome[contig] / header->target_len[i]);
                }

                {
                    std::lock_guard<std::mutex> input_lock(output_mutex);
                    output["contigs"].append(contig_info);
                    if ((has_read_length && read_length != dp_info->read_length) || !dp_info->read_length_unique)
                    {
                        read_length_unique = false;
                    }
                    read_length = std::max(dp_info->read_length, read_length);
                    has_read_length = true;
                    per_chromosome_depths[contig] = std::move(dp_info);
                }

                auto t1 = std::chrono::high_resolution_clock::now();
                auto ticks = Milliseconds(t1 - t0);
                LOG()->info(
                    "Thread {} done estimating depth for {} ; DP = {} after {} us", tid.str(), contig,
                    contig_info["depth"].asDouble(), ticks.count());
            }
        }
    };

    std::vector<std::thread> workers;
    for (int i = 0; i < parameters.threads(); ++i)
    {
        workers.emplace_back(worker);
    }

    for (auto& w : workers)
    {
        w.join();
    }

    if (has_read_length)
    {
        output["read_length"] = (Json::UInt64)read_length;
    }

    if (!read_length_unique)
    {
        output["read_length_unique"] = false;
    }

    if (!autosome.empty())
    {
        Json::Value sc_info;
        double sc_depth = 0;
        size_t sc_length = 0;
        size_t sc_non_n_length = 0;
        sc_info["contigs"] = Json::arrayValue;
        for (const auto& contig : autosome)
        {
            if (bam_chromosomes.count(contig) == 0)
            {
                continue;
            }
            // TODO technically we should use mean but we don't compute it here.
            // TODO we can add further analysis for the case where median != mean also.
            sc_depth += reference.contigSize(contig) * per_chromosome_depths[contig]->depth_median;
            sc_length += reference.contigSize(contig);
            sc_non_n_length += reference.contigNonNSize(contig);
            sc_info["contigs"].append(contig);
        }
        sc_info["depth"] = sc_depth / sc_length;
        output["autosome"] = sc_info;
    }

    if (!sex_chromosomes.empty())
    {
        Json::Value sc_info;
        sc_info["contigs"] = Json::arrayValue;
        for (const auto& contig : sex_chromosomes)
        {
            if (bam_chromosomes.count(contig) == 0)
            {
                continue;
            }
            sc_info["contigs"].append(contig);
        }
        output["sex_chromosomes"] = sc_info;
    }

    return output;
}
}
