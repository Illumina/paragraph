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
 * Functions to run graph alignment in grmpy
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */
#include <iostream>

#include "common/ReadExtraction.hh"
#include "grmpy/AlignSamples.hh"

#include "paragraph/Disambiguation.hh"

#include <thread>

#include "common/Error.hh"

namespace grmpy
{

/**
 * Run single sample alignment
 * @param sample sample data structure
 */
void alignSingleSample(const Parameters& parameters, genotyping::SampleInfo& sample)
{
    auto logger = LOG();
    // set up paragraph aligner
    const auto output_options = paragraph::Parameters::NODE_READ_COUNTS | paragraph::Parameters::EDGE_READ_COUNTS
        | paragraph::Parameters::PATH_READ_COUNTS | paragraph::Parameters::DETAILED_READ_COUNTS;
    paragraph::Parameters paragraph_parameters(
        parameters.max_reads(),
        parameters.max_reads() + 1, // disable variants: min reads for a variant > max reads we read
        0.01, parameters.bad_align_frac(), output_options, true);
    paragraph_parameters.set_threads(parameters.alignment_threads());

    logger->info("Loading parameters for sample %s", sample.sample_name());
    paragraph_parameters.load(parameters.graph_path(), parameters.reference_path());
    logger->info("Done loading parameters");

    common::ReadBuffer all_reads;
    common::extractReads(
        sample.filename(), parameters.reference_path(), paragraph_parameters.target_regions(), parameters.max_reads(),
        all_reads);
    Json::Value output = paragraph::alignAndDisambiguate(paragraph_parameters, all_reads);
    output["bam"] = sample.filename();
    sample.set_alignment_data(output);
}

/**
 * main function for alignment
 *
 * @param parameters parameters for genotyping
 */
void alignSamples(Parameters& parameters)
{
    struct AlignTask
    {
        explicit AlignTask(genotyping::SampleInfo& si_)
            : si(si_)
        {
        }
        genotyping::SampleInfo& si;
    };

    std::list<AlignTask> tasks;
    for (auto& si : parameters.getSamples())
    {
        if (!si.get_alignment_data().isNull())
        {
            continue;
        }
        tasks.emplace_back(si);
    }

    auto next = tasks.begin();
    std::mutex tasks_mutex;

    auto align_worker = [&parameters, &tasks, &next, &tasks_mutex]() {
        bool work_left = true;
        while (work_left)
        {
            auto this_task = tasks.end();
            {
                std::lock_guard<std::mutex> tasks_lock(tasks_mutex);
                if (next == tasks.end())
                {
                    work_left = false;
                }
                else
                {
                    this_task = next++;
                }
            }
            if (this_task != tasks.end())
            {
                alignSingleSample(parameters, this_task->si);
            }
        }
    };

    std::list<std::thread> workers;
    for (int thread = 0; thread < std::min((int)tasks.size(), std::min(1, parameters.sample_threads())); ++thread)
    {
        workers.emplace_back(align_worker);
    }
    for (auto& worker : workers)
    {
        worker.join();
    }
}
}
