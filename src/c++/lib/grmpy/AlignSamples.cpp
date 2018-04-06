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

#include "genotyping/SampleInfo.hh"
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
void alignSingleSample(
    const Parameters& parameters, const std::string& graphPath, const std::string& referencePath,
    common::BamReader& reader, genotyping::SampleInfo& sample)
{
    auto logger = LOG();
    // set up paragraph aligner
    const auto output_options = paragraph::Parameters::NODE_READ_COUNTS | paragraph::Parameters::EDGE_READ_COUNTS
        | paragraph::Parameters::PATH_READ_COUNTS | paragraph::Parameters::DETAILED_READ_COUNTS;
    paragraph::Parameters paragraph_parameters(
        parameters.max_reads(),
        parameters.max_reads() + 1, // disable variants: min reads for a variant > max reads we read
        0.01, parameters.bad_align_frac(), output_options, parameters.exact_sequence_matching(),
        parameters.graph_sequence_matching(), parameters.kmer_sequence_matching(), false);
    paragraph_parameters.set_threads(parameters.threads());
    paragraph_parameters.set_kmer_len(parameters.bad_align_uniq_kmer_len());

    logger->info("Loading parameters for sample {}", sample.sample_name());
    paragraph_parameters.load(graphPath, referencePath);
    logger->info("Done loading parameters");

    common::ReadBuffer all_reads;

    common::extractReads(reader, paragraph_parameters.target_regions(), parameters.max_reads(), all_reads);
    Json::Value output = paragraph::alignAndDisambiguate(paragraph_parameters, all_reads);
    output["bam"] = sample.filename();
    sample.set_alignment_data(output);
}
}
