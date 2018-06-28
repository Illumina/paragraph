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
 * \brief Align / seed reads via graph kmer index
 *
 * \file PathAligner.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "grm/PathAligner.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "common/Alignment.hh"

#include <algorithm>
#include <list>

#include <boost/range/adaptor/filtered.hpp>

#include "common/Error.hh"

namespace grm
{

struct ExactMatch
{
    size_t qpos;
    graphtools::Path path;
    bool isReverse;
};

struct PathAligner::Impl
{
    int32_t kmerSize = 32;
    std::unique_ptr<graphtools::KmerIndex> pKmerIndex;
};

PathAligner::PathAligner(int32_t kmer_size)
    : impl_(new Impl())
{
    impl_->kmerSize = kmer_size;
}
PathAligner::~PathAligner() = default;
PathAligner::PathAligner(PathAligner&& rhs) noexcept = default;
PathAligner& PathAligner::operator=(PathAligner&& rhs) noexcept = default;

void PathAligner::setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const&)
{
    impl_->pKmerIndex.reset(new graphtools::KmerIndex(*g, impl_->kmerSize));
}

void PathAligner::alignRead(common::Read& read)
{
    ++attempted_;

    const auto kmer_length = impl_->pKmerIndex->kmerLength();

    const size_t read_length = read.bases().size();
    if (read_length < kmer_length)
    {
        return;
    }

    std::list<ExactMatch> matches;
    for (int strand = 0; strand < 2; ++strand)
    {
        const bool is_reverse_strand = strand != 0;
        const std::string& read_bases = is_reverse_strand ? graphtools::reverseComplement(read.bases()) : read.bases();

        for (size_t pos = 0; pos + kmer_length <= read_bases.size(); ++pos)
        {
            const std::string kmer = read_bases.substr(pos, kmer_length);

            if (impl_->pKmerIndex->numPaths(kmer) == 1)
            {
                size_t qpos = pos;
                const auto extended
                    = graphtools::extendPathMatching(impl_->pKmerIndex->getPaths(kmer).front(), read_bases, qpos);
                matches.push_back(ExactMatch{ qpos, extended, is_reverse_strand });
                pos = matches.back().qpos + matches.back().path.length();
            }
        }
    }

    if (!matches.empty())
    {
        ++anchored_;
    }
    using boost::adaptors::filtered;
    auto filtered_mems
        = matches | filtered([read_length](ExactMatch const& mem) -> bool { return mem.path.length() == read_length; });

    if (filtered_mems.empty())
    {
        return;
    }

    auto const& mem_to_translate = filtered_mems.front();

    if (mem_to_translate.isReverse)
    {
        read.set_bases(graphtools::reverseComplement(read.bases()));
        read.set_is_graph_reverse_strand(true);
    }
    else
    {
        read.set_is_graph_reverse_strand(false);
    }

    std::string cigar;
    if (mem_to_translate.qpos > 0)
    {
        cigar += std::to_string(mem_to_translate.qpos) + "S";
    }
    cigar += std::to_string(mem_to_translate.path.length()) + "M";
    if (mem_to_translate.qpos + mem_to_translate.path.length() < read_length)
    {
        cigar += std::to_string(read_length - mem_to_translate.qpos - mem_to_translate.path.length()) + "S";
    }
    graphtools::Alignment linearAlignment(0, cigar);
    graphtools::GraphAlignment graphAlignment
        = graphtools::projectAlignmentOntoGraph(linearAlignment, mem_to_translate.path);

    read.set_graph_alignment_score(mem_to_translate.path.length());
    read.set_graph_cigar(graphAlignment.generateCigar());
    read.set_graph_pos(mem_to_translate.path.startPosition());

    read.set_graph_mapping_status(common::Read::MAPPED);
    if (std::next(filtered_mems.begin()) != filtered_mems.end())
    {
        read.set_is_graph_alignment_unique(false);
        read.set_graph_mapq(0);
    }
    else
    {
        read.set_is_graph_alignment_unique(true);
        read.set_graph_mapq(60);
    }

    ++mapped_;
}
}