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
 * \brief KmerAligner implementation
 *
 * \file KmerAligner.cpp
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#include "grm/KlibAligner.hh"
#include "common/Alignment.hh"
#include "common/Error.hh"
#include "common/Klib.hh"
#include "oligo/KmerGenerator.hh"

#include "graphutils/SequenceOperations.hh"

#include <iterator>
#include <numeric>

namespace grm
{

namespace klibAligner
{

    typedef std::size_t PositionType;
    struct Candidate
    {
        static const int DUMMY_PATH_ID_ = -1;
        Candidate()
            : pathId_(DUMMY_PATH_ID_)
            , position_(0)
            , reverse_(false)
            , score_(0)
        {
        }
        Candidate(const int pathId, PositionType position, bool reverse, int score, const std::vector<uint32_t>& cigar)
            : pathId_(pathId)
            , position_(position)
            , reverse_(reverse)
            , score_(score)
            , cigar_(cigar)
        {
        }
        int pathId_;
        PositionType position_;
        bool reverse_;
        int score_;
        std::vector<uint32_t> cigar_;

        bool operator==(const Candidate& that)
        {
            return score_ == that.score_ && position_ == that.position_ && reverse_ == that.reverse_
                && pathId_ == that.pathId_;
        }

        // push/pop_heap will result in the element with worst score evicted.
        static bool betterScore(const Candidate& left, const Candidate& right) { return left.score_ > right.score_; }

        friend std::ostream& operator<<(std::ostream& os, const Candidate& c)
        {
            return os << "Candidate(" << c.position_ << "p," << c.reverse_ << "r," << c.score_ << "s," << c.pathId_
                      << "pid)";
        }
    };

    struct BasicPath
    {
        typedef std::map<size_t, graphtools::NodeId> NodeStarts;

        int pathId_;
        NodeStarts starts;

        const graphtools::Path& path_;

        BasicPath(const int pathId, const graphtools::Path& p, const graphtools::Graph* g)
            : pathId_(pathId)
            , path_(p)
        {
            std::size_t nodeStart = 0;
            for (const auto& node_id : p.nodeIds())
            {
                starts.emplace(nodeStart, node_id);
                nodeStart += g->nodeSeq(node_id).length();
            }
        }

        typename NodeStarts::const_iterator findStartNode(std::size_t pos) const
        {
            // lower bound returns first node after start unless pos matches exactly, make sure path is not empty
            assert(starts.size() >= 1);
            typename NodeStarts::const_iterator ret = starts.lower_bound(pos);
            if (ret == starts.end())
            {
                ret = std::prev(starts.end());
            }
            else if (ret->first > pos)
            {
                ret = std::prev(ret);
            }

            return ret;
        }

        friend std::ostream& operator<<(std::ostream& os, const BasicPath& p)
        {
            return os << "Path(" << p.pathId_ << "id " << p.path_.seq() << ")";
        }
    };

} // namespace klibAligner

using namespace klibAligner;

struct KlibAlignerImpl
{
    // The defaults are here to match those in GraphAligner.cpp
    static const int match_ = 1;
    static const int mismatch_ = -4;
    // The ksw implementation counts the first base of indel as a combination of open and extension, so 5 instead of 6!
    static const int gapOpen_ = 5;
    static const int gapExtension_ = 1;

    // position, mismatchCount

    typedef std::vector<Candidate> Candidates;

    typedef BasicPath Path;
    std::vector<Path> paths_;

    // NOTE, these transient buffers make the class thread-unsafe but allow avoiding dynamic
    // memory allocations
    // best candidate alignments by count of mismatches.
    // capacity limits the number of equivalent candidates to keep.
    // This has to be number of paths + 1 because we want to know if any of the paths
    // have more than 1 candidate for the best alignment
    mutable Candidates candidates_;
    // this is just a transient buffer for single sequence/path seed alignments before
    // the duplicates are removed
    mutable std::vector<Candidate> seedCandidates_;
    mutable std::vector<common::KlibAlignment> klibAligners_;

    graphtools::Graph const* graph_ = 0;
    void alignRead(common::Read& read) const;
    template <bool reverse> void align(const std::string& sequence, const Path& path, Candidates& candidates) const;

    void setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths);

private:
    template <typename CigarIT>
    int buildGraphCigar(
        int start, std::string::const_iterator itSeq, typename Path::NodeStarts::const_iterator node, const Path& path,
        CigarIT pathCigarBegin, CigarIT pathCigarEnd, std::string& cigar, int& score) const;
    void pickBest(
        const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const;
    void updateAlignment(
        const Candidate& candidate, const std::string& rvBases, const std::string& bases, const Path& path,
        common::Read& read) const;
};

/**
 * \brief defines match for the purpose of the alignment.
 */
inline bool isMatch(const char readBase, const char referenceBase) { return readBase == referenceBase; }

void KlibAlignerImpl::setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths)
{
    paths_.clear();
    klibAligners_.clear();
    graph_ = g;
    const common::AlignmentParameters ap(match_, mismatch_, gapOpen_, gapExtension_);
    for (const graphtools::Path& p : paths)
    {
        paths_.emplace_back(Path(paths_.size(), p, g));

        klibAligners_.push_back(common::KlibAlignment());
        klibAligners_.back().setParameters(ap);
        klibAligners_.back().setRef(paths_.back().path_.seq().c_str());
    }

    // This has to be number of paths + 1 because we want to know if any of the paths
    // have more than 1 candidate for the best alignment
    // + 1 for heap push/pop
    candidates_.reserve(paths_.size() + 1 + 1);
}

template <typename CigarIT>
int KlibAlignerImpl::buildGraphCigar(
    int pathPos, std::string::const_iterator itSeq, typename Path::NodeStarts::const_iterator node, const Path& path,
    CigarIT pathCigarIt, const CigarIT pathCigarEnd, std::string& cigar, int& matches) const
{
    matches = 0;
    size_t nodePos = pathPos - node->first;
    const int ret = nodePos;
    cigar = std::to_string(node->second) + "[";
    while (pathCigarEnd != pathCigarIt)
    {
        const common::CigarCoder coder(*pathCigarIt);
        switch (coder.getCode())
        {
        case common::ALIGN:
        {
            std::size_t alignLength = coder.getLength();
            while (alignLength)
            {
                const std::size_t nodeAlignLength
                    = std::min<uint32_t>(alignLength, graph_->nodeSeq(node->second).size() - nodePos);
                assert(nodePos + nodeAlignLength <= graph_->nodeSeq(node->second).size());
                assert(nodePos + node->first <= path.path_.seq().size());

                const std::string& path_sequence = path.path_.seq();
                std::string::const_iterator itRef = path_sequence.begin() + nodePos + node->first;
                cigar += common::makeCigarBit(itRef, itSeq, nodeAlignLength, matches);

                alignLength -= nodeAlignLength;
                if (alignLength)
                {
                    node = std::next(node);
                    assert(path.starts.end() != node); // ran out of nodes with non I/S cigar components
                    nodePos = 0;
                    cigar += "]" + std::to_string(node->second) + "[";
                }
                else
                {
                    nodePos += nodeAlignLength;
                }
            }
            break;
        }

        case common::DELETE:
        {
            std::size_t delLength = coder.getLength();
            while (delLength)
            {
                const std::size_t nodeAlignLength
                    = std::min<std::size_t>(delLength, graph_->nodeSeq(node->second).size() - nodePos);
                assert(nodePos + nodeAlignLength <= graph_->nodeSeq(node->second).size());
                assert(nodePos + node->first <= path.path_.seq().size());

                if (nodeAlignLength)
                {
                    cigar += std::to_string(nodeAlignLength) + "D";
                }

                delLength -= nodeAlignLength;
                if (delLength)
                {
                    node = std::next(node);
                    assert(path.starts.end() != node); // ran out of nodes with non I/S cigar components
                    nodePos = 0;
                    cigar += "]" + std::to_string(node->second) + "[";
                }
                else
                {
                    nodePos += nodeAlignLength;
                }
            }
            break;
        }

        case common::INSERT:
        {
            cigar += std::to_string(coder.getLength()) + "I";
            itSeq += coder.getLength();
            break;
        }

        case common::SOFT_CLIP:
        {
            cigar += std::to_string(coder.getLength()) + "S";
            itSeq += coder.getLength();
            break;
        }

        default:
        {
            assert(false); // unexpected CIGAR component
        }
        }
        ++pathCigarIt;
    }

    cigar += "]";
    LOG()->trace("buildGraphCigar RET cigar={}", cigar);
    LOG()->trace("buildGraphCigar ret={}", ret);
    return ret;
}

void KlibAlignerImpl::updateAlignment(
    const Candidate& candidate, const std::string& rvBases, const std::string& bases, const Path& path,
    common::Read& read) const
{
    std::string pathCigar;
    int alignment_score = 0;
    // buildCigar might introduce soft clips at the start
    std::string::const_iterator itSeq = (candidate.reverse_ ? rvBases : bases).begin();
    typename Path::NodeStarts::const_iterator startNode = path.findStartNode(candidate.position_);
    const int start = buildGraphCigar(
        candidate.position_, itSeq, startNode, path, candidate.cigar_.begin(), candidate.cigar_.end(), pathCigar,
        alignment_score);

    if (candidate.reverse_)
    {
        // we don't need to update reverse_strand_match_pos: the match
        // already is on the reverse strand because we didn't find one forward
        // and second matches to the reverse strand are checked below.
        read.set_bases(rvBases);
        read.set_is_graph_reverse_strand(!read.is_reverse_strand());
    }
    else
    {
        // check if read also matches on the reverse strand
        read.set_is_graph_reverse_strand(read.is_reverse_strand());
    }

    read.set_graph_pos(start);
    read.set_graph_cigar(pathCigar);
    read.set_graph_alignment_score(alignment_score);
    read.set_graph_mapq(60);
    read.set_is_graph_alignment_unique(true);
    read.set_graph_mapping_status(common::Read::MAPPED);
}

/**
 * \return best path id
 * \precondition non-empty candidates list
 */
void KlibAlignerImpl::pickBest(
    const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const
{
    assert(!candidates.empty());
    const auto bestIt = std::min_element(candidates.begin(), candidates.end(), Candidate::betterScore);
    const Candidate& best = *bestIt;
    const Path& path = paths_[best.pathId_];

    updateAlignment(best, rvBases, bases, path, read);

    for (auto secondBestIt = std::min_element(bestIt + 1, candidates.end(), Candidate::betterScore);
         candidates.end() != secondBestIt;
         secondBestIt = std::min_element(secondBestIt + 1, candidates.end(), Candidate::betterScore))
    {
        const Candidate& secondBest = *secondBestIt;

        if (secondBest.score_ == best.score_)
        {
            common::Read secondBestRead = read;
            const Path& secondBestPath = paths_[secondBest.pathId_];

            updateAlignment(secondBest, rvBases, bases, secondBestPath, secondBestRead);

            if (secondBestRead.graph_cigar() != read.graph_cigar() || secondBestRead.graph_pos() != read.graph_pos())
            {
                read.set_graph_mapq(0);
                read.set_is_graph_alignment_unique(false);
                read.set_graph_mapping_status(common::Read::BAD_ALIGN);
                break;
            }
        }
        else
        {
            // no more as good candidates
            break;
        }
    }
}

template <bool reverse>
void KlibAlignerImpl::align(const std::string& sequence, const Path& path, Candidates& candidates) const
{
    common::KlibAlignment& klibAligner = klibAligners_.at(path.pathId_);
    klibAligner.setQuery(sequence.c_str());
    int r0 = -1, r1 = -1, a0 = -1, a1 = -1;
    int nCigar = -1;
    uint32_t* cigarIt = 0;
    klibAligner.getCigar(r0, r1, a0, a1, nCigar, cigarIt);
    LOG()->trace("KlibAlignerImpl::align r0={}, r1={}, a0={}, a1={}", r0, r1, a0, a1);

    if (r1 < r0)
    {
        // fully soft clipped alignment. Ignore
        return;
    }
    common::Cigar pathCigar;
    if (a0)
    {
        pathCigar.push_back(common::CigarCoder(a0, common::SOFT_CLIP).getValue());
    }
    pathCigar.insert(pathCigar.end(), cigarIt, cigarIt + nCigar);
    const uint32_t rightClip = sequence.length() - a1 - 1;
    if (rightClip)
    {
        const common::CigarCoder clip(rightClip, common::SOFT_CLIP);
        pathCigar.push_back(clip.getValue());
    }

    candidates.push_back(Candidate(path.pathId_, r0, reverse, klibAligner.getScore(), pathCigar));
    std::push_heap(candidates.begin(), candidates.end(), Candidate::betterScore);
    if (candidates.capacity() == candidates.size())
    {
        std::pop_heap(candidates.begin(), candidates.end(), Candidate::betterScore);
        candidates.pop_back();
    }
}

void KlibAlignerImpl::alignRead(common::Read& read) const
{
    read.set_graph_mapping_status(common::Read::UNMAPPED);
    candidates_.clear();
    const std::string bases = read.bases();
    const std::string rvBases = graphtools::reverseComplement(bases);
    for (const auto& path : paths_)
    {
        align<false>(bases, path, candidates_);
        align<true>(rvBases, path, candidates_);
    }

    if (!candidates_.empty())
    {
        pickBest(bases, rvBases, candidates_, read);
    }
}

KlibAligner::KlibAligner()
    : impl_(new KlibAlignerImpl())
{
}

KlibAligner::~KlibAligner() = default;

KlibAligner::KlibAligner(KlibAligner&& rhs) noexcept
    : impl_(std::move(rhs.impl_))
{
}

KlibAligner& KlibAligner::operator=(KlibAligner&& rhs) noexcept
{
    impl_ = std::move(rhs.impl_);
    return *this;
}

/**
 * Set the graph to align to
 * @param g a graph
 * @param paths list of paths
 */
void KlibAligner::setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths)
{
    impl_->setGraph(g, paths);
}

/**
 * Align a read to the graph and update the graph_* fields.
 *
 * We will align both the forward and the reverse strand and return
 * the alignment which is unique, or with the better score, or default
 * to the forward strand.
 *
 * @param read read structure
 */
void KlibAligner::alignRead(common::Read& read)
{
    ++attempted_;
    impl_->alignRead(read);
    mapped_ += common::Read::MAPPED == read.graph_mapping_status();
}
}
