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

#include "grm/KmerAligner.hh"
#include "common/Error.hh"
#include "common/Klib.hh"
#include "oligo/KmerGenerator.hh"

#include "graphutils/SequenceOperations.hh"

#include <iterator>
#include <numeric>

namespace grm
{

namespace kmerAligner
{

    typedef std::size_t PositionType;
    struct Candidate
    {
        static const int DUMMY_PATH_ID_ = -1;
        Candidate()
            : pathId_(DUMMY_PATH_ID_)
            , position_(0)
            , reverse_(false)
            , mismatchCount_(0)
        {
        }
        Candidate(const int pathId, PositionType position, bool reverse, unsigned mismatchCount)
            : pathId_(pathId)
            , position_(position)
            , reverse_(reverse)
            , mismatchCount_(mismatchCount)
        {
        }
        int pathId_;
        PositionType position_;
        bool reverse_;
        unsigned mismatchCount_;

        bool operator==(const Candidate& that)
        {
            return mismatchCount_ == that.mismatchCount_ && position_ == that.position_ && reverse_ == that.reverse_
                && pathId_ == that.pathId_;
        }

        // push/pop_heap will result in the element with largest number of mismatches
        // evicted.
        static bool lessMismatches(const Candidate& left, const Candidate& right)
        {
            // Compare mismatchCount_ only as otherwise finding best/second best is
            // going to malfunction.
            return left.mismatchCount_ < right.mismatchCount_;
        }

        friend std::ostream& operator<<(std::ostream& os, const Candidate& c)
        {
            return os << "Candidate(" << c.position_ << "p," << c.reverse_ << "r," << c.mismatchCount_ << "mm,"
                      << c.pathId_ << "pid)";
        }
    };

    struct KmerPosition
    {
        typedef unsigned KmerType;

        KmerPosition()
            : kmer_(0)
            , position_(0)
        {
        }
        KmerPosition(KmerType kmer, PositionType position)
            : kmer_(kmer)
            , position_(position)
        {
        }
        KmerType kmer_;
        PositionType position_;

        bool operator<(const KmerPosition& that) const
        {
            return kmer_ < that.kmer_ || (kmer_ == that.kmer_ && position_ < that.position_);
        }

        friend std::ostream& operator<<(std::ostream& os, const KmerPosition& kp)
        {
            return os << "KmerPosition(" << oligo::bases(kp.kmer_) << "," << kp.position_ << "p)";
        }
    };

    typedef std::vector<KmerPosition> KmerPositions;

    template <unsigned KMER_LENGTH>
    void makeKmers(std::string::const_iterator begin, std::string::const_iterator end, KmerPositions& kmerPositions)
    {
        kmerPositions.clear();
        oligo::KmerGenerator<KMER_LENGTH, unsigned, std::string::const_iterator> kmerGenerator(begin, end);
        unsigned kmer = 0;
        std::string::const_iterator position;
        while (kmerGenerator.next(kmer, position))
        {
            kmerPositions.push_back(KmerPosition(kmer, static_cast<PositionType>(std::distance(begin, position))));
        }

        std::sort(kmerPositions.begin(), kmerPositions.end());
    }

    template <unsigned KMER_LENGTH> struct BasicPath
    {
        typedef std::map<size_t, graphtools::NodeId> NodeStarts;

        int pathId_;
        NodeStarts starts;

        KmerPositions kmerPositions_;

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
            // note seq() returns string by value. Make sure it does not get destroyed...
            const std::string& path_sequence = path_.seq();
            makeKmers<KMER_LENGTH>(path_sequence.begin(), path_sequence.end(), kmerPositions_);
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
    };

} // namespace kmerAligner

using namespace kmerAligner;

template <unsigned KMER_LENGTH> struct KmerAligner<KMER_LENGTH>::KmerAlignerImpl
{
    // position, mismatchCount

    typedef std::vector<Candidate> Candidates;

    typedef BasicPath<KMER_LENGTH> Path;
    std::vector<Path> paths_;
    // NOTE, these transient buffers make the class thread-unsafe but allow avoiding dynamic
    // memory allocations
    mutable KmerPositions fwKmerPositions_;
    mutable KmerPositions rvKmerPositions_;
    // best candidate alignments by count of mismatches.
    // capacity limits the number of equivalent candidates to keep.
    // This has to be number of paths + 1 because we want to know if any of the paths
    // have more than 1 candidate for the best alignment
    mutable Candidates candidates_;
    // this is just a transient buffer for single sequence/path seed alignments before
    // the duplicates are removed
    mutable std::vector<Candidate> seedCandidates_;

    graphtools::Graph const* graph_;
    void alignRead(common::Read& read) const;

    void align(
        const std::string& bases, const KmerPositions& fwKmerPositions, const std::string& rvBases,
        const KmerPositions& rvKmerPositions, const Path& path, Candidates& candidates) const;

    template <bool reverse>
    void align(
        const std::string& bases, const KmerPositions& sequenceKmerPositions, const Path& path,
        Candidates& candidates) const;

    void setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths);

private:
    bool updateAlignment(
        const Path& path, std::size_t pos, bool is_reverse_match, const std::string& bases,
        const std::string& rev_bases, common::Read& read) const;
    int buildCigar(
        int start, std::string::const_iterator itSeq, std::size_t lengthLeft,
        typename Path::NodeStarts::const_iterator startNode, const Path& path, std::size_t leftClip,
        std::size_t rightClip, std::string& cigar, int& alignment_score) const;
    void pickBest(
        const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const;
};

/**
 * \brief defines match for the purpose of the alignment.
 */
inline bool isMatch(const char readBase, const char referenceBase) { return readBase == referenceBase; }

inline bool isMismatch(const char readBase, const char referenceBase) { return !isMatch(readBase, referenceBase); }

unsigned countMismatches(const std::string& sequence, const std::string& reference, const int offset)
{
    assert(0 <= offset && reference.size() >= offset + sequence.size());

    const unsigned mismatches = std::inner_product(
        sequence.begin(), sequence.end(), reference.begin() + offset, 0U, std::plus<unsigned>(), &isMismatch);
    return mismatches;
}

template <unsigned KMER_LENGTH>
template <bool reverse>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::align(
    const std::string& bases, const KmerPositions& sequenceKmerPositions, const Path& path,
    Candidates& candidates) const
{
    seedCandidates_.clear();
    typename KmerPositions::const_iterator ppIt = path.kmerPositions_.begin();
    for (const auto& sp : sequenceKmerPositions)
    {
        while (path.kmerPositions_.end() != ppIt && ppIt->kmer_ < sp.kmer_)
        {
            ++ppIt;
        }

        while (path.kmerPositions_.end() != ppIt && ppIt->kmer_ == sp.kmer_)
        {
            const int offset = int(ppIt->position_) - int(sp.position_);
            // ignore candidates that overhang the path. If they are relevant
            // the path flanks should be made longer.
            if (0 <= offset && path.path_.seq().size() >= offset + bases.size())
            {
                seedCandidates_.push_back(Candidate(path.pathId_, offset, reverse, -1U));
            }
            ++ppIt;
        }
    }

    std::sort(seedCandidates_.begin(), seedCandidates_.end(), [](const Candidate& left, const Candidate& right) {
        return left.position_ < right.position_;
    });
    seedCandidates_.erase(
        std::unique(
            seedCandidates_.begin(), seedCandidates_.end(),
            [](const Candidate& left, const Candidate& right) { return left.position_ == right.position_; }),
        seedCandidates_.end());

    for (const auto& ac : seedCandidates_)
    {
        candidates.push_back(
            Candidate(ac.pathId_, ac.position_, ac.reverse_, countMismatches(bases, path.path_.seq(), ac.position_)));
        std::push_heap(candidates.begin(), candidates.end(), Candidate::lessMismatches);
        if (candidates.capacity() == candidates.size())
        {
            std::pop_heap(candidates.begin(), candidates.end(), Candidate::lessMismatches);
            candidates.pop_back();
        }
    }
}

template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::align(
    const std::string& bases, const KmerPositions& fwKmerPositions, const std::string& rvBases,
    const KmerPositions& rvKmerPositions, const Path& path, Candidates& candidates) const
{
    align<false>(bases, fwKmerPositions, path, candidates);
    align<true>(rvBases, rvKmerPositions, path, candidates);
}

template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::setGraph(
    graphtools::Graph const* g, std::list<graphtools::Path> const& paths)
{
    graph_ = g;
    for (const auto& p : paths)
    {
        paths_.emplace_back(Path(paths_.size(), p, g));
    }

    // This has to be number of paths + 1 because we want to know if any of the paths
    // have more than 1 candidate for the best alignment
    // + 1 for heap push/pop
    candidates_.reserve(paths_.size() + 1 + 1);
}

char getCigarOp(const char s, const char r) { return s == r ? 'M' : s == 'N' ? 'N' : r == 'N' ? 'N' : 'X'; }

/**
 * \brief Soft clip the ends that fall on N-containing source and sink
 */
template <typename It> int calculateSoftClip(It itRef, It itSeq, std::size_t length)
{
    return std::distance(itRef, std::find_if(itRef, itRef + length, [](const char c) { return 'N' != c; }));
}

/**
 * \param  matches will be incremented by the number of matches found
 * \return CIGAR with mismatches clipped on left if first is set
 *         or on right if last is set. If the entire sequence is
 *         clipped, an empty CIGAR is returned
 */
std::string
makeCigarBit(std::string::const_iterator itRef, std::string::const_iterator& itSeq, std::size_t length, int& matches)
{
    std::string ret;

    char lastOp = 0;
    std::size_t lastOpLen = 0;
    for (; length; ++itSeq, ++itRef, --length, ++lastOpLen)
    {
        const char op = getCigarOp(*itRef, *itSeq);
        if (op != lastOp)
        {
            if (lastOpLen)
            {
                ret += std::to_string(lastOpLen) + lastOp;
                if ('M' == lastOp)
                {
                    matches += lastOpLen;
                }
            }
            lastOp = op;
            lastOpLen = 0;
        }
    }
    if (lastOpLen)
    {
        ret += std::to_string(lastOpLen) + lastOp;
        if ('M' == lastOp)
        {
            matches += lastOpLen;
        }
    }

    return ret;
}

template <unsigned KMER_LENGTH>
int KmerAligner<KMER_LENGTH>::KmerAlignerImpl::buildCigar(
    int start, std::string::const_iterator itSeq, std::size_t length_left,
    typename Path::NodeStarts::const_iterator start_node, const Path& path, std::size_t leftClip,
    const std::size_t rightClip, std::string& cigar, int& alignment_score) const
{
    alignment_score = 0;
    auto this_start = (size_t)((start));
    while (start_node != path.starts.end() && length_left > 0)
    {
        auto this_length = length_left;
        auto next_step = std::next(start_node);
        if (next_step != path.starts.end())
        {
            this_length = std::min(length_left, next_step->first - start_node->first - this_start);
        }
        if (this_length > 0)
        {
            const auto this_node_length = graph_->nodeSeq(start_node->second).size();
            assert(this_start + this_length <= this_node_length);
            // note seq() returns string by value. Make sure it does not get destroyed...
            const std::string& path_sequence = path.path_.seq();
            std::string::const_iterator itRef = path_sequence.begin() + this_start + start_node->first;
            int matches = 0;
            const std::string bit = makeCigarBit(itRef, itSeq, this_length, matches);
            cigar += std::to_string(start_node->second) + "[";
            if (leftClip)
            {
                cigar += std::to_string(leftClip) + "S";
                leftClip = 0;
            }

            cigar += bit;

            if (rightClip && this_length == length_left)
            {
                cigar += std::to_string(rightClip) + "S";
            }

            cigar += "]";
            alignment_score += matches;
        }
        length_left -= this_length;
        start_node = next_step;
        this_start = 0;
    }
    LOG()->debug(cigar);
    return start;
}

template <unsigned KMER_LENGTH>
bool KmerAligner<KMER_LENGTH>::KmerAlignerImpl::updateAlignment(
    const Path& path, std::size_t pos, bool is_reverse_match, const std::string& bases, const std::string& rev_bases,
    common::Read& read) const
{
    std::string::const_iterator itSeq = (is_reverse_match ? rev_bases : bases).begin();
    // note seq() returns string by value. Make sure it does not get destroyed...
    const std::string& path_sequence = path.path_.seq();
    std::string::const_iterator itRef = path_sequence.begin() + pos;
    std::size_t length = bases.size();
    const int leftClip = calculateSoftClip(itRef, itSeq, length);
    pos += leftClip;

    std::reverse_iterator<std::string::const_iterator> ritRef(itRef + length);
    std::reverse_iterator<std::string::const_iterator> ritSeq(itSeq + length);
    const int rightClip = calculateSoftClip(ritRef, ritSeq, length - leftClip);

    typename Path::NodeStarts::const_iterator startNode = path.findStartNode(pos);
    assert(startNode != path.starts.end());

    std::string graphCigar;
    int alignment_score = 0;
    // buildCigar might introduce soft clips at the start
    int start = buildCigar(
        pos - startNode->first, itSeq + leftClip, length - leftClip - rightClip, startNode, path, leftClip, rightClip,
        graphCigar, alignment_score);

    read.set_graph_pos(start);

    if (is_reverse_match)
    {
        // we don't need to update reverse_strand_match_pos: the match
        // already is on the reverse strand because we didn't find one forward
        // and second matches to the reverse strand are checked below.
        read.set_bases(rev_bases);
        read.set_is_graph_reverse_strand(!read.is_reverse_strand());
    }
    else
    {
        // check if read also matches on the reverse strand
        read.set_is_graph_reverse_strand(read.is_reverse_strand());
    }

    read.set_graph_cigar(graphCigar);
    read.set_graph_alignment_score(alignment_score);
    read.set_graph_mapq(60);
    read.set_is_graph_alignment_unique(true);
    read.set_graph_mapping_status(common::Read::MAPPED);
    return true;
}

/**
 * \return best path id
 * \precondition non-empty candidates list
 */
template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::pickBest(
    const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const
{
    assert(!candidates.empty());
    const auto bestIt = std::min_element(candidates.begin(), candidates.end(), Candidate::lessMismatches);
    const Candidate& best = *bestIt;
    if (best.mismatchCount_ > 2)
    {
        return;
    }
    const Path& path = paths_[best.pathId_];
    updateAlignment(path, best.position_, best.reverse_, bases, rvBases, read);

    for (auto secondBestIt = std::min_element(bestIt + 1, candidates.end(), Candidate::lessMismatches);
         candidates.end() != secondBestIt;
         secondBestIt = std::min_element(secondBestIt + 1, candidates.end(), Candidate::lessMismatches))
    {
        const Candidate& secondBest = *secondBestIt;

        if (secondBest.mismatchCount_ == best.mismatchCount_)
        {
            common::Read secondBestRead = read;
            const Path& secondBestPath = paths_[secondBest.pathId_];
            updateAlignment(secondBestPath, secondBest.position_, secondBest.reverse_, bases, rvBases, secondBestRead);
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

template <unsigned KMER_LENGTH> void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::alignRead(common::Read& read) const
{
    read.set_graph_mapping_status(common::Read::UNMAPPED);
    candidates_.clear();
    const std::string bases = read.bases();
    fwKmerPositions_.clear();
    makeKmers<KMER_LENGTH>(bases.begin(), bases.end(), fwKmerPositions_);
    const auto rvBases = graphtools::reverseComplement(bases);
    rvKmerPositions_.clear();
    makeKmers<KMER_LENGTH>(rvBases.begin(), rvBases.end(), rvKmerPositions_);
    for (const auto& path : paths_)
    {
        align(bases, fwKmerPositions_, rvBases, rvKmerPositions_, path, candidates_);
    }

    if (!candidates_.empty())
    {
        pickBest(bases, rvBases, candidates_, read);
    }
}

template <unsigned KMER_LENGTH>
KmerAligner<KMER_LENGTH>::KmerAligner()
    : impl_(new KmerAlignerImpl())
{
}

template <unsigned KMER_LENGTH> KmerAligner<KMER_LENGTH>::~KmerAligner() = default;

template <unsigned KMER_LENGTH>
KmerAligner<KMER_LENGTH>::KmerAligner(KmerAligner&& rhs) noexcept
    : impl_(std::move(rhs.impl_))
{
}

template <unsigned KMER_LENGTH>
KmerAligner<KMER_LENGTH>& KmerAligner<KMER_LENGTH>::operator=(KmerAligner&& rhs) noexcept
{
    impl_ = std::move(rhs.impl_);
    return *this;
}

/**
 * Set the graph to align to
 * @param g a graph
 * @param paths list of paths
 */
template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::setGraph(graphtools::Graph const* g, std::list<graphtools::Path> const& paths)
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
template <unsigned KMER_LENGTH> void KmerAligner<KMER_LENGTH>::alignRead(common::Read& read)
{
    ++attempted_;
    impl_->alignRead(read);
    mapped_ += common::Read::MAPPED == read.graph_mapping_status();
}

template class KmerAligner<16>;
// template class KmerAligner<32>;
// for unit tests
template class KmerAligner<10>;
}
