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
#include "common/Genetics.hh"
#include "oligo/KmerGenerator.hh"
#include <iterator>
#include <numeric>

namespace grm
{

typedef std::size_t PositionType;
struct Candidate
{
    static const std::string DUMMY_PATH_ID_;
    Candidate()
        : pathId_(DUMMY_PATH_ID_)
        , position_(0)
        , reverse_(false)
        , mismatchCount_(0)
    {
    }
    Candidate(const std::string& pathId, PositionType position, bool reverse, unsigned mismatchCount)
        : pathId_(pathId)
        , position_(position)
        , reverse_(reverse)
        , mismatchCount_(mismatchCount)
    {
    }
    std::reference_wrapper<const std::string> pathId_;
    PositionType position_;
    bool reverse_;
    unsigned mismatchCount_;

    bool operator==(const Candidate& that)
    {
        return mismatchCount_ == that.mismatchCount_ && position_ == that.position_ && reverse_ == that.reverse_
            && pathId_.get() == that.pathId_.get();
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
                  << c.pathId_.get() << "pid)";
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

struct Path
{
    std::string path_id;
    std::string sequence_id;
    std::string path_sequence;
    typedef std::map<size_t, size_t> NodeStarts;
    NodeStarts starts;
    std::vector<std::string> node_names;
    std::vector<uint64_t> node_ids;

    KmerPositions kmerPositions_;

    friend std::ostream& operator<<(std::ostream& os, const Path& p)
    {
        return os << "Path(" << p.path_id << "," << p.sequence_id << "," << p.path_sequence << ","
                  << accumulate(p.path_sequence.begin(), p.path_sequence.end(), std::string("")) << ","
                  << ")";
    }
};

template <unsigned KMER_LENGTH> struct KmerAligner<KMER_LENGTH>::KmerAlignerImpl
{
    // position, mismatchCount

    typedef std::vector<Candidate> Candidates;

    std::list<Path> paths_;
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

    std::unordered_map<uint64_t, size_t> node_lengths_;
    void alignRead(common::Read& read) const;

    void align(
        const std::string& bases, const KmerPositions& fwKmerPositions, const std::string& rvBases,
        const KmerPositions& rvKmerPositions, const Path& path, Candidates& candidates) const;

    template <bool reverse>
    void align(
        const std::string& bases, const KmerPositions& seqeuenceKmerPositions, const Path& path,
        Candidates& candidates) const;

    void setGraph(graphs::Graph const& g, Json::Value const& paths);

private:
    template <bool overlap>
    static void
    makeKmers(std::string::const_iterator begin, std::string::const_iterator end, KmerPositions& kmerPositions);
    static Path makePath(
        const Json::Value& p, const graphs::Graph& g, const std::unordered_map<std::string, uint64_t>& node_id_map);

    static void dumpPaths(
        graphs::WalkableGraph const& wg, uint64_t currentId, uint64_t endId, Path currentPath, bool ref,
        std::list<Path>& path);
    bool updateAlignment(
        const Path& path, std::size_t pos, bool is_reverse_match, const std::string& bases,
        const std::string& rev_bases, common::Read& read) const;
    int buildCigar(
        int start, bool is_reverse_match, const std::string& rev_bases, const std::string& bases,
        typename Path::NodeStarts::const_iterator start_node, const Path& path, std::string& cigar,
        int& alignment_score) const;
    void pickBest(
        const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const;
};

const std::string Candidate::DUMMY_PATH_ID_;

template <unsigned KMER_LENGTH>
template <bool overlap>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::makeKmers(
    std::string::const_iterator begin, std::string::const_iterator end, KmerPositions& kmerPositions)
{
    kmerPositions.clear();
    oligo::KmerGenerator<KMER_LENGTH, unsigned, std::string::const_iterator> kmerGenerator(begin, end);
    unsigned kmer = 0;
    std::string::const_iterator position;
    while (kmerGenerator.next(kmer, position))
    {
        kmerPositions.push_back(KmerPosition(kmer, std::distance(begin, position)));
        if (!overlap)
        {
            kmerGenerator.skip(KMER_LENGTH - 1);
        }
    }

    std::sort(kmerPositions.begin(), kmerPositions.end());
}

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
            if (0 <= offset && path.path_sequence.size() >= offset + bases.size())
            {
                seedCandidates_.push_back(Candidate(path.path_id, offset, reverse, -1U));
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
            Candidate(ac.pathId_, ac.position_, ac.reverse_, countMismatches(bases, path.path_sequence, ac.position_)));
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
Path KmerAligner<KMER_LENGTH>::KmerAlignerImpl::makePath(
    const Json::Value& p, const graphs::Graph& g, const std::unordered_map<std::string, uint64_t>& node_id_map)
{
    assert(p.isMember("nodes"));
    assert(p.isMember("path_id"));
    assert(p.isMember("sequence"));
    Path path;
    path.path_id = p["path_id"].asString();
    path.sequence_id = p["sequence"].asString();
    for (const auto& n : p["nodes"])
    {
        const auto node_name = n.asString();
        const auto node_id = node_id_map.at(node_name);
        path.starts.emplace(path.path_sequence.size(), path.node_ids.size());
        auto node = g.nodes.find(node_id);
        assert(node != g.nodes.cend());
        path.path_sequence += node->second->sequence();
        path.node_ids.push_back(node_id);
        path.node_names.push_back(node_name);
    }
    makeKmers<true>(path.path_sequence.begin(), path.path_sequence.end(), path.kmerPositions_);
    return path;
}

template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::setGraph(graphs::Graph const& g, Json::Value const& paths)
{
    std::unordered_map<std::string, uint64_t> node_id_map;

    for (const auto& n : g.nodes)
    {
        assert(node_id_map.count(n.second->name()) == 0);
        node_id_map[n.second->name()] = n.first;
        node_lengths_[n.first] = n.second->sequence().size();
    }

    for (const auto& p : paths)
    {
        paths_.emplace_back(makePath(p, g, node_id_map));
    }

    // This has to be number of paths + 1 because we want to know if any of the paths
    // have more than 1 candidate for the best alignment
    // + 1 for heap push/pop
    candidates_.reserve(paths_.size() + 1 + 1);
}

char getCigarOp(const char s, const char r) { return s == r ? 'M' : s == 'N' ? 'N' : r == 'N' ? 'N' : 'X'; }

template <typename It> int calculateSoftClip(It& itRef, It& itSeq, std::size_t& length)
{
    int clip = 0;
    static const std::size_t MATCHES_IN_A_ROW_NEEDED = 5;
    int matchesNeeded = std::min(length, MATCHES_IN_A_ROW_NEEDED);
    int matchesInARow = 0;
    while (length && matchesInARow < matchesNeeded)
    {
        if (*itRef != *itSeq)
        {
            matchesInARow = 0;
        }
        else
        {
            ++matchesInARow;
        }
        ++itRef;
        ++itSeq;
        --length;
        ++clip;
    }
    if (matchesNeeded == matchesInARow)
    {
        clip -= matchesNeeded;
        itRef -= matchesNeeded;
        itSeq -= matchesNeeded;
        length += matchesNeeded;
    }
    return clip;
}

/**
 * \return CIGAR with mismatches clipped on left if first is set
 *         or on right if last is set. If the entire sequence is
 *         clipped, an empty CIGAR is returned
 */
std::string makeCigarBit(
    std::string::const_iterator itRef, std::string::const_iterator& itSeq, std::size_t length, const bool first,
    const bool last, int& pos)
{
    std::string ret;
    if (first)
    {
        int clip = calculateSoftClip(itRef, itSeq, length);
        if (clip)
        {
            ret += std::to_string(clip) + 'S';
            pos += clip;
        }
    }

    int endClip = 0;
    if (last)
    {
        std::reverse_iterator<std::string::const_iterator> ritRef(itRef + length);
        std::reverse_iterator<std::string::const_iterator> ritSeq(itSeq + length);
        endClip = calculateSoftClip(ritRef, ritSeq, length);
    }

    if (!length)
    {
        if (endClip)
        {
            ret += std::to_string(endClip) + 'S';
        }
        return ret;
    }

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
            }
            lastOp = op;
            lastOpLen = 0;
        }
    }
    if (lastOpLen)
    {
        ret += std::to_string(lastOpLen) + lastOp;
    }

    if (endClip)
    {
        ret += std::to_string(endClip) + 'S';
    }

    return ret;
}

template <unsigned KMER_LENGTH>
int KmerAligner<KMER_LENGTH>::KmerAlignerImpl::buildCigar(
    int start, bool is_reverse_match, const std::string& rev_bases, const std::string& bases,
    typename Path::NodeStarts::const_iterator start_node, const Path& path, std::string& cigar,
    int& alignment_score) const
{
    alignment_score = 0;
    auto length_left = bases.size();
    auto this_start = (size_t)((start));
    std::list<uint64_t> graph_nodes_traversed;
    std::string::const_iterator itSeq = (is_reverse_match ? rev_bases : bases).begin();
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
            const auto this_node_length = node_lengths_.find(path.node_ids[start_node->second])->second;
            assert(this_start + this_length <= this_node_length);
            graph_nodes_traversed.push_back(start_node->second);
            const bool first = cigar.empty();
            const std::string bit = makeCigarBit(
                path.path_sequence.begin() + this_start + start_node->first, itSeq, this_length, first,
                this_length == length_left, start);
            cigar += std::to_string(path.node_ids[start_node->second]) + "[" + bit + "]";
        }
        alignment_score += this_length;
        length_left -= this_length;
        start_node = next_step;
        this_start = 0;
    }
    return start;
}

template <unsigned KMER_LENGTH>
bool KmerAligner<KMER_LENGTH>::KmerAlignerImpl::updateAlignment(
    const Path& path, std::size_t pos, bool is_reverse_match, const std::string& bases, const std::string& rev_bases,
    common::Read& read) const
{
    typename Path::NodeStarts::const_iterator start_node = path.starts.lower_bound(pos);
    if (start_node == path.starts.end())
    {
        assert(path.starts.size() >= 1);
        start_node = std::prev(path.starts.end());
    }
    else if (start_node->first > pos)
    {
        start_node = std::prev(start_node);
    }

    assert(start_node != path.starts.end());
    std::string cigar;
    int alignment_score = 0;
    // buildCigar might introduce soft clips at the start
    int start = buildCigar(
        pos - start_node->first, is_reverse_match, rev_bases, bases, start_node, path, cigar, alignment_score);

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

    read.set_graph_cigar(cigar);
    read.set_graph_alignment_score(alignment_score);
    read.set_graph_mapq(60);
    read.set_is_graph_alignment_unique(true);
    read.set_graph_mapping_status(reads::MAPPED);
    return true;
}

template <unsigned KMER_LENGTH>
void KmerAligner<KMER_LENGTH>::KmerAlignerImpl::pickBest(
    const std::string& bases, const std::string& rvBases, const Candidates& candidates, common::Read& read) const
{
    const auto bestIt = std::min_element(candidates.begin(), candidates.end(), Candidate::lessMismatches);
    const Candidate& best = *bestIt;
    const std::string& alignPathId = best.pathId_;
    const Path& path = *std::find_if(
        paths_.begin(), paths_.end(), [&alignPathId](const Path& p) { return p.path_id == alignPathId; });
    updateAlignment(path, best.position_, best.reverse_, bases, rvBases, read);

    for (auto secondBestIt = std::min_element(bestIt + 1, candidates.end(), Candidate::lessMismatches);
         candidates.end() != secondBestIt;
         secondBestIt = std::min_element(secondBestIt + 1, candidates.end(), Candidate::lessMismatches))
    {
        const Candidate& secondBest = *secondBestIt;

        if (secondBest.mismatchCount_ == best.mismatchCount_)
        {
            common::Read secondBestRead = read;
            const std::string& secondBestPathId = secondBest.pathId_;
            const Path& secondBestPath
                = *std::find_if(paths_.begin(), paths_.end(), [&secondBestPathId](const Path& p) {
                      return p.path_id == secondBestPathId;
                  });
            updateAlignment(secondBestPath, secondBest.position_, secondBest.reverse_, bases, rvBases, secondBestRead);
            if (secondBestRead.graph_cigar() != read.graph_cigar() || secondBestRead.graph_pos() != read.graph_pos())
            {
                read.set_graph_mapq(0);
                read.set_is_graph_alignment_unique(false);
                read.set_graph_mapping_status(reads::BAD_ALIGN);
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
    read.set_graph_mapping_status(reads::UNMAPPED);
    candidates_.clear();
    const std::string bases = read.bases();
    fwKmerPositions_.clear();
    makeKmers<true>(bases.begin(), bases.end(), fwKmerPositions_);
    const auto rvBases = common::reverseComplement(bases);
    rvKmerPositions_.clear();
    makeKmers<true>(rvBases.begin(), rvBases.end(), rvKmerPositions_);
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
void KmerAligner<KMER_LENGTH>::setGraph(graphs::Graph const& g, Json::Value const& paths)
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
    mapped_ += reads::MAPPED == read.graph_mapping_status();
}

template class KmerAligner<16>;
// template class KmerAligner<32>;
}
