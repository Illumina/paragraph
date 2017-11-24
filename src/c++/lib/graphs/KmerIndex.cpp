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

#include "graphs/KmerIndex.hh"

#include <unordered_set>

#include "common/Error.hh"

namespace graphs
{

struct SimpleGraphMatchPos : GraphMatchPos
{
    SimpleGraphMatchPos() = default;

    SimpleGraphMatchPos(SimpleGraphMatchPos const& rhs)
        : node_(rhs.node_)
        , pos_(rhs.pos_)
        , length_(rhs.length_)
    {
    }

    SimpleGraphMatchPos(uint64_t node, int pos, size_t length)
        : node_(node)
        , pos_(pos)
        , length_(length)
    {
    }

    SimpleGraphMatchPos& operator=(SimpleGraphMatchPos const& rhs)
    {
        node_ = rhs.node_;
        pos_ = rhs.pos_;
        length_ = rhs.length_;
        return *this;
    }

    uint64_t node() const override { return node_; }
    int pos() const override { return pos_; }
    size_t length() const override { return length_; }

    uint64_t node_ = (uint64_t)-1;
    int pos_ = -1;
    size_t length_ = 0;
};

typedef std::list<SimpleGraphMatchPos> SimpleMatchList;
typedef std::unordered_map<std::string, SimpleMatchList> SimpleHashIndex;

static std::list<std::string> enumerateKmers(WalkableGraph const& g, int k, uint64_t node, int start_pos)
{
    std::list<std::string> results;
    if (k == 0)
    {
        results = { "" };
        return results;
    }

    auto this_node = g.node(node);
    if (start_pos + k <= (int)this_node->sequence().size())
    {
        results.emplace_back(
            this_node->sequence().substr(static_cast<unsigned long>(start_pos), static_cast<unsigned long>(k)));
    }
    else
    {
        for (const auto& s : g.succ(node))
        {
            auto suffixes = enumerateKmers(g, static_cast<int>(k - this_node->sequence().size() + start_pos), s, 0);
            for (const auto& suffix : suffixes)
            {
                results.emplace_back(this_node->sequence().substr(static_cast<unsigned long>(start_pos)) + suffix);
            }
        }
    }
    return results;
}

static void uniqueMatchPositions(
    WalkableGraph const& g, int k, uint64_t node, SimpleHashIndex& target,
    std::unordered_map<std::string, int>& all_kmer_counts)
{
    auto this_node = g.node(node);
    for (int start_pos = 0; start_pos < (int)this_node->sequence().size(); ++start_pos)
    {
        auto kmers = enumerateKmers(g, k, node, start_pos);
        std::unordered_map<std::string, int> kmer_counts;
        for (auto const& kmer : kmers)
        {
            auto k_it = kmer_counts.find(kmer);
            if (k_it == kmer_counts.end())
            {
                kmer_counts[kmer] = 1;
            }
            else
            {
                k_it->second++;
            }
        }
        for (auto const& kc : kmer_counts)
        {
            target[kc.first].emplace_back(node, start_pos, kc.first.size());
            // count kmers globally
            auto k_it = all_kmer_counts.find(kc.first);
            if (k_it == all_kmer_counts.end())
            {
                k_it = all_kmer_counts.emplace(kc.first, 1).first;
            }
            else
            {
                k_it->second = (k_it->second >= 0) ? k_it->second + kc.second : k_it->second - kc.second;
            }
            // if we have more than one match, make this a "bad" / multi-mapping kmer
            if (kc.second > 1)
            {
                k_it->second = -abs(k_it->second);
            }
        }
    }
}

static void buildKmerIndex(
    WalkableGraph const& g, int k, SimpleHashIndex& target, std::unordered_map<std::string, int>& all_kmer_counts)
{
    target.clear();

    std::unordered_set<uint64_t> visited;
    std::set<uint64_t> active = { g.source() };
    while (!active.empty())
    {
        auto current = *active.begin();
        uniqueMatchPositions(g, k, current, target, all_kmer_counts);
        visited.insert(current);
        for (auto const& next : g.succ(current))
        {
            if (visited.count(next) == 0 && active.count(next) == 0)
            {
                active.insert(next);
            }
        }
        active.erase(current);
    }

    for (auto const& a : g.allNodes())
    {
        assert(visited.count(a) != 0);
    }
}

struct KmerIndex::KmerIndexImpl
{
    KmerIndexImpl(WalkableGraph const& g, int k_)
        : graph(g)
        , k(k_)
    {
        buildKmerIndex(graph, k, index, kmer_counts);
    }
    WalkableGraph graph;
    int k;

    SimpleHashIndex index;
    std::unordered_map<std::string, int> kmer_counts;

    std::string searched;
    std::list<SimpleHashIndex::iterator> matches;
};

KmerIndex::KmerIndex(WalkableGraph const& g, int k)
    : _impl(new KmerIndexImpl(g, k))
{
}

KmerIndex::~KmerIndex() = default;

KmerIndex::KmerIndex(KmerIndex&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

KmerIndex& KmerIndex::operator=(KmerIndex&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

void KmerIndex::search(std::string const& str)
{
    // kmer index with hashing -> cannot search for elements shorter than k
    // we restrict search to length k but this can be extended to > k
    // by searching incrementally through the kmers of the string
    assert(str.size() == (unsigned)_impl->k);

    _impl->matches.clear();
    auto start_result = _impl->index.find(str.substr(0, (unsigned long)_impl->k));

    if (start_result == _impl->index.end())
    {
        return;
    }

    _impl->matches.push_back(start_result);
    _impl->searched = str;
}

uint64_t KmerIndex::count()
{
    assert(!_impl->searched.empty());
    size_t pos_count = 0;
    for (auto const& m : _impl->matches)
    {
        pos_count += m->second.size();
    }
    return pos_count;
}

std::list<std::unique_ptr<GraphMatchPos>> KmerIndex::matches()
{
    std::list<std::unique_ptr<GraphMatchPos>> matchlist;
    for (auto const& m : _impl->matches)
    {
        for (const auto& km : m->second)
        {
            matchlist.emplace_back(new SimpleGraphMatchPos(km));
        }
    }
    return matchlist;
}

int KmerIndex::kmerCount(std::string const& kmer) const
{
    auto it = _impl->kmer_counts.find(kmer);
    if (it != _impl->kmer_counts.end())
    {
        return abs(it->second);
    }
    return 0;
}

std::list<std::string> KmerIndex::kmers() const
{
    std::list<std::string> kmerlist;
    for (const auto& k : _impl->kmer_counts)
    {
        if (k.second >= 0)
        {
            kmerlist.push_back(k.first);
        }
    }
    return kmerlist;
}

std::list<std::string> KmerIndex::badkmers() const
{
    std::list<std::string> kmerlist;
    for (const auto& k : _impl->kmer_counts)
    {
        if (k.second < 0)
        {
            kmerlist.push_back(k.first);
        }
    }
    return kmerlist;
}
}
