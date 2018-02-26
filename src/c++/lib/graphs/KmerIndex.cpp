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

#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>

#include "common/Error.hh"
#include "common/HashHelper.hh"

using std::list;
using std::string;
using std::unordered_set;
using std::vector;

namespace graphs
{

struct KmerIndex::KmerIndexImpl
{
    explicit KmerIndexImpl(StringToPathsMap kmer_to_paths_map);
    explicit KmerIndexImpl(const std::shared_ptr<WalkableGraph>& wgraph_ptr, int32_t kmer_len);
    void addKmerPathsStartingAtNode(const std::shared_ptr<WalkableGraph>& wgraph_ptr, int32_t node_id);
    void addKmerPaths(const std::list<GraphPath>& kmer_paths);
    void updateNodeKmerIndex();
    int32_t kmer_len;
    StringToPathsMap kmer_to_paths_map;
    std::unordered_map<uint32_t, size_t> node_kmer_counts;
    std::unordered_map<std::pair<uint32_t, uint32_t>, size_t> edge_kmer_counts;
};

KmerIndex::KmerIndexImpl::KmerIndexImpl(StringToPathsMap kmer_to_paths_map_)
    : kmer_to_paths_map(std::move(kmer_to_paths_map_))
{
    kmer_len = 0;
    for (const auto& kv : kmer_to_paths_map)
    {
        const string& kmer = kv.first;
        kmer_len = static_cast<int32_t>(kmer.length());
        break;
    }
    updateNodeKmerIndex();
}

KmerIndex::KmerIndexImpl::KmerIndexImpl(const std::shared_ptr<WalkableGraph>& wgraph_ptr, int32_t kmer_len_)
    : kmer_len(kmer_len_)
{
    const list<uint64_t> node_ids = wgraph_ptr->allNodes();
    for (uint64_t node_id : node_ids)
    {
        addKmerPathsStartingAtNode(wgraph_ptr, static_cast<int32_t>(node_id));
    }
    updateNodeKmerIndex();
}

void KmerIndex::KmerIndexImpl::addKmerPathsStartingAtNode(
    const std::shared_ptr<WalkableGraph>& wgraph_ptr, int32_t node_id)
{
    const string node_seq = wgraph_ptr->node(static_cast<uint64_t>(node_id))->sequence();
    vector<int32_t> node_list;
    node_list.push_back(node_id);
    for (size_t pos = 0; pos != node_seq.length(); ++pos)
    {
        GraphPath path(wgraph_ptr, pos, node_list, pos);
        addKmerPaths(path.extendBy(0, kmer_len - 1));
    }
}

void KmerIndex::KmerIndexImpl::addKmerPaths(const list<GraphPath>& kmer_paths)
{
    for (const GraphPath& kmer_path : kmer_paths)
    {
        kmer_to_paths_map[kmer_path.seq()].push_back(kmer_path);
    }
}

void KmerIndex::KmerIndexImpl::updateNodeKmerIndex()
{
    node_kmer_counts.clear();
    edge_kmer_counts.clear();
    for (const auto& kmer_and_paths : kmer_to_paths_map)
    {
        // kmer is unique
        if (kmer_and_paths.second.size() == 1)
        {
            bool has_previous = false;
            uint32_t previous_node = 0;
            for (auto const& path_node_id : kmer_and_paths.second.front().node_ids())
            {
                node_kmer_counts[path_node_id] += 1;
                if (has_previous)
                {
                    edge_kmer_counts[std::make_pair(previous_node, path_node_id)] += 1;
                }
                has_previous = true;
                previous_node = path_node_id;
            }
        }
    }
}

KmerIndex::KmerIndex(const StringToPathsMap& kmer_to_paths_map)
    : _impl(new KmerIndexImpl(kmer_to_paths_map))
{
}

KmerIndex::KmerIndex(std::shared_ptr<WalkableGraph> wgraph_ptr, int32_t kmer_len)
    : _impl(new KmerIndexImpl(wgraph_ptr, kmer_len))
{
}

KmerIndex::KmerIndex(const KmerIndex& other)
    : _impl(new KmerIndexImpl(*other._impl))
{
}

KmerIndex::KmerIndex(KmerIndex&& other) noexcept
    : _impl(std::move(other._impl))
{
}

KmerIndex& KmerIndex::operator=(const KmerIndex& other)
{
    if (this != &other)
    {
        _impl.reset(new KmerIndexImpl(*other._impl));
    }
    return *this;
}

KmerIndex& KmerIndex::operator=(KmerIndex&& other) noexcept
{
    _impl = std::move(other._impl);
    return *this;
}

KmerIndex::~KmerIndex() = default;

bool KmerIndex::operator==(const KmerIndex& other) const
{
    return (_impl->kmer_to_paths_map == other._impl->kmer_to_paths_map && _impl->kmer_len == other._impl->kmer_len);
}

static string encodePaths(const list<GraphPath>& paths)
{
    list<string> path_encodings;
    for (const auto& path : paths)
    {
        path_encodings.push_back(path.encode());
    }
    return boost::algorithm::join(path_encodings, ",");
}

string KmerIndex::encode() const
{
    list<string> kv_encodings;
    for (const auto& kv : _impl->kmer_to_paths_map)
    {
        const string encoding_of_paths = encodePaths(kv.second);
        const string kv_encoding = "{" + kv.first + "->" + encoding_of_paths + "}";
        kv_encodings.push_back(kv_encoding);
    }
    return boost::algorithm::join(kv_encodings, ",");
}

bool KmerIndex::contains(const std::string& kmer) const
{
    return _impl->kmer_to_paths_map.find(kmer) != _impl->kmer_to_paths_map.end();
}

size_t KmerIndex::numPaths(const std::string& kmer) const
{
    if (!contains(kmer))
    {
        return 0;
    }
    return _impl->kmer_to_paths_map.at(kmer).size();
}

const list<GraphPath>& KmerIndex::getPaths(const std::string& kmer) const { return _impl->kmer_to_paths_map.at(kmer); }

unordered_set<string> KmerIndex::getKmersWithNonzeroCount() const
{
    unordered_set<string> kmers_with_nonzero_count;
    for (const auto& kv : _impl->kmer_to_paths_map)
    {
        kmers_with_nonzero_count.insert(kv.first);
    }
    return kmers_with_nonzero_count;
}

size_t KmerIndex::numUniqueKmersOverlappingNode(uint32_t node_id) const
{
    auto node_it = _impl->node_kmer_counts.find(node_id);
    if (node_it != _impl->node_kmer_counts.end())
    {
        return node_it->second;
    }
    return 0;
}

size_t KmerIndex::numUniqueKmersOverlappingEdge(uint32_t from, uint32_t to) const
{
    auto edge_it = _impl->edge_kmer_counts.find(std::make_pair(from, to));
    if (edge_it != _impl->edge_kmer_counts.end())
    {
        return edge_it->second;
    }
    return 0;
}

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index)
{
    os << kmer_index.encode();
    return os;
}
}
