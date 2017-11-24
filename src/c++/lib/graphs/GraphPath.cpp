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

#include "graphs/GraphPath.hh"

#include <cassert>
#include <sstream>

#include "common/Error.hh"

using std::list;
using std::ostream;
using std::shared_ptr;
using std::string;
using std::to_string;

namespace graphs
{
struct GraphPath::Impl
{
    Impl(shared_ptr<WalkableGraph> wgraph_ptr, int32_t start_position, const list<int32_t>& nodes, int32_t end_position)
        : wgraph_ptr_(wgraph_ptr)
        , start_position_(start_position)
        , end_position_(end_position)
        , nodes_(nodes)
    {
    }
    bool isNodePositionValid(int32_t node_id, int32_t position) const;
    bool areNodesOrdered() const;
    bool arePositionsOrdered() const;
    bool isPathEmpty() const;
    bool isFirstNodePosValid() const;
    bool isLastNodePosValid() const;
    bool isPathConected() const;

    bool operator==(const Impl& other) const
    {
        return (wgraph_ptr_ == other.wgraph_ptr_) && (start_position_ == other.start_position_)
            && (end_position_ == other.end_position_) && (nodes_ == other.nodes_);
    }

    shared_ptr<WalkableGraph> wgraph_ptr_;
    int32_t start_position_;
    int32_t end_position_;
    list<int32_t> nodes_;
};

bool GraphPath::Impl::areNodesOrdered() const
{
    int32_t cur_node_id = nodes_.front();
    list<int32_t>::const_iterator node_id_iter = nodes_.begin();
    ++node_id_iter; // Assuming the path contains at least one node.
    while (node_id_iter != nodes_.end())
    {
        const int32_t next_node_id = *node_id_iter;
        if (cur_node_id >= next_node_id)
        {
            return false;
        }
        cur_node_id = next_node_id;
        ++node_id_iter;
    }
    return true;
}

bool GraphPath::Impl::arePositionsOrdered() const
{
    if (nodes_.size() == 1 && start_position_ > end_position_)
    {
        return false;
    }
    return true;
}

bool GraphPath::Impl::isNodePositionValid(int32_t node_id, int32_t position) const
{
    if (position < 0)
    {
        return false;
    }
    Node* node = wgraph_ptr_->node(node_id);
    const string& node_seq = node->sequence();
    if ((int32_t)node_seq.length() <= position)
    {
        return false;
    }
    return true;
}

bool GraphPath::Impl::isPathEmpty() const { return nodes_.empty(); }

bool GraphPath::Impl::isFirstNodePosValid() const
{
    const int32_t first_node_id = nodes_.front();
    return isNodePositionValid(first_node_id, start_position_);
}

bool GraphPath::Impl::isLastNodePosValid() const
{
    const int32_t last_node_id = nodes_.back();
    return isNodePositionValid(last_node_id, end_position_);
}

bool GraphPath::Impl::isPathConected() const
{
    list<int32_t>::const_iterator start_iter;
    list<int32_t>::const_iterator end_iter;
    for (start_iter = nodes_.begin(); start_iter != std::prev(nodes_.end()); ++start_iter)
    {
        end_iter = std::next(start_iter);
        if (!wgraph_ptr_->hasEdge(*start_iter, *end_iter))
        {
            return false;
        }
    }
    return true;
}

GraphPath::GraphPath(
    shared_ptr<WalkableGraph> wgraph_ptr, int32_t start_position, const list<int32_t>& nodes, int32_t end_position)
    : pimpl_(new Impl(wgraph_ptr, start_position, nodes, end_position))
{
}

GraphPath::~GraphPath() = default;

GraphPath::GraphPath(const GraphPath& other)
    : pimpl_(new Impl(*other.pimpl_))
{
}

GraphPath::GraphPath(GraphPath&& other) noexcept
    : pimpl_(std::move(other.pimpl_))
{
}

GraphPath& GraphPath::operator=(const GraphPath& other)
{
    if (this != &other)
    {
        pimpl_.reset(new Impl(*other.pimpl_));
    }
    return *this;
}

GraphPath& GraphPath::operator=(GraphPath&& other) noexcept
{
    pimpl_ = std::move(other.pimpl_);
    return *this;
}

int32_t GraphPath::start_position() const { return pimpl_->start_position_; }
int32_t GraphPath::end_position() const { return pimpl_->end_position_; }

string GraphPath::seq() const
{
    string path_seq;
    size_t node_index = 0;
    for (int32_t node_id : pimpl_->nodes_)
    {
        const Node* node = pimpl_->wgraph_ptr_->node(node_id);
        string node_seq = node->sequence();
        if (node_index == 0)
        {
            node_seq = node_seq.substr(pimpl_->start_position_);
        }

        if (node_index == pimpl_->nodes_.size() - 1)
        {
            const int32_t end_node_start = pimpl_->nodes_.size() == 1 ? pimpl_->start_position_ : 0;
            const int32_t segment_len = pimpl_->end_position_ - end_node_start + 1;
            node_seq = node_seq.substr(0, segment_len);
        }

        path_seq += node_seq;
        ++node_index;
    }
    return path_seq;
}

bool GraphPath::isValid() const
{
    return (
        !pimpl_->isPathEmpty() && pimpl_->isFirstNodePosValid() && pimpl_->isLastNodePosValid()
        && pimpl_->areNodesOrdered() && pimpl_->arePositionsOrdered() && pimpl_->isPathConected());
}

bool GraphPath::operator==(const GraphPath& other) const { return *pimpl_ == *other.pimpl_; }

string GraphPath::encode() const
{
    string path_encoding;

    size_t node_index = 0;
    const size_t last_index = pimpl_->nodes_.size() - 1;
    for (int32_t node_id : pimpl_->nodes_)
    {
        string node_encoding;
        if (node_index == 0) // Encoding first node.
        {
            node_encoding = "(" + to_string(node_id) + "@" + to_string(pimpl_->start_position_) + ")";
        }

        if (node_index == last_index) // Encoding last node.
        {
            node_encoding += "-(" + to_string(node_id) + "@" + to_string(pimpl_->end_position_) + ")";
        }

        if (node_index != 0 && node_index != last_index) // Encoding intermediate node.
        {
            node_encoding = "-(" + to_string(node_id) + ")";
        }
        path_encoding += node_encoding;
        ++node_index;
    }

    return path_encoding;
}

ostream& operator<<(ostream& os, const GraphPath& path) { return os << path.encode(); }

GraphPath GraphPath::extendStartPosition(int32_t extension_len) const
{
    const int32_t extended_first_node_pos = pimpl_->start_position_ - extension_len;
    GraphPath extended_path(pimpl_->wgraph_ptr_, extended_first_node_pos, pimpl_->nodes_, pimpl_->end_position_);
    if (!extended_path.isValid())
    {
        error("ERROR: Cannot extend %s left by %i", encode().c_str(), extension_len);
    }
    return extended_path;
}

GraphPath GraphPath::extendEndPosition(int32_t extension_len) const
{
    int32_t extended_last_node_pos = pimpl_->end_position_ + extension_len;
    GraphPath extended_path(pimpl_->wgraph_ptr_, pimpl_->start_position_, pimpl_->nodes_, extended_last_node_pos);
    if (!extended_path.isValid())
    {
        error("ERROR: Cannot extend %s right by %i", encode().c_str(), extension_len);
    }
    return extended_path;
}

GraphPath GraphPath::extendStartNodeTo(int32_t node_id) const
{
    list<int32_t> extended_nodes = pimpl_->nodes_;
    extended_nodes.insert(extended_nodes.begin(), node_id);
    Node* new_node = pimpl_->wgraph_ptr_->node(node_id);
    int32_t new_node_seq_len = new_node->sequence().length();
    GraphPath extended_path(pimpl_->wgraph_ptr_, new_node_seq_len - 1, extended_nodes, pimpl_->end_position_);
    if (!extended_path.isValid())
    {
        error("ERROR: Cannot extend %s left to node %i", encode().c_str(), node_id);
    }

    return extended_path;
}

GraphPath GraphPath::extendEndNodeTo(int32_t node_id) const
{
    list<int32_t> extended_nodes = pimpl_->nodes_;
    extended_nodes.push_back(node_id);
    GraphPath extended_path(pimpl_->wgraph_ptr_, pimpl_->start_position_, extended_nodes, 0);
    if (!extended_path.isValid())
    {
        error("ERROR: Cannot extend %s right to node %i", encode().c_str(), node_id);
    }

    return extended_path;
}

list<GraphPath> GraphPath::extendStartBy(int32_t extension_len) const
{
    list<GraphPath> extended_paths;

    const int32_t start_node_id = pimpl_->nodes_.front();
    // const int32_t start_node_length = wgraph_ptr_->node(start_node_id)->sequence().length();

    // Start position gives the maximum extension.
    if (extension_len <= pimpl_->start_position_)
    {
        extended_paths.push_back(extendStartPosition(extension_len));
    }
    else
    {
        const list<uint64_t> pred_node_ids = pimpl_->wgraph_ptr_->pred(start_node_id);
        const int32_t leftover_length = extension_len - pimpl_->start_position_ - 1;
        for (uint64_t pred_node_id : pred_node_ids)
        {
            const GraphPath path_with_this_node = extendStartNodeTo(pred_node_id);
            list<GraphPath> extensions_of_path_with_this_node = path_with_this_node.extendStartBy(leftover_length);
            extended_paths.splice(extended_paths.end(), extensions_of_path_with_this_node);
        }
    }

    return extended_paths;
}

list<GraphPath> GraphPath::extendEndBy(int32_t extension_len) const
{
    list<GraphPath> extended_paths;

    const int32_t end_node_id = pimpl_->nodes_.back();
    const int32_t end_node_length = pimpl_->wgraph_ptr_->node(end_node_id)->sequence().length();
    const int32_t max_extension_at_end_node = end_node_length - pimpl_->end_position_ - 1;

    if (extension_len <= max_extension_at_end_node)
    {
        extended_paths.push_back(extendEndPosition(extension_len));
    }
    else
    {
        const list<uint64_t> succ_node_ids = pimpl_->wgraph_ptr_->succ(end_node_id);
        const int32_t leftover_length = extension_len - max_extension_at_end_node - 1;
        for (uint64_t succ_node_id : succ_node_ids)
        {
            const GraphPath path_with_this_node = extendEndNodeTo(succ_node_id);
            list<GraphPath> extensions_of_path_with_this_node = path_with_this_node.extendEndBy(leftover_length);
            extended_paths.splice(extended_paths.end(), extensions_of_path_with_this_node);
        }
    }

    return extended_paths;
}

list<GraphPath> GraphPath::extendBy(int32_t start_extension_len, int32_t end_extension_len) const
{
    list<GraphPath> extended_paths;
    list<GraphPath> start_extended_paths = extendStartBy(start_extension_len);
    for (const GraphPath& path : start_extended_paths)
    {
        list<GraphPath> end_extended_paths = path.extendEndBy(end_extension_len);
        extended_paths.splice(extended_paths.end(), end_extended_paths);
    }
    return extended_paths;
}
}