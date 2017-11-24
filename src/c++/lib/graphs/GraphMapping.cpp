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

#include "graphs/GraphMapping.hh"
#include <typeindex>

using std::map;
using std::string;

namespace graphs
{

Operation::Operation(std::string cigar, string query, string reference)
{
    string length_encoding = cigar;
    length_encoding.pop_back();
    char type_encoding = cigar.back();

    decodeOperation(type_encoding, std::stoi(length_encoding), std::move(query), std::move(reference));
    validate();
}

void Operation::decodeOperation(char type_encoding, int length, string query, string reference)
{
    query_ = std::move(query);
    reference_ = std::move(reference);

    switch (type_encoding)
    {
    case 'M':
        type_ = Type::kMatch;
        break;
    case 'N':
        type_ = Type::kMissingBases;
        break;
    case 'X':
        type_ = Type::kMismatch;
        break;
    case 'I':
        type_ = Type::kInsertionToRef;
        break;
    case 'D':
        type_ = Type::kDeletionFromRef;
        break;
    case 'S':
        type_ = Type::kSoftClipping;
        break;
    default:
        error("Error: %c is unknown CIGAR operation", type_encoding);
    }
    length_ = length;
}

void Operation::validate() const
{
    const bool full_length_query = query_.length() == (size_t)length_;
    const bool full_length_ref = reference_.length() == (size_t)length_;
    switch (type_)
    {
    case Type::kMatch:
        if (full_length_query && query_ == reference_)
            return;
        break;
    case Type::kMismatch:
        if (full_length_query && query_.length() == reference_.length())
        {
            bool found_matching_base = false;
            for (size_t index = 0; index != query_.length(); ++index)
            {
                if (query_[index] == reference_[index])
                    found_matching_base = true;
            }
            if (!found_matching_base)
                return;
        }
        break;
    case Type::kMissingBases:
        if (query_.length() == reference_.length() && full_length_query)
        {
            bool found_non_n_base_in_query_and_ref = false;
            for (size_t index = 0; index != query_.length(); ++index)
            {
                if (query_[index] != 'N' && reference_[index] != 'N')
                {
                    found_non_n_base_in_query_and_ref = true;
                }
            }
            if (!found_non_n_base_in_query_and_ref)
                return;
        }
        break;
    case Type::kDeletionFromRef:
        if (query_.empty() && !reference_.empty() && full_length_ref)
            return;
        break;
    case Type::kInsertionToRef:
        if (!query_.empty() && reference_.empty() && full_length_query)
            return;
        break;
    case Type::kSoftClipping:
        if (!query_.empty() && reference_.empty() && full_length_query)
            return;
        break;
    }

    error("Error: %s and %s are incompatible with operation %c", query_.c_str(), reference_.c_str(), asSymbol());
}

int32_t Operation::querySpan() const
{
    switch (type_)
    {
    case Type::kMatch:
    case Type::kMismatch:
    case Type::kMissingBases:
    case Type::kInsertionToRef:
    case Type::kSoftClipping:
        return length_;
    case Type::kDeletionFromRef:
        return 0;
    }
    return -1;
}

int32_t Operation::referenceSpan() const
{
    switch (type_)
    {
    case Type::kMatch:
    case Type::kMismatch:
    case Type::kMissingBases:
    case Type::kDeletionFromRef:
        return length_;
    case Type::kInsertionToRef:
    case Type::kSoftClipping:
        return 0;
    }
    return -1;
}

char Operation::asSymbol() const
{
    static const map<Operation::Type, char> kOperationToChar
        = { { Operation::Type::kMatch, 'M' },          { Operation::Type::kMismatch, 'X' },
            { Operation::Type::kInsertionToRef, 'I' }, { Operation::Type::kDeletionFromRef, 'D' },
            { Operation::Type::kSoftClipping, 'S' },   { Operation::Type::kMissingBases, 'N' } };
    return kOperationToChar.at(type_);
}

std::ostream& operator<<(std::ostream& os, const Operation& operation)
{

    os << operation.length() << operation.asSymbol() << "(" << operation.reference() << "->" << operation.query()
       << ")";
    return os;
}

Mapping::Mapping(
    int32_t reference_start, const std::string& cigar, const std::string& query, const std::string& reference)
{
    decodeOperations(reference_start, cigar, query, reference);
    updateCounts();
}

void Mapping::updateCounts()
{
    clipped_ = 0;
    matched_ = 0;
    mismatched_ = 0;
    missing_ = 0;
    inserted_ = 0;
    deleted_ = 0;
    for (size_t i = 0; i < num_operations(); ++i)
    {
        auto const& m = operations_[i];
        switch (m.type())
        {
        case Operation::Type::kSoftClipping:
            clipped_ += m.length();
            break;
        case Operation::Type::kMatch:
            matched_ += m.length();
            break;
        case Operation::Type::kMismatch:
            mismatched_ += m.length();
            break;
        case Operation::Type::kMissingBases:
            missing_ += m.length();
            break;
        case Operation::Type::kInsertionToRef:
            inserted_ += m.length();
            break;
        case Operation::Type::kDeletionFromRef:
            deleted_ += m.length();
            break;
        }
    }
}

void Mapping::decodeOperations(
    int32_t reference_start, const std::string& cigar, const std::string& query, const std::string& reference)
{
    auto logger = LOG();
    int32_t ref_pos = reference_start;
    int32_t query_pos = 0;
    string length_encoding;
    for (char c : cigar)
    {
        if (isalpha(c) != 0)
        {
            string query_piece;
            string reference_piece;
            int32_t operation_length = std::stoi(length_encoding);
            switch (c)
            {
            case 'M':
            case 'X':
            case 'N':
            {
                query_piece = query.substr((size_t)query_pos, (size_t)operation_length);
                reference_piece = reference.substr((size_t)ref_pos, (size_t)operation_length);
                query_pos += operation_length;
                ref_pos += operation_length;
#ifdef _DEBUG
                // Check for bugs in the aligner
                for (size_t pos_ = 0; pos_ < (size_t)operation_length; ++pos_)
                {
                    char this_op = 'M';
                    if (query_piece[pos_] == 'N' || reference_piece[pos_] == 'N')
                    {
                        this_op = 'N';
                    }
                    else if (query_piece[pos_] != reference_piece[pos_])
                    {
                        this_op = 'X';
                    }
                    if (c != this_op)
                    {
                        error(
                            "When decoding CIGAR string %s for %s vs. %s -- "
                            "at ref_pos %i and query_pos %i, the operation %c doesn't "
                            "match what we observe in the ref/query strings %s and %s!",
                            cigar.c_str(), reference.c_str(), query.c_str(), ref_pos, query_pos, c,
                            reference_piece.c_str(), query_piece.c_str());
                    }
                }
#endif
                break;
            }
            case 'I':
            case 'S':
                query_piece = query.substr((size_t)query_pos, (size_t)operation_length);
                query_pos += operation_length;
                break;
            case 'D':
                reference_piece = reference.substr((size_t)ref_pos, (size_t)operation_length);
                ref_pos += operation_length;
                break;
            default:
                error("Error: %c is unknown CIGAR operation", c);
            }
            operations_.emplace_back(c, operation_length, query_piece, reference_piece);
            length_encoding.clear();
        }
        else
        {
            if (isdigit(c) == 0)
            {
                error("Error: %s is malformed CIGAR string", cigar.c_str());
            }
            length_encoding += c;
        }
    }
}

string Mapping::query() const
{
    string query;
    for (const auto& operation : operations_)
    {
        if (operation.type() != Operation::Type::kSoftClipping)
        {
            query += operation.query();
        }
    }
    return query;
}

string Mapping::reference() const
{
    string reference;
    for (const auto& operation : operations_)
    {
        reference += operation.reference();
    }
    return reference;
}

int32_t Mapping::querySpan() const
{
    int32_t query_span = 0;
    for (const auto& operation : operations_)
    {
        query_span += operation.querySpan();
    }
    return query_span;
}

int32_t Mapping::referenceSpan() const
{
    int32_t reference_span = 0;
    for (const auto& operation : operations_)
    {
        reference_span += operation.referenceSpan();
    }
    return reference_span;
}

std::ostream& operator<<(std::ostream& os, const Mapping& mapping)
{
    for (size_t index = 0; index != mapping.num_operations(); ++index)
    {
        os << mapping[index];
    }

    return os;
}

NodeMapping::NodeMapping(int32_t reference_start, const string& node_cigar, const string& query, const Graph& graph)
{
    std::string nodeid_encoding;
    for (size_t index = 0; index != node_cigar.length(); ++index)
    {
        if (node_cigar[index] == '[')
        {
            node_id_ = std::stoull(nodeid_encoding);
            std::string cigar = node_cigar.substr(index + 1);
            cigar.pop_back();
            const string& reference = graph.nodes.at(node_id_)->sequence();
            decodeOperations(reference_start, cigar, query, reference);
            updateCounts();
            break;
        }
        if (isdigit(node_cigar[index]) == 0)
        {
            error("Error: %s is a malformed node CIGAR", node_cigar.c_str());
        }
        nodeid_encoding += node_cigar[index];
    }
}

std::ostream& operator<<(std::ostream& os, const NodeMapping& node_mapping)
{
    os << node_mapping.node_id() << '[' << static_cast<const Mapping&>(node_mapping) << ']';
    return os;
}

GraphMapping::GraphMapping(
    int32_t first_node_reference_start, const std::string& graph_cigar, const std::string& query, const Graph& graph)
    : first_node_reference_start_(first_node_reference_start)
{
    int32_t query_pos = 0;
    std::string node_cigar;
    for (size_t index = 0; index != graph_cigar.length(); ++index)
    {
        node_cigar += graph_cigar[index];
        if (node_cigar.back() == ']')
        {
            string query_piece = query.substr((size_t)query_pos);
            int32_t ref_pos = node_mappings_.empty() ? first_node_reference_start : 0;
            NodeMapping node_mapping(ref_pos, node_cigar, query_piece, graph);
            node_mappings_.push_back(node_mapping);
            query_pos += node_mapping.querySpan();
            node_cigar.clear();
        }
    }
}

string GraphMapping::query() const
{
    string query;
    for (const auto& node_mapping : node_mappings_)
    {
        query += node_mapping.query();
    }
    return query;
}

string GraphMapping::reference() const
{
    string reference;
    for (const auto& node_mapping : node_mappings_)
    {
        reference += node_mapping.reference();
    }
    return reference;
}

int32_t GraphMapping::querySpan() const
{
    int32_t query_span = 0;
    for (const auto& node_mapping : node_mappings_)
    {
        query_span += node_mapping.querySpan();
    }
    return query_span;
}

int32_t GraphMapping::referenceSpan() const
{
    int32_t reference_span = 0;
    for (const auto& node_mapping : node_mappings_)
    {
        reference_span += node_mapping.referenceSpan();
    }
    return reference_span;
}
}
