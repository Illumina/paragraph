// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Paragraph
// Copyright (c) 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//		http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
// See the License for the specific language governing permissions and limitations
//
//

/**
 *  \brief BAM alignment data structure
 *
 * \file BamAlignment.hh
 * \author Egor Dolzhenko
 * \email edolzhenko@illumina.com
 *
 */

#pragma once

#include <memory>
#include <string>

#include "json/json.h"

namespace common
{

// Holds a BAM alignment and decodes them from HTSLib structs.
class Read
{
public:
    enum MappingStatus
    {
        UNMAPPED = 0,
        MAPPED = 1,
        BAD_ALIGN = 2,
    };

    Read() = default;
    Read(Read const& rhs) = default;

    Read(std::string const& fragment_id, std::string const& bases, std::string const& quals)
    {
        set_fragment_id(fragment_id);
        set_bases(bases);
        set_quals(quals);
    }

    void setCoreInfo(std::string const& fragment_id, std::string const& bases, std::string const& quals)
    {
        set_fragment_id(fragment_id);
        set_bases(bases);
        set_quals(quals);
    }

    bool is_initialized() const { return !bases().empty(); }

    std::string const& fragment_id() const { return fragment_id_; };
    void set_fragment_id(std::string const& value) { fragment_id_ = value; };
    std::string const& bases() const { return bases_; };
    void set_bases(std::string const& value) { bases_ = value; };
    std::string const& quals() const { return quals_; };
    void set_quals(std::string const& value) { quals_ = value; };
    int32_t chrom_id() const { return chrom_id_; };
    void set_chrom_id(int32_t value) { chrom_id_ = value; };
    int32_t pos() const { return pos_; };
    void set_pos(int32_t value) { pos_ = value; };
    uint8_t mapq() const { return mapq_; };
    void set_mapq(uint8_t value) { mapq_ = value; };

    bool is_reverse_strand() const { return is_reverse_strand_; };
    void set_is_reverse_strand(bool value) { is_reverse_strand_ = value; };
    bool is_mate_reverse_strand() const { return is_mate_reverse_strand_; };
    void set_is_mate_reverse_strand(bool value) { is_mate_reverse_strand_ = value; };
    bool is_mapped() const { return is_mapped_; };
    void set_is_mapped(bool value) { is_mapped_ = value; };
    bool is_first_mate() const { return is_first_mate_; };
    void set_is_first_mate(bool value) { is_first_mate_ = value; };
    bool is_mate_mapped() const { return is_mate_mapped_; };
    void set_is_mate_mapped(bool value) { is_mate_mapped_ = value; };
    int32_t mate_chrom_id() const { return mate_chrom_id_; };
    void set_mate_chrom_id(int32_t value) { mate_chrom_id_ = value; };
    int32_t mate_pos() const { return mate_pos_; };
    void set_mate_pos(int32_t value) { mate_pos_ = value; };

    int32_t graph_pos() const { return graph_pos_; };
    void set_graph_pos(int32_t value) { graph_pos_ = value; };
    std::string const& graph_cigar() const { return graph_cigar_; };
    void set_graph_cigar(std::string value) { graph_cigar_ = value; };
    int32_t graph_mapq() const { return graph_mapq_; };
    void set_graph_mapq(int32_t value) { graph_mapq_ = value; };
    int32_t graph_alignment_score() const { return graph_alignment_score_; };
    void set_graph_alignment_score(int32_t value) { graph_alignment_score_ = value; };
    bool is_graph_alignment_unique() const { return is_graph_alignment_unique_; };
    void set_is_graph_alignment_unique(bool value) { is_graph_alignment_unique_ = value; };
    bool is_graph_reverse_strand() const { return is_graph_reverse_strand_; };
    void set_is_graph_reverse_strand(bool value) { is_graph_reverse_strand_ = value; };

    std::vector<std::string> const& graph_nodes_supported() const { return graph_nodes_supported_; };
    void add_graph_nodes_supported(std::string const& value) { graph_nodes_supported_.push_back(value); };
    std::string const& graph_nodes_supported(size_t pos) const { return graph_nodes_supported_[pos]; };
    void clear_graph_nodes_supported() { graph_nodes_supported_.clear(); };

    std::vector<std::string> const& graph_edges_supported() const { return graph_edges_supported_; };
    std::string const& graph_edges_supported(size_t pos) const { return graph_edges_supported_[pos]; };
    void add_graph_edges_supported(std::string const& value) { graph_edges_supported_.push_back(value); };
    void clear_graph_edges_supported() { graph_edges_supported_.clear(); };

    std::vector<std::string> const& graph_sequences_supported() const { return graph_sequences_supported_; };
    std::string const& graph_sequences_supported(size_t pos) const { return graph_sequences_supported_[pos]; };
    void add_graph_sequences_supported(std::string const& value) { graph_sequences_supported_.push_back(value); };
    void clear_graph_sequences_supported() { graph_sequences_supported_.clear(); };

    std::vector<std::string> const& graph_sequences_broken() const { return graph_sequences_broken_; };
    std::string const& graph_sequences_broken(size_t pos) const { return graph_sequences_broken_[pos]; };
    void add_graph_sequences_broken(std::string const& value) { graph_sequences_broken_.push_back(value); };
    void clear_graph_sequences_broken() { graph_sequences_broken_.clear(); };

    MappingStatus graph_mapping_status() const { return graph_mapping_status_; }
    void set_graph_mapping_status(MappingStatus status) { graph_mapping_status_ = status; }

    bool operator==(const Read& other) const
    {
        return fragment_id() == other.fragment_id() && bases() == other.bases() && quals() == other.quals()
            && chrom_id() == other.chrom_id() && pos() == other.pos() && mapq() == other.mapq()
            && is_mapped() == other.is_mapped() && is_first_mate() == other.is_first_mate()
            && is_mate_mapped() == other.is_mate_mapped() && mate_chrom_id() == other.mate_chrom_id()
            && mate_pos() == other.mate_pos();
    }

    Json::Value toJson() const
    {
        Json::Value val;

        if (!fragment_id_.empty())
            val["fragmentId"] = fragment_id_;
        if (!bases_.empty())
            val["bases"] = bases_;
        if (!quals_.empty())
            val["quals"] = quals_;
        if (chrom_id_)
            val["chromId"] = chrom_id_;
        if (pos_)
            val["pos"] = pos_;
        if (mapq_)
            val["mapq"] = mapq_;

        if (is_reverse_strand_)
            val["isReverseStrand"] = true;
        if (is_mate_reverse_strand_)
            val["isMateReverseStrand"] = true;
        if (is_mapped_)
            val["isMapped"] = true;
        if (is_first_mate_)
            val["isFirstMate"] = true;
        if (is_mate_mapped_)
            val["isMateMapped"] = true;
        if (mate_chrom_id_)
            val["mateChromId"] = mate_chrom_id_;
        if (mate_pos_)
            val["matePos"] = mate_pos_;

        if (graph_pos_)
            val["graphPos"] = graph_pos_;
        if (!graph_cigar_.empty())
            val["graphCigar"] = graph_cigar_;
        if (graph_mapq_)
            val["graphMapq"] = graph_mapq_;
        if (graph_alignment_score_)
            val["graphAlignmentScore"] = graph_alignment_score_;
        if (is_graph_alignment_unique_)
            val["isGraphAlignmentUnique"] = true;
        if (is_graph_reverse_strand_)
            val["isGraphReverseStrand"] = true;

        if (!graph_nodes_supported_.empty())
        {
            val["graphNodesSupported"] = Json::arrayValue;
            for (auto const& s : graph_nodes_supported_)
            {
                val["graphNodesSupported"].append(s);
            }
        }
        if (!graph_edges_supported_.empty())
        {
            val["graphEdgesSupported"] = Json::arrayValue;
            for (auto const& s : graph_edges_supported_)
            {
                val["graphEdgesSupported"].append(s);
            }
        }
        if (!graph_sequences_supported_.empty())
        {
            val["graphSequencesSupported"] = Json::arrayValue;
            for (auto const& s : graph_sequences_supported_)
            {
                val["graphSequencesSupported"].append(s);
            }
        }
        if (!graph_sequences_broken_.empty())
        {
            val["graphSequencesBroken"] = Json::arrayValue;
            for (auto const& s : graph_sequences_broken_)
            {
                val["graphSequencesBroken"].append(s);
            }
        }

        switch (graph_mapping_status_)
        {
        case BAD_ALIGN:
            val["graphMappingStatus"] = "BAD_ALIGN";
            break;
        case MAPPED:
            val["graphMappingStatus"] = "MAPPED";
            break;
        case UNMAPPED:
        default:
            break;
        }
        return val;
    }

private:
    std::string fragment_id_;
    std::string bases_;
    std::string quals_;
    int32_t chrom_id_ = -1;
    int32_t pos_ = -1;
    uint8_t mapq_ = 0;

    bool is_reverse_strand_ = false;
    bool is_mate_reverse_strand_ = false;
    bool is_mapped_ = false;
    bool is_first_mate_ = true;
    bool is_mate_mapped_ = false;
    int32_t mate_chrom_id_ = -1;
    int32_t mate_pos_ = -1;

    int32_t graph_pos_ = 0;
    std::string graph_cigar_;
    int32_t graph_mapq_ = 0;
    int32_t graph_alignment_score_ = 0;
    bool is_graph_alignment_unique_ = false;
    bool is_graph_reverse_strand_ = false;

    std::vector<std::string> graph_nodes_supported_;
    std::vector<std::string> graph_edges_supported_;
    std::vector<std::string> graph_sequences_supported_;
    std::vector<std::string> graph_sequences_broken_;

    MappingStatus graph_mapping_status_ = UNMAPPED;
};

typedef std::unique_ptr<Read> p_Read;
typedef std::vector<p_Read> ReadBuffer;

template <typename read_container> static inline ReadBuffer toReadBuffer(read_container const& reads)
{
    ReadBuffer rb;
    for (Read const& read : reads)
    {
        rb.emplace_back(new Read(read));
    }
    return rb;
}
}
