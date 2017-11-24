//
// Copyright (c) 2016 Illumina, Inc.
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

#include "genotyping/GraphGenotyper.hh"
#include "common/Error.hh"
#include "genotyping/PopStats.hh"
#include "graphs/Graph.hh"
#include <fstream>
#include <iostream>
#include <utility>

using std::string;
using std::vector;

namespace genotyping
{
struct GraphGenotyper::GraphGenotyperImpl
{
    GraphGenotyperImpl(
        const double genotype_error_rate, const int min_overlap_bases, const int max_read_times,
        const std::string reference_allele_name)
        : genotype_error_rate(genotype_error_rate)
        , min_overlap_bases(static_cast<int32_t>(min_overlap_bases))
        , max_read_times(static_cast<int32_t>(max_read_times))
        , reference_allele_name(reference_allele_name){};

    // name of sequences (alleles)
    std::vector<std::string> sequence_names;

    // basic information for final output
    Json::Value basic_json_information;

    Idxdepth sample_idxdepth;

    /**
     *  genotyping result storage
     *  order is the same as sample_names and idx_stats_
     */
    std::map<std::string, BreakpointStat> breakpoints_map; // BP key -> info & gt & stats

    /**
     *  whole-variant genotypes by BP vote
     */
    std::vector<VariantGenotype> graph_genotypes;

    /**
     *  population-scale stats for whole-variant genotypes across all samples
     */
    double call_rate;
    double zero_count_rate;
    double pass_rate;
    int num_valid_samples;
    std::vector<double> allele_frequencies; // order the same oas sequence_names_
    double hwe_pval;

    // parameters
    const double genotype_error_rate;
    const int32_t min_overlap_bases;
    const int32_t max_read_times;
    const std::string reference_allele_name;
};

GraphGenotyper::GraphGenotyper(
    const double genotype_error_rate, const int min_overlap_bases, const int max_read_times,
    const std::string reference_allele_name)
    : _impl(new GraphGenotyperImpl(genotype_error_rate, min_overlap_bases, max_read_times, reference_allele_name))
{
}

GraphGenotyper::~GraphGenotyper() = default;

void GraphGenotyper::genotypeGraph(
    const string& paragraph_input_path, const string& reference_path, const string& manifest_path, bool use_em)
{
    auto logger = LOG();

    // read sample info
    _impl->sample_idxdepth.load(manifest_path);
    logger->info("Done loading bam stats from manifest.");

    // read graph and read count info
    loadGraphAndCounts(paragraph_input_path, reference_path);
    logger->info("Loaded paragraph input JSON.");

    // breakpoint genotyping
    logger->info("Running breakpoint genotyping...");
    try
    {
        computeBreakpointGenotypes(use_em);
    }
    catch (const std::exception& e)
    {
        throw e.what();
    }
    logger->info("Breakpoint genotyping completed. Calculating whole-variant variant genotypes...");

    // whole variant genotyping
    try
    {
        computeVariantGenotypes();
    }
    catch (const std::exception& e)
    {
        throw e.what();
    }
    logger->info("Whole-variant genotyping completed.");
}

void GraphGenotyper::loadGraphAndCounts(const string& paragraph_input_path, const string& reference_path)
{
    // open & read json
    Json::Reader reader;
    Json::Value paragraph_json;
    std::ifstream graph_desc(paragraph_input_path);
    reader.parse(graph_desc, paragraph_json);
    if (!paragraph_json.isMember("samples"))
    {
        throw std::logic_error(
            "Missing sample information in paragraph json. Raw paragraph output needs to be merged through "
            "merge_paragraph_json.py before going into the genotyper");
    }
    graph_desc.close();

    // build graph
    graphs::Graph graph;
    graphs::fromJson(paragraph_json, reference_path, graph);
    graphs::WalkableGraph wgraph(graph);

    // load basic json info for output
    _impl->basic_json_information["eventinfo"] = paragraph_json["eventinfo"];
    _impl->basic_json_information["graphinfo"] = Json::objectValue;
    vector<string> copied_keys = { "ID", "target_regions", "sequencenames"};
    for (auto& key : copied_keys)
    {
        _impl->basic_json_information["graphinfo"][key] = paragraph_json[key];
    }
    _impl->basic_json_information["graphinfo"]["nodes"] = Json::arrayValue;
    for (auto const& n : paragraph_json["nodes"])
    {
        Json::Value node = Json::objectValue;
        node["name"] = n["name"];
        if (n.isMember("sequences"))
        {
            node["sequences"] = n["sequences"];
        }
        _impl->basic_json_information["graphinfo"]["nodes"].append(node);
    }
    _impl->basic_json_information["graphinfo"]["edges"] = Json::arrayValue;
    for (auto const& e : paragraph_json["edges"])
    {
        Json::Value edge = Json::objectValue;
        edge["name"] = e["from"].asString() + "_" + e["to"].asString();
        if (e.isMember("sequences"))
        {
            edge["sequences"] = e["sequences"];
        }
        _impl->basic_json_information["graphinfo"]["edges"].append(edge);
    }

    // load allele info from JSON
    loadSeqNames(paragraph_json);
    if (_impl->sequence_names.empty())
    {
        throw std::logic_error("Unable to load allele (sequencenames) from paragraph json. Please check your input!");
    }

    // load breakpoints & counts
    loadBreakpointInfo(wgraph, paragraph_json);
    loadEdgeCounts(paragraph_json);
}

void GraphGenotyper::loadSeqNames(Json::Value& paragraph_json)
{
    vector<string> raw_seq_names;
    for (auto& seq_name : paragraph_json["sequencenames"])
    {
        raw_seq_names.push_back(seq_name.asString());
    }
    // now try to put reference in position 0
    int ref_index = -1;
    for (int i = 0; i < (int)raw_seq_names.size(); i++)
    {
        if (raw_seq_names[i] == _impl->reference_allele_name)
        {
            ref_index = i;
            break;
        }
    }
    if (ref_index == -1)
    {
        std::cerr << "Warning: no reference allele specified. Use input allele order. This might be problematic "
                     "for brekapoint genotyping"
                  << std::endl;
        _impl->sequence_names = raw_seq_names;
    }
    else
    {
        _impl->sequence_names.push_back(_impl->reference_allele_name);
        for (int i = 0; i < (int)raw_seq_names.size(); i++)
        {
            if (i == ref_index)
            {
                continue;
            }
            _impl->sequence_names.push_back(raw_seq_names[i]);
        }
    }
}

void GraphGenotyper::loadBreakpointInfo(graphs::WalkableGraph& wgraph, Json::Value& paragraph_json)
{
    bool source_exist = wgraph.nodeName(wgraph.source()) == "source";
    bool sink_exist = wgraph.nodeName(wgraph.sink()) == "sink";
    if ((source_exist && !sink_exist) || (!source_exist && sink_exist))
    {
        throw std::logic_error("Bad graph: only have source or sink node.");
    }

    std::map<string, vector<uint64_t>> edge_to_seq_indexes = generateEgeNameToSeqIndexMap(paragraph_json, source_exist);

    for (auto node : wgraph.allNodes())
    {
        if (source_exist)
        {
            if (node == wgraph.source() || node == wgraph.sink())
            {
                continue;
            }
        }

        vector<const vector<uint64_t>*> raw_neighbors;
        raw_neighbors.resize(2);
        const auto pred = wgraph.pred(node);
        const vector<uint64_t> pred_vec = { pred.begin(), pred.end() };
        const auto succ = wgraph.succ(node);
        const vector<uint64_t> succ_vec = { succ.begin(), succ.end() };
        raw_neighbors[0] = &pred_vec;
        raw_neighbors[1] = &succ_vec;

        string node_name = wgraph.nodeName(node);
        for (size_t index = 0; index < 2; index++)
        {
            const vector<uint64_t> neighbors = source_exist
                ? removeSourceSink(*raw_neighbors[index], wgraph.source(), wgraph.sink())
                : *raw_neighbors[index];
            if (neighbors.size() > 1)
            {
                bool from_node_fixed = (index != 0);
                string node_key = from_node_fixed ? (node_name + "_") : ("_" + node_name);
                vector<string> neighbor_names;
                for (auto& neighbor : neighbors)
                {
                    neighbor_names.push_back(wgraph.nodeName(neighbor));
                }
                BreakpointStat bp_element(node_name, from_node_fixed, neighbor_names, edge_to_seq_indexes);
                _impl->breakpoints_map.insert(std::make_pair(node_key, bp_element));
            }
        }
    }
}

std::map<string, vector<uint64_t>>
GraphGenotyper::generateEgeNameToSeqIndexMap(Json::Value& paragraph_json, bool source_exist)
{
    std::map<string, vector<uint64_t>> edge_name_to_seq_indexes; // edge name --> sequence index

    std::map<string, uint64_t> seq_name_to_index; // sequence name --> sequence index
    for (size_t i = 0; i < _impl->sequence_names.size(); i++)
    {
        seq_name_to_index[_impl->sequence_names[i]] = (uint64_t)i;
    }

    for (auto edge : paragraph_json["edges"])
    {
        if (!edge.isMember("sequences"))
        {
            if (!source_exist)
            {
                throw std::logic_error("Missing sequence label at non-sink or non-source node!");
            }
        }
        string edge_name = edge["from"].asString() + "_" + edge["to"].asString();
        for (auto& seq_name : edge["sequences"])
        {
            edge_name_to_seq_indexes[edge_name].push_back(seq_name_to_index[seq_name.asString()]);
        }
    }
    return edge_name_to_seq_indexes;
}

vector<uint64_t> GraphGenotyper::removeSourceSink(const vector<uint64_t>& node_vec, uint64_t source, uint64_t sink)
{
    std::vector<uint64_t> new_vec;
    for (auto v : node_vec)
    {
        if (v != source && v != sink)
        {
            new_vec.push_back(v);
        }
    }
    return new_vec;
}

void GraphGenotyper::loadEdgeCounts(Json::Value& paragraph_json)
{
    for (auto& bp_stat : _impl->breakpoints_map)
    {
        for (int i = 0; i < (int)_impl->sample_idxdepth.sampleSize(); i++)
        {
            bp_stat.second.addBreakpointEdgeCountForOneSample(paragraph_json, _impl->sample_idxdepth.getSampleName(i));
        }
    }
}

void GraphGenotyper::computeBreakpointGenotypes(bool use_em)
{
    for (auto& bp_stat : _impl->breakpoints_map)
    {
        bp_stat.second.genotype(
            _impl->genotype_error_rate, (int)_impl->sample_idxdepth.sampleSize(), _impl->max_read_times,
            _impl->min_overlap_bases, _impl->sample_idxdepth, use_em);
    }
}

void GraphGenotyper::computeVariantGenotypes()
{
    // whole variant genotyping
    for (size_t sample = 0; sample < _impl->sample_idxdepth.sampleSize(); sample++)
    {
        VariantGenotype sample_genotype = VariantGenotype();
        for (auto& breakpoint : _impl->breakpoints_map)
        {
            sample_genotype.addInfoFromSingleBreakpoint(breakpoint.second.getGenotype(sample));
        }
        sample_genotype.genotype();
        _impl->graph_genotypes.push_back(sample_genotype);
    }
    // update stats
    if (_impl->graph_genotypes.size() > 1)
    {
        updateVariantStats();
    }
}

void GraphGenotyper::updateVariantStats()
{
    int num_missing = 0;
    int num_zero_count = 0;
    int num_pass = 0;
    for (auto& sample_genotype : _impl->graph_genotypes)
    {
        if (sample_genotype.empty())
        {
            num_missing++;
            if (sample_genotype.zeroCount())
            {
                num_zero_count++;
            }
            continue;
        }
        if (sample_genotype.pass())
        {
            num_pass++;
        }
    }
    _impl->num_valid_samples = (int)_impl->graph_genotypes.size() - num_missing;
    _impl->call_rate = (double)_impl->num_valid_samples / _impl->graph_genotypes.size();
    _impl->zero_count_rate = (double)num_zero_count / _impl->graph_genotypes.size();
    _impl->pass_rate = (double)num_pass / _impl->graph_genotypes.size();

    PopStats pop_stats(_impl->graph_genotypes);
    std::map<uint64_t, int> raw_allele_counts = pop_stats.alleleCounts();
    _impl->allele_frequencies.resize(_impl->sequence_names.size(), 0);
    for (auto& allele : raw_allele_counts)
    {
        _impl->allele_frequencies[allele.first] = (double)allele.second / _impl->num_valid_samples / 2;
    }
    _impl->hwe_pval = pop_stats.getHWE();
}

void GraphGenotyper::toCsv(std::ostream* out)
{
    *out << "#Sample\tGenotype" << std::endl;
    for (size_t sample = 0; sample < _impl->sample_idxdepth.sampleSize(); sample++)
    {
        string sample_name = _impl->sample_idxdepth.getSampleName((int)sample);
        *out << sample_name << "," << string(_impl->graph_genotypes[sample]) << std::endl;
    }
}

void GraphGenotyper::toJson(std::ostream* out)
{
    Json::Value output_json = _impl->basic_json_information;
    output_json["sequencenames"] = Json::arrayValue;
    for (auto& seq_name : _impl->sequence_names)
    {
        output_json["sequencenames"].append(seq_name);
    }

    if (_impl->graph_genotypes.size() > 1)
    {
        output_json["call_rate"] = _impl->call_rate;
        output_json["zero_count_rate"] = _impl->zero_count_rate;
        output_json["pass_rate"] = _impl->pass_rate;
        output_json["AF"] = Json::Value();
        for (size_t i = 0; i < _impl->allele_frequencies.size(); i++)
        {
            output_json["AF"][_impl->sequence_names[i]] = _impl->allele_frequencies[i];
        }
        output_json["hwe_p"] = _impl->hwe_pval;
        output_json["breakpoints"] = Json::Value();
        for (auto& bp_stat : _impl->breakpoints_map)
        {
            output_json["breakpoints"][bp_stat.first] = bp_stat.second.descriptionsToJson();
        }
    }

    // sample-specific information
    output_json["samples"] = Json::Value();
    vector<string> edge_names_to_print = getUniqueEdgeNames();

    for (int sample = 0; sample < (int)_impl->sample_idxdepth.sampleSize(); sample++)
    {
        string sample_name = _impl->sample_idxdepth.getSampleName((int)sample);
        output_json["samples"][sample_name] = Json::Value();
        output_json["samples"][sample_name] = _impl->graph_genotypes[sample].toJson(_impl->sequence_names);
        for (auto& key_name : { "read_counts_by_edge", "breakpoints" })
        {
            output_json["samples"][sample_name][key_name] = Json::Value();
        }

        for (auto& bp_stat : _impl->breakpoints_map)
        {
            for (auto& edge_name : edge_names_to_print)
            {
                auto current_edge_count = bp_stat.second.getEdgeCount(sample, edge_name);
                if (current_edge_count > 0)
                {
                    output_json["samples"][sample_name]["read_counts_by_edge"][edge_name]
                        = static_cast<int>(current_edge_count);
                }
            }
        }

        for (auto& bp_stat : _impl->breakpoints_map)
        {
            output_json["samples"][sample_name]["breakpoints"][bp_stat.first]
                = bp_stat.second.sampleGenotypeToJson(sample, _impl->sequence_names);
        }
    }
    Json::StyledStreamWriter w;
    w.write(*out, output_json);
}

vector<string> GraphGenotyper::getUniqueEdgeNames()
{
    std::map<string, bool> unique_edge_map;
    for (auto& bp_stat : _impl->breakpoints_map)
    {
        vector<string> edge_names = bp_stat.second.exportEdgeNames();
        for (auto& e_name : edge_names)
        {
            if (unique_edge_map.find(e_name) == unique_edge_map.end())
            {
                unique_edge_map[e_name] = true;
            }
        }
    }
    vector<string> unique_edge_names;
    for (auto& element : unique_edge_map)
    {
        unique_edge_names.push_back(element.first);
    }
    return unique_edge_names;
}
}
