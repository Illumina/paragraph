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
 * Genotyper for the graph-represented whole variants across many samples
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/GraphGenotyper.hh"
#include "GraphGenotyperImpl.hh"
#include "genotyping/BreakpointFinder.hh"
#include "genotyping/BreakpointStatistics.hh"
#include "genotyping/GenotypeSet.hh"
#include "genotyping/PopulationStatistics.hh"

#include "common/Error.hh"

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using std::list;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::unordered_map;
using std::vector;

namespace genotyping
{

using graphtools::Graph;

GraphGenotyper::GraphGenotyper()
    : _impl(new GraphGenotyperImpl())
{
}

GraphGenotyper::~GraphGenotyper() = default;

/**
 * Set the graph we genotype on
 * @param graph our graph to genotype
 */
void GraphGenotyper::reset(Graph const* graph)
{
    // reset all counts
    _impl.reset(new GraphGenotyperImpl());
    _impl->graph = graph;

    // work out allele and edge names
    const auto bp_map = createBreakpointMap(*_impl->graph);
    set<string> allele_names;
    for (const auto& bp : bp_map)
    {
        _impl->breakpointnames.push_back(bp.first);
        auto const& bp_info = bp.second;
        for (auto const& an : bp_info.canonicalAlleleNames())
        {
            allele_names.insert(an);
        }
    }
    _impl->allelenames.resize(allele_names.size());
    std::copy(allele_names.begin(), allele_names.end(), _impl->allelenames.begin());
}

/**
 * @return the graph (asserts if no graph is set)
 */
Graph const& GraphGenotyper::getGraph() const
{
    assert(_impl->graph);
    return *_impl->graph;
}

/**
 * Add a alignment, depth and summary statistics for a sample
 */
void GraphGenotyper::addAlignment(SampleInfo const& sampleinfo)
{
    const std::string& samplename = sampleinfo.sample_name();
    const Json::Value& alignment = sampleinfo.get_alignment_data();
    const double depth = sampleinfo.autosome_depth();
    const int read_length = sampleinfo.read_length();

    _impl->samplenames.push_back(samplename);
    _impl->samplenameindex[samplename] = _impl->samplenames.size() - 1;

    // add breakpoint map and counts
    _impl->breakpoint_maps.push_back(createBreakpointMap(*_impl->graph));
    for (auto& breakpoint : _impl->breakpoint_maps.back())
    {
        breakpoint.second.addCounts(alignment);
    }
    _impl->depths.emplace_back(depth, read_length);
    _impl->sexes.emplace_back(sampleinfo.sex());

    // extract extra information and check we have the same event
    if (alignment.isMember("eventinfo"))
    {
        if (_impl->basic_info.isMember("eventinfo"))
        {
            assert(alignment["eventinfo"] == _impl->basic_info["eventinfo"]);
        }
        else
        {
            _impl->basic_info["eventinfo"] = alignment["eventinfo"];
        }
    }

    if (!_impl->basic_info.isMember("graphinfo"))
    {
        _impl->basic_info["graphinfo"] = Json::objectValue;

        // load event ID
        if (alignment.isMember("ID"))
        {
            _impl->basic_info["graphinfo"]["ID"] = alignment["ID"];
        }
        else if (alignment.isMember("vcf_records"))
        {
            // comma separate multiple vcf record ID together
            std::string event_id = "";
            for (auto const& rec : alignment["vcf_records"])
            {
                if (rec.isMember("id"))
                {
                    if (!event_id.empty())
                    {
                        event_id += ",";
                    }
                    event_id += rec["id"].asString();
                }
            }
            _impl->basic_info["graphinfo"]["ID"] = event_id;
        }

        if (!_impl->basic_info.isMember("breakpointinfo"))
        {
            _impl->basic_info["breakpointinfo"] = Json::arrayValue;

            // write edge + allele map
            const auto& breakpoint_map = _impl->breakpoint_maps.back();
            for (const auto& breakpoint : breakpoint_map)
            {
                Json::Value value = Json::objectValue;

                value["name"] = breakpoint.first;
                value["mapped_alleles"] = Json::objectValue;
                for (const auto& allele : breakpoint.second.allAlleleNames())
                {
                    const auto& canonical_allele = breakpoint.second.getCanonicalAlleleName(allele);
                    if (canonical_allele != allele)
                    {
                        value["mapped_alleles"][allele] = canonical_allele;
                    }
                }

                _impl->basic_info["breakpointinfo"].append(value);
            }
        }

        // load basic json info for output
        vector<string> copied_keys = { "target_regions", "sequencenames" };
        for (auto& key : copied_keys)
        {
            _impl->basic_info["graphinfo"][key] = alignment[key];
        }
        _impl->basic_info["graphinfo"]["nodes"] = Json::arrayValue;
        for (auto const& n : alignment["nodes"])
        {
            Json::Value node = Json::objectValue;
            node["name"] = n["name"];
            if (n.isMember("sequences"))
            {
                node["sequences"] = n["sequences"];
            }
            _impl->basic_info["graphinfo"]["nodes"].append(node);
        }
        _impl->basic_info["graphinfo"]["edges"] = Json::arrayValue;
        for (auto const& e : alignment["edges"])
        {
            Json::Value edge = Json::objectValue;
            edge["name"] = e["from"].asString() + "_" + e["to"].asString();
            if (e.isMember("sequences"))
            {
                edge["sequences"] = e["sequences"];
            }
            _impl->basic_info["graphinfo"]["edges"].append(edge);
        }
    }

    // graph alignment statistics summary
    if (!_impl->basic_info.isMember("samples"))
    {
        _impl->basic_info["samples"] = Json::Value();
    }
    _impl->basic_info["samples"][samplename] = alignment["alignment_statistics"];
    auto& alignment_stat_json = _impl->basic_info["samples"][samplename];
    for (auto& k : alignment["fragment_statistics"].getMemberNames())
    {
        if (k != "linear_histogram" && k != "graph_histogram") // skip histogram output because it is too lengthy
        {
            alignment_stat_json[k] = alignment["fragment_statistics"][k];
        }
    }
}

/**
 * @return set of genotypes for all alignments that were added
 */
Json::Value GraphGenotyper::getGenotypes()
{
    // produce genotypes
    runGenotyping();

    Json::Value result = _impl->basic_info;
    auto& samples = result["samples"];

    for (const auto& samplename : _impl->samplenames)
    {
        samples[samplename]["breakpoints"] = Json::objectValue;
    }

    map<string, GenotypeSet> genotypeSets; // sample->breakpoints to breakpoint->samples. for popluation statistics

    for (size_t isample = 0; isample < _impl->samplenames.size(); ++isample)
    {
        const string& samplename = _impl->samplenames[isample];
        const BreakpointMap& breakpoints = _impl->breakpoint_maps[isample];

        // initialize blank GT
        static const std::vector<std::string> no_alleles;
        static const Genotype empty_genotype = Genotype();

        // print breakpoint genotypes (breakpoint_maps doesn't have "" breakpoint)
        for (const auto& breakpoint : breakpoints)
        {
            const std::string& breakpointname = breakpoint.first;
            auto& this_set = genotypeSets[breakpointname];

            auto gt_it = _impl->graph_genotypes.find(std::make_pair(samplename, breakpointname));
            if (gt_it != _impl->graph_genotypes.end())
            {
                const auto& allele_names = alleleNames();
                this_set.add(allele_names, gt_it->second);
                samples[samplename]["breakpoints"][breakpointname] = Json::objectValue;
                auto& breakpoint_json = samples[samplename]["breakpoints"][breakpointname];
                breakpoint_json["gt"] = gt_it->second.toJson(allele_names);

                // output read counts
                const auto breakpoint_it = _impl->breakpoint_maps[isample].find(breakpointname);
                if (breakpoint_it != _impl->breakpoint_maps[isample].end())
                {
                    breakpoint_json["counts"] = Json::objectValue;
                    breakpoint_json["counts"]["edges"] = Json::objectValue;
                    breakpoint_json["counts"]["alleles"] = Json::objectValue;
                    for (const auto& bp_edgename : breakpoint_it->second.edgeNames())
                    {
                        breakpoint_json["counts"]["edges"][bp_edgename] = breakpoint_it->second.getCount(bp_edgename);
                    }
                    for (const auto& bp_allelename : breakpoint_it->second.canonicalAlleleNames())
                    {
                        breakpoint_json["counts"]["alleles"][bp_allelename]
                            = breakpoint_it->second.getCount(bp_allelename);
                    }
                }
            }
            else
            {
                this_set.add(no_alleles, empty_genotype);
            }
        }

        // print whole variant genotypes
        auto gt_it = _impl->graph_genotypes.find(std::make_pair(samplename, ""));
        auto& this_set = genotypeSets[""];
        if (gt_it != _impl->graph_genotypes.end())
        {
            const auto& allele_names = alleleNames();
            this_set.add(allele_names, gt_it->second);
            samples[samplename]["gt"] = gt_it->second.toJson(allele_names);
        }
        else
        {
            this_set.add(no_alleles, empty_genotype);
            samples[samplename]["gt"] = empty_genotype.toJson(no_alleles);
        }
    }

    // print population statistics for more than one sample
    if (_impl->samplenames.size() > 1)
    {
        result["population"] = Json::objectValue;
        auto& pop = result["population"];
        for (auto& iset : genotypeSets)
        {
            PopulationStatistics ps(iset.second);
            if (iset.first.empty())
            {
                pop = ps.toJson();
            }
            else
            {
                if (!pop.isMember("breakpoints"))
                {
                    pop["breakpoints"] = Json::objectValue;
                }
                pop["breakpoints"][iset.first] = ps.toJson();
            }
        }
    }

    return result;
}

/**
 * @return a list of allele names
 */
std::vector<std::string> const& GraphGenotyper::alleleNames() const { return _impl->allelenames; }

/**
 * Set the genotype for a particular sample
 *
 * @param samplename sample name
 * @param breakpointname name of breakpoint ("" for combined GT)
 * @param alleles names of the alleles for the genotype
 * @param genotype variant genotype(s)
 */
void GraphGenotyper::setGenotype(const std::string& samplename, const std::string& breakpointname, Genotype genotype)
{
    _impl->graph_genotypes[make_pair(samplename, breakpointname)] = std::move(genotype);
}

/**
 * Get the genotype for a particular sample
 *
 * @param samplename sample name
 * @param breakpointname name of breakpoint ("" for combined GT)
 * @return the genotype
 */
Genotype GraphGenotyper::getGenotype(const std::string& samplename, const std::string& breakpointname) const
{
    auto gt_it = _impl->graph_genotypes.find(make_pair(samplename, breakpointname));
    if (gt_it == _impl->graph_genotypes.end())
    {
        return Genotype();
    }
    return gt_it->second;
}

/**
 * @return a list of sample names
 */
std::vector<std::string> const& GraphGenotyper::sampleNames() const { return _impl->samplenames; }

/**
 * @return a list of breakpoint names
 */
std::list<std::string> const& GraphGenotyper::breakpointNames() const { return _impl->breakpointnames; }

/**
 * Get the alignment read counts
 * @param sample_index index of sample (name is in sampleNames[sample_index])
 * @param edge_or_allele_name name of edge or allele
 * @return alignment result for sample
 */
int32_t
GraphGenotyper::getCount(size_t sample_index, string const& breakpoint, std::string const& edge_or_allele_name) const
{
    assert(sample_index < _impl->breakpoint_maps.size());
    auto bp_it = _impl->breakpoint_maps[sample_index].find(breakpoint);
    assert(bp_it != _impl->breakpoint_maps[sample_index].end());
    return bp_it->second.getCount(edge_or_allele_name);
}

/**
 * Get the depth data for a sample
 * @param sample_index index of sample (name is in sampleNames[sample_index])
 * @return pair of expected mean depth and read length
 */
std::pair<double, int> const& GraphGenotyper::getDepthAndReadlength(size_t sample_index) const
{
    return _impl->depths[sample_index];
}

/**
 * Get sex integer for a sample
 */
SampleInfo::Sex GraphGenotyper::getSampleSex(size_t sample_index) const { return _impl->sexes[sample_index]; };
}
