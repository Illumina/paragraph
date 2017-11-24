#include "genotyping/BreakpointStat.hh"
#include "genotyping/PopStats.hh"
#include <algorithm>

using std::string;
using std::vector;

namespace genotyping
{
BreakpointStat::BreakpointStat(
    string& node_name, bool from_node_fixed, vector<string>& neighbor_names,
    std::map<string, vector<uint64_t>>& edge_to_seq_indexes)
{
    for (auto& neighbor_name : neighbor_names)
    {
        string edge_name = from_node_fixed ? (node_name + "_" + neighbor_name) : (neighbor_name + "_" + node_name);
        edge_names.push_back(edge_name);
        haplotype_indexes.push_back(edge_to_seq_indexes[edge_name]);
        std::map<uint64_t, bool> dup_index_check_hash; // check for duplicate sequence indexes from one node
        for (auto& index : edge_to_seq_indexes[edge_name])
        {
            if (dup_index_check_hash.find(index) != dup_index_check_hash.end())
            {
                throw std::logic_error("More than one edge from the same node have a same sequence name. "
                                       "Please double check your graph description!");
            }
            dup_index_check_hash[index] = true;
        }
    }
}

void BreakpointStat::addBreakpointEdgeCountForOneSample(Json::Value& paragraph_json, const string& sample_name)
{
    if (!paragraph_json["samples"].isMember(sample_name))
    {
        edge_counts.push_back({});
        return;
    }
    Json::Value* sample_json = &paragraph_json["samples"][sample_name];
    if (!sample_json->isMember("read_counts_by_edge"))
    {
        string message = "cannot find key read_counts_by_edge in sample " + sample_name;
        throw std::logic_error(message);
    }
    vector<int32_t> current_sample_edge_counts;
    for (auto& edge_name : edge_names)
    {
        int single_edge_count = (*sample_json)["read_counts_by_edge"].isMember(edge_name)
            ? (*sample_json)["read_counts_by_edge"][edge_name].asInt()
            : 0;
        current_sample_edge_counts.push_back(static_cast<int32_t>(single_edge_count));
    }
    edge_counts.push_back(current_sample_edge_counts);
}

void BreakpointStat::genotype(
    double genotype_error_rate, int sample_size, int max_read_times, int min_overlap_bases, Idxdepth& idxdepth,
    bool use_em)
{
    EMbreakpointGenotyper em_genotyper(
        genotype_error_rate, static_cast<size_t>(sample_size), edge_counts[0].size(),
        static_cast<int32_t>(max_read_times), static_cast<int32_t>(min_overlap_bases));
    genotypes = em_genotyper.genotype(edge_counts, idxdepth, use_em);
    recodeGenotypesWithSequenceIndex();

    num_iterations = em_genotyper.numIteration();
    PopStats pop_stats(genotypes);
    call_rate = pop_stats.callRate();
    hwe_p = pop_stats.getHWE();
}

void BreakpointStat::recodeGenotypesWithSequenceIndex()
{
    for (auto& raw_genotype : genotypes)
    {
        if (!raw_genotype.gt.empty())
        {
            raw_genotype.recode(haplotype_indexes);
        }
    }
}

Json::Value BreakpointStat::descriptionsToJson()
{
    Json::Value output_json;
    if (genotypes.size() > 1)
    {
        if (num_iterations >= 0)
        {
            output_json["num_iterations"] = num_iterations;
        }
        output_json["call_rate"] = call_rate;
        output_json["hwe_p"] = hwe_p;
    }
    return output_json;
}

Json::Value BreakpointStat::sampleGenotypeToJson(int sample, std::vector<std::string>& allele_names)
{
    return genotypes[sample].toJson(allele_names);
}

int32_t BreakpointStat::getEdgeCount(int sample, std::string& edge_name)
{
    int edge_count;
    auto p = std::find(edge_names.begin(), edge_names.end(), edge_name);
    if (p != edge_names.end())
    {
        if (!edge_counts[sample].empty())
        {
            int index = p - edge_names.begin();
            edge_count = edge_counts[sample][index];
        }
        else
        {
            edge_count = 0;
        }
    }
    else
    {
        edge_count = 0;
    }
    return edge_count;
}
}
