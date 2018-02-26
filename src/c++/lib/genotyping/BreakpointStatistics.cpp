#include "genotyping/BreakpointStatistics.hh"
#include <algorithm>

#include "common/Error.hh"

using std::string;
using std::vector;

namespace genotyping
{

/**
 * Create breakpoint statistics on graph
 * @param wgraph graph information
 * @param node_name name of node
 * @param forward true to use successor node/edges, false to use predecessors
 */
BreakpointStatistics::BreakpointStatistics(
    graphs::WalkableGraph const& wgraph, std::string const& node_name, bool forward)
{
    const auto node_id = wgraph.nodeId(node_name);
    const auto allele_nodes = forward ? wgraph.succ(node_id) : wgraph.pred(node_id);

    // note that sequences / sequence labels that have the same name as an edge
    // will not work well with this.
    assert(allele_nodes.size() > 1);
    for (auto const& an : allele_nodes)
    {
        const auto an_name = wgraph.nodeName(an);
        const std::string edge_name = forward ? (node_name + "_" + an_name) : (an_name + "_" + node_name);
        const auto edge = forward ? wgraph.edge(node_id, an) : wgraph.edge(an, node_id);
        for (const auto& sequence : edge->sequence_ids())
        {
            const auto allele_name = wgraph.header().sequencenames((int)sequence);
            auto a_it = std::find(allele_names.begin(), allele_names.end(), allele_name);
            if (a_it == allele_names.end())
            {
                allele_names.push_back(allele_name);
                allele_name_to_index[allele_name] = allele_names.size() - 1;
                edgename_to_alleles[edge_name].push_back(allele_names.size() - 1);
            }
            else
            {
                edgename_to_alleles[edge_name].push_back(static_cast<unsigned long&&>(a_it - allele_names.begin()));
            }
        }
    }
    for (auto const& edgename : edgename_to_alleles)
    {
        edge_names.push_back(edgename.first);
        edge_name_to_index[edgename.first] = edge_names.size() - 1;
    }

    // put REF allele first
    auto ref_it = allele_name_to_index.find("REF");
    if (ref_it != allele_name_to_index.end() && ref_it->second != 0)
    {
        const auto allele_to_swap_name = allele_names[0];
        const size_t allele_to_swap_index = ref_it->second;
        // swap names
        allele_names[0] = "REF";
        allele_names[ref_it->second] = allele_to_swap_name;
        // swap index positions
        allele_name_to_index[allele_to_swap_name] = ref_it->second;
        allele_name_to_index["REF"] = 0;
        // change the mapping edgename->allele_index
        for (auto& edgename : edgename_to_alleles)
        {
            for (auto& original_allele : edgename.second)
            {
                if (original_allele == 0)
                {
                    original_allele = allele_to_swap_index;
                }
                else if (original_allele == allele_to_swap_index)
                {
                    original_allele = 0;
                }
            }
        }
    }
}

/**
 * Add edge counts from paragraph output
 * @param paragraph_json JSON output from alignAndDisambiguate
 */
void BreakpointStatistics::addCounts(Json::Value const& paragraph_json)
{
    if (!paragraph_json.isMember("read_counts_by_edge"))
    {
        error("Cannot find key read_counts_by_edge in JSON");
    }
    for (auto const& edge_name : edge_names)
    {
        const auto e_it = edge_name_to_index.find(edge_name);
        assert(e_it != edge_name_to_index.end());
        const int this_edge_count = paragraph_json["read_counts_by_edge"].isMember(edge_name)
            ? paragraph_json["read_counts_by_edge"][edge_name].asInt()
            : 0;

        if (this_edge_count == 0)
        {
            continue;
        }

        if (edge_counts.size() <= e_it->second)
        {
            edge_counts.resize(edge_names.size(), 0);
            assert(edge_counts.size() > e_it->second);
        }

        edge_counts[e_it->second] += this_edge_count;

        for (const auto& allele : edgename_to_alleles[edge_name])
        {
            if (allele_counts.size() <= allele)
            {
                allele_counts.resize(allele_names.size(), 0);
                assert(allele_counts.size() > allele);
            }
            allele_counts[allele] += this_edge_count;
        }
    }
}

int32_t BreakpointStatistics::getCount(std::string const& edge_or_allele_name) const
{
    const auto e_it = edge_name_to_index.find(edge_or_allele_name);
    const auto a_it = allele_name_to_index.find(edge_or_allele_name);
    int32_t result = 0;
    if (e_it != edge_name_to_index.end() && a_it != allele_name_to_index.end())
    {
        error("Allele / sequence name %s is ambiguous with an edge name.", edge_or_allele_name.c_str());
    }
    else if (e_it != edge_name_to_index.end())
    {
        result = e_it->second >= edge_counts.size() ? 0 : edge_counts[e_it->second];
    }
    else if (a_it != allele_name_to_index.end())
    {
        result = a_it->second >= allele_counts.size() ? 0 : allele_counts[a_it->second];
    }
    else
    {
        error("Unknown edge or allele: %s", edge_or_allele_name.c_str());
    }
    return result;
}
}
