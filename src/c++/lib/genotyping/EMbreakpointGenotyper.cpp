#include "genotyping/EMbreakpointGenotyper.hh"
#include "genotyping/BreakpointGenotyper.hh"
#include <cmath>
#include <limits>
#include <numeric>

using std::string;
using std::vector;

namespace genotyping
{

string EMbreakpointGenotyper::ZEROCOUNT = "ZEROCOUNT";
string EMbreakpointGenotyper::MAXCOUNT = "MAXCOUNT";

struct EMbreakpointGenotyper::EMgenotyperImpl
{
    EMgenotyperImpl(
        double genotype_error_rate, size_t sample_size, size_t num_haplotypes, int32_t max_read_times,
        int32_t min_overlap_bases)
        : starting_error_rate(genotype_error_rate)
        , sample_size(sample_size)
        , num_haplotypes(num_haplotypes)
        , max_read_times(max_read_times)
        , min_overlap_bases(min_overlap_bases)
        , num_genotypes(num_haplotypes * (num_haplotypes + 1) / 2)
        , default_prior((double)1 / num_genotypes)
        , min_difference(0.001)
        , max_num_iterations(100)
        , reference_haplotype_index(0)
        , reference_haplotype_fraction(0.5)
        , default_hap_fraction(0.45){};

    std::vector<Genotype> genotypes;

    /**
     *  quick access to genotype index given 2 haplotype indexes
     */
    std::map<GenotypeVector, int> genotype_index_map;

    /**
     *  iteratively updated parameters
     */
    std::vector<double> error_rates;
    std::vector<double> haplotype_fractions;
    std::vector<double> genotype_fractions;

    /**
     *  updated stats from EM result
     */
    bool is_converged;
    int num_valid_samples;
    int num_iterations_to_converge;

    /**
     *  const stats for EM and fixed genotyping
     */
    const double starting_error_rate;
    const size_t sample_size;
    const size_t num_haplotypes;
    const int32_t max_read_times;
    const int32_t min_overlap_bases;
    const size_t num_genotypes;
    const double default_prior;
    const double min_difference;
    const int max_num_iterations;
    const uint64_t reference_haplotype_index;
    const double reference_haplotype_fraction;
    const double default_hap_fraction;
};

EMbreakpointGenotyper::EMbreakpointGenotyper(
    double genotype_error_rate, size_t sample_size, size_t num_haplotypes, int32_t max_read_times,
    int32_t min_overlap_bases)
    : _impl(new EMgenotyperImpl(genotype_error_rate, sample_size, num_haplotypes, max_read_times, min_overlap_bases))
{
}

EMbreakpointGenotyper::~EMbreakpointGenotyper() = default;

int EMbreakpointGenotyper::numIteration() const { return _impl->num_iterations_to_converge; }

vector<Genotype> EMbreakpointGenotyper::genotype(vector<EdgeCounts>& counts, Idxdepth& idxdepth, bool use_em)
{
    // stats initialize
    int genotype_index = 0;
    for (uint64_t i = 0; i < _impl->num_haplotypes; i++)
    {
        for (uint64_t j = i; j < _impl->num_haplotypes; j++)
        {
            _impl->genotype_index_map[{ i, j }] = genotype_index;
        }
        genotype_index++;
    }
    _impl->error_rates.resize(_impl->num_genotypes, _impl->starting_error_rate);
    _impl->haplotype_fractions.resize(_impl->num_haplotypes, _impl->reference_haplotype_fraction);
    _impl->genotype_fractions.resize(_impl->num_genotypes, (double)1 / _impl->num_genotypes);

    _impl->genotypes.reserve(_impl->sample_size);
    for (size_t i = 0; i < _impl->sample_size; i++)
    {
        _impl->genotypes.push_back(genotyping::Genotype());
        _impl->genotypes[i].gt = {};
    }

    // label missing genotypes
    std::map<size_t, bool> skip_indexes;
    for (size_t index = 0; index < _impl->sample_size; index++)
    {
        if (counts[index].empty())
        {
            skip_indexes[index] = true;
            _impl->genotypes[index] = {};
            _impl->genotypes[index].filter = ZEROCOUNT;
            continue;
        }
        int total_count = std::accumulate(counts[index].begin(), counts[index].end(), 0);
        if (total_count == 0)
        {
            skip_indexes[index] = true;
            _impl->genotypes[index] = {};
            _impl->genotypes[index].filter = ZEROCOUNT;
            continue;
        }
        SingleIdxdepth sid = idxdepth.getIdxStats(index);
        auto max_edge_count = static_cast<const int>(_impl->max_read_times * sid.autosome_depth);
        if (total_count > max_edge_count)
        {
            skip_indexes[index] = true;
            _impl->genotypes[index] = {};
            _impl->genotypes[index].filter = MAXCOUNT;
            continue;
        }
    }
    _impl->num_valid_samples = _impl->sample_size - static_cast<int>(skip_indexes.size());

    // check if all missing
    if (skip_indexes.size() == _impl->sample_size)
    {
        return _impl->genotypes;
    }

    // initialize genotypes using fixed parameters
    deriveGenotypes(counts, idxdepth, skip_indexes);

    // runEM
    if (use_em)
    {
        updateGenotypesViaEM(counts, idxdepth, skip_indexes);
    }
    else
    {
        _impl->is_converged = false;
        _impl->num_iterations_to_converge = -1;
    }

    return _impl->genotypes;
}

void EMbreakpointGenotyper::deriveGenotypes(
    vector<EdgeCounts>& counts, Idxdepth& idxdepth, std::map<size_t, bool>& skip_indexes)
{
    for (size_t i = 0; i < idxdepth.sampleSize(); i++)
    {
        if (skip_indexes.find(i) != skip_indexes.end())
        {
            continue;
        }
        _impl->genotypes[i] = genotypeSingleSample(counts[i], idxdepth.getIdxStats(i));
    }
}

Genotype EMbreakpointGenotyper::genotypeSingleSample(EdgeCounts& counts, const SingleIdxdepth& sid)
{
    BreakpointGenotyper genotyper(
        sid.autosome_depth, sid.read_length, _impl->starting_error_rate, _impl->min_overlap_bases);
    Genotype genotype = genotyper.genotype(counts, _impl->haplotype_fractions, _impl->error_rates);
    return genotype;
}

void EMbreakpointGenotyper::updateGenotypesViaEM(
    vector<EdgeCounts>& counts, Idxdepth& idxdepth, std::map<size_t, bool>& skip_indexes)
{
    // update parameters first
    vector<Genotype> prev_genotypes;

    // breakpoint genotyping through loops
    int num_iterations = 0;
    while (num_iterations < _impl->max_num_iterations)
    {
        vector<double> prev_genotype_fractions = _impl->genotype_fractions;
        updateParameters(counts);
        deriveGenotypes(counts, idxdepth, skip_indexes);
        double difference_this_loop = getGenotypeFractionMse(prev_genotype_fractions);
        if (difference_this_loop <= _impl->min_difference)
        {
            _impl->is_converged = true;
            _impl->num_iterations_to_converge = num_iterations;
            break;
        }
        num_iterations++;
    }
}

double EMbreakpointGenotyper::getGenotypeFractionMse(vector<double>& prev_fracs)
{
    double mse = 0;
    for (size_t i = 0; i < _impl->num_genotypes; i++)
    {
        double diff = prev_fracs[i] - _impl->genotype_fractions[i];
        mse += diff * diff;
    }
    mse /= _impl->num_genotypes;
    return mse;
}

void EMbreakpointGenotyper::updateParameters(vector<EdgeCounts>& counts)
{
    // update error rate and _impl->genotype_fractions
    for (uint64_t h1 = 0; h1 < _impl->num_haplotypes; h1++)
    {
        for (uint64_t h2 = h1; h2 < _impl->num_haplotypes; h2++)
        {
            int genotype_index = _impl->genotype_index_map[{ h1, h2 }];
            vector<double> error_rate_vec = getGenotypeErrorVec(counts, h1, h2);
            if (error_rate_vec.empty())
            {
                _impl->error_rates[genotype_index] = _impl->starting_error_rate;
                _impl->genotype_fractions[genotype_index] = _impl->default_prior;
            }
            else
            {
                _impl->error_rates[genotype_index] = std::accumulate(error_rate_vec.begin(), error_rate_vec.end(), 0.0)
                    / (double)error_rate_vec.size();
                _impl->genotype_fractions[genotype_index]
                    = (double)_impl->error_rates.size() / _impl->num_valid_samples;
            }
        }
    }
    // update hap fraction
    for (uint64_t h = 0; h < _impl->num_haplotypes; h++)
    {
        _impl->haplotype_fractions[h] = getHaplotypeFraction(counts, h);
    }
}

vector<double> EMbreakpointGenotyper::getGenotypeErrorVec(vector<EdgeCounts>& counts, uint64_t h1, uint64_t h2)
{
    vector<double> error_rate_vec;
    GenotypeVector target_gt = { h1, h2 };
    for (size_t s = 0; s < _impl->sample_size; s++)
    {
        if (_impl->genotypes[s].gt != target_gt)
        {
            continue;
        }
        int num_spurious = 0;
        for (size_t i = 0; i < _impl->num_haplotypes; i++)
        {
            if (i == h1 || i == h2)
            {
                continue;
            }
            num_spurious += counts[s][i];
        }
        error_rate_vec.push_back((double)num_spurious / _impl->genotypes[s].num_reads);
    }
    return error_rate_vec;
}

double EMbreakpointGenotyper::getHaplotypeFraction(vector<EdgeCounts>& counts, uint64_t h)
{
    if (h == _impl->reference_haplotype_index)
    {
        return _impl->reference_haplotype_fraction;
    }
    GenotypeVector target_gt
        = { std::min(h, _impl->reference_haplotype_index), std::max(h, _impl->reference_haplotype_index) };
    vector<double> allele_fractions;
    for (size_t s = 0; s < _impl->sample_size; s++)
    {
        if (_impl->genotypes[s].gt == target_gt)
        {
            allele_fractions.push_back(counts[s][h] / (double)_impl->genotypes[s].num_reads);
        }
    }
    double hap_fraction;
    if (allele_fractions.empty())
    {
        hap_fraction = _impl->default_hap_fraction;
    }
    else
    {
        hap_fraction
            = std::accumulate(allele_fractions.begin(), allele_fractions.end(), 0.0) / (double)allele_fractions.size();
    }
    return hap_fraction;
}
}