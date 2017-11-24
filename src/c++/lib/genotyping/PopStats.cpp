#include "genotyping/PopStats.hh"
#include <boost/math/distributions/chi_squared.hpp>

using std::vector;

namespace genotyping
{
PopStats::PopStats(vector<Genotype>& breakpoint_genotypes)
{
    num_valid_samples = 0;
    num_total_samples = static_cast<int>(breakpoint_genotypes.size());
    for (auto& sample_genotype : breakpoint_genotypes)
    {
        if (sample_genotype.gt.empty())
        {
            continue;
        }
        num_valid_samples++;
        if (genotype_counts.find(sample_genotype.gt) == genotype_counts.end())
        {
            genotype_counts[sample_genotype.gt] = 1;
        }
        else
        {
            genotype_counts[sample_genotype.gt]++;
        }
        for (int i = 0; i < 2; i++)
        {
            if (allele_counts.find(sample_genotype.gt[i]) == allele_counts.end())
            {
                allele_counts[sample_genotype.gt[i]] = 1;
            }
            else
            {
                allele_counts[sample_genotype.gt[i]]++;
            }
        }
    }
    // add missingness
    for (auto& al1 : allele_counts)
    {
        for (auto& al2 : allele_counts)
        {
            if (al1.first > al2.first)
            {
                continue;
            }
            if (genotype_counts.find({ al1.first, al2.first }) == genotype_counts.end())
            {
                genotype_counts[{ al1.first, al2.first }] = 0;
            }
        }
    }
}

PopStats::PopStats(vector<VariantGenotype>& variant_genotypes)
{
    num_valid_samples = 0;
    num_total_samples = static_cast<int>(variant_genotypes.size());
    for (auto& sample_genotype : variant_genotypes)
    {
        if (sample_genotype.empty())
        {
            continue;
        }
        num_valid_samples++;
        GenotypeVector gt = sample_genotype.getGenotype();
        if (genotype_counts.find(gt) == genotype_counts.end())
        {
            genotype_counts[gt] = 1;
        }
        else
        {
            genotype_counts[gt]++;
        }
        for (int i = 0; i < 2; i++)
        {
            if (allele_counts.find(gt[i]) == allele_counts.end())
            {
                allele_counts[gt[i]] = 1;
            }
            else
            {
                allele_counts[gt[i]]++;
            }
        }
    }
    // add missingness
    for (auto& al1 : allele_counts)
    {
        for (auto& al2 : allele_counts)
        {
            if (al1.first > al2.first)
            {
                continue;
            }
            if (genotype_counts.find({ al1.first, al2.first }) == genotype_counts.end())
            {
                genotype_counts[{ al1.first, al2.first }] = 0;
            }
        }
    }
}

double PopStats::getHWE()
{
    if (genotype_counts.size() <= 1)
    {
        return 1;
    }

    double chisq_val = 0;
    std::map<GenotypeVector, int> expected_counts;
    for (auto& gv : genotype_counts)
    {
        uint64_t h1 = gv.first[0];
        uint64_t h2 = gv.first[1];
        int e_count;
        if (h1 == h2)
        {
            e_count = std::round(
                ((double)allele_counts[h1] / num_valid_samples / 2)
                * ((double)allele_counts[h1] / num_valid_samples / 2) * num_valid_samples);
        }
        else
        {
            e_count = std::round(
                2 * ((double)allele_counts[h1] / num_valid_samples / 2)
                * ((double)allele_counts[h2] / num_valid_samples / 2) * num_valid_samples);
        }

        int diff = e_count - gv.second;
        // the case where expected is zero
        // this should be replaced with fisher's exact test for rare alleles
        if (e_count == 0)
        {
            e_count = 1;
        }
        double norm_diff_square = (double)(diff * diff) / e_count;
        chisq_val += norm_diff_square;
    }
    boost::math::chi_squared chisq_distribution(1);
    double hwe_pval = 1 - boost::math::cdf(chisq_distribution, chisq_val);
    return hwe_pval;
}
}