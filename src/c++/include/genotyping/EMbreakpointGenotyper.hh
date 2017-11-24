#pragma once

#include "genotyping/Genotype.hh"
#include "genotyping/Idxdepth.hh"
#include <memory>

namespace genotyping
{

typedef std::vector<int32_t> EdgeCounts;

class EMbreakpointGenotyper
{
public:
    EMbreakpointGenotyper(
        double genotype_error_rate, size_t sample_size, size_t num_haplotypes, int32_t max_read_times,
        int32_t min_overlap_bases);
    ~EMbreakpointGenotyper();

    /**
     * given read count and bam stats, provides final genotypes
     * if use_em, run EM. if not, run fixed
     * return a vector of derived genotypes
     */
    std::vector<Genotype> genotype(std::vector<EdgeCounts>& counts, Idxdepth& idxdepth, bool use_em);

    /**
     * simple getters for EM result
     */
    int numIteration() const;

private:
    /**
     * genotyping functions to call breakpoint genotyper in all samples with given read counts
     */
    void deriveGenotypes(std::vector<EdgeCounts>& counts, Idxdepth& idxdepth, std::map<size_t, bool>& skip_indexes);

    /**
     *  Invode breakpoint genotyper in one sample
     */
    Genotype genotypeSingleSample(EdgeCounts& counts, const SingleIdxdepth& sid);

    /**
     * EM iterations
     */
    void
    updateGenotypesViaEM(std::vector<EdgeCounts>& counts, Idxdepth& idxdepth, std::map<size_t, bool>& skip_indexes);

    /**
     *  update parameters based on this round of genotyping
     */
    void updateParameters(std::vector<EdgeCounts>& counts);

    /**
     * mean squared error of genotype fractions between two iterations
     * M step convergence criteria
     */
    double getGenotypeFractionMse(std::vector<double>& prev_fracs);

    /**
     *  return error rate for each sample in this genotype
     * */
    std::vector<double> getGenotypeErrorVec(std::vector<EdgeCounts>& counts, uint64_t h1, uint64_t h2);

    /**
     *  update haplotype fraction
     */
    double getHaplotypeFraction(std::vector<EdgeCounts>& counts, uint64_t h);

    struct EMgenotyperImpl;
    std::unique_ptr<EMgenotyperImpl> _impl;

    /**
     * filter flags
     */
    static std::string ZEROCOUNT;
    static std::string MAXCOUNT;
};
};
