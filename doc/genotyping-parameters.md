# How to set genotyping paramters externally

This file shows how to set customized genotyping parameters externally in a JSON file.

External genotyping parameter can be passed to multigrmpy.py and grmpy as JSON file with command line option `-G`

Below includes all allowed parameter fields:

```javascript
{
     // ploidy of samples
    "ploidy": 2

    // Minimal overlapping bases between a read and a mapped edge.
    //      Used in estimating distribution parameter in genotyper
    "min_overlap_bases": 16,


    // Minimal depth p value for a breakpoint to pass depth test
    "coverage_test_cutoff": 0.00001

    // Allele names in graph(s).
    //     If other alleles were observed in graph, they will be excluded from analysis.
    "allele_names": [
        "REF",
        "ALT"
    ],


    // Error rate for each allele.
    //     Alleles must be in the same orrder as "allele_names"
    "allele_error_rates": [
        0.04,
        0.04
    ],
    // These two fields will be effective only if "allele_error_rates" doesn't exist.
    "reference_allele_error_rate": 0.05, //  error rate of reference allele
    "other_allele_error_rate": 0.05, // uniformed error rate of all alternative alleles


    // Fraction of the non-reference allele in a reference-alternative heterozygotes.
    //     Alleles must be in the same order as "allele_names"
    "het_haplotype_fractions": [
        0.5,
        0.45
    ],
    // het_haplotype_fraction can be a float to represent a uniform fraction, such as:
    // "het_haplotype_fractions": 0.5




    // Fraction of each genotype in the population.
    //     Genotype indices are zero-based, relative to the "allele_names" array
    "genotype_fractions": {
        "0/0": 0.5,
        "0/1": 0.3,
        "1/1": 0.2
    },
    // This field will be effective only if "genotype_fractions" doesn't exist
    "other_genotype_fraction": 0.33, // uniform genotype fraction for each genotype
}
```
