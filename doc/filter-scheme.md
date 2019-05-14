# Filters used in genotyper output

## Breakpoint level filters

* **GQ**

Low genotype quality for this breakpoint

* **NO_READS**

No reads in this breakpoint

* **BP_DEPTH**

Total number of reads on this breakpoint (from all alleles) fail the coverage test

## Variant level filters

All breakpoint level filters can be included in variant level filters as well.

* **PASS**

Variant PASS all filters

* **CONFLICT**

Variant has genotype conflicts in one or more breakpoints

* **BP_NO_GT**

Exist one or more breakpoint with missing genotypes

* **NO_VALID_GT**

All breakpoints have missing genotypes