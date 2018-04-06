# Filters used in genotyper output

## Variant level filters

* **PASS**

Variant PASS all filters

* **CONFLICT**

Variant has genotype conflicts in one or more breakpoints

* **EXIST_BAD_BP**

Varaint has one or more breakpoint that fails breakpoint-level filter

* **ALL_BAD_BP**

All breakpoints in this variant fail breakpoint-level filter

* **MISSING**

Variant has one or more breakpoints with no spanning read

## Breakpoint level filters

* **DEPTH**

Total number of reads on this breakpoint (from all alleles) fail the coverage test