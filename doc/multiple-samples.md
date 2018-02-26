# Genotyping large set of variants across multiple samples

A simple way to parallelize genotyping of a large variant set across many
samples is to use [Snakemake](https://snakemake.readthedocs.io/en/stable/).
For example, suppose that the working directory `analysis/` has the following
contents (see also directory [pg-example](pg-example)).

```bash
analysis
├── manifests
│   ├── NA12877.txt
│   └── NA12878.txt
├── Snakemake
└── variants.vcf
```

where directory `manifests` contains the manifest files (one per sample),
`variants.vcf` contains a list of variants, and the file `Snakefile` is as
follows.

```python
import os

reference='/path/to/reference/genome.fa'
vcf='variants.vcf'

# Collect sample names.
samples = [fname.replace('.txt', '') for fname in os.listdir('manifests/')]

rule targets:
    input:
        expand('results/{sample}/genotype.json.gz', sample=samples)

rule graph_typing:
    input:
        'manifests/{sample}.txt'
    output:
        'results/{sample}/genotypes.json.gz'
    params:
        ref=reference,
        vcf=vcf,
        sample='{sample}'
    shell:
        r'''
            export PARAGRAPH=/path/to/paragraph/installation
            python3 ${{PARAGRAPH}}/multigrmpy.py \
                -i {params.vcf} \
                -m {input} \
                -r {params.ref} \
                -o results/{params.sample}
        '''
```

With these files in place, the analysis can be started by executing `snakemake` command with the appropriate parameters specifying the number of concurrent jobs, computing cluster configuration (if appropriate), etc. Please see [Snakemake](https://snakemake.readthedocs.io/en/stable/) documentation for the details.