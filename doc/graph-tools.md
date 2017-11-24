
## <a name='Furtherinformationonexecutablesinbin'></a>Further information on executables in bin/

<!-- vscode-markdown-toc -->
* [Further information on executables in bin/](#Furtherinformationonexecutablesinbin)
* [Bam stats estimation](#Bamstatsestimation)
	* [idxdepth](#idxdepth)
* [Read counts on graph](#Readcountsongraph)
	* [multiparagraph.py](#multiparagraph.py)
	* [paragraph](#paragraph)
* [Genotyper](#Genotyper)
	* [multigrmpy.py](#multigrmpy.py)
	* [grmpy](#grmpy)
* [Other tools](#Othertools)
	* [vcf2paragraph.py](#vcf2paragraph.py)
	* [compare_alignments.py](#compare_alignments.py)
	* [findgrm.py](#findgrm.py)
	* [msa2vcf.py](#msa2vcf.py)
	* [paragraph2dot.py](#paragraph2dot.py)

<!-- vscode-markdown-toc-config
	numbering=false
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

## <a name='Bamstatsestimation'></a>Bam stats estimation

### <a name='idxdepth'></a>idxdepth

A brief index-based depth estimation for BAM and CRAM.

Basic usage:

```bash
bin/idxdepth -b \<bam/cram> -r \<reference fasta> -o \<output>
```

The output will be a json file with depth across whole genome and on each chromosome:
```javascript
{
    "autosome": {
        "contigs": [
            "chr1"
        ],
        "depth": 1
    },
    "bam_path": "fake_path.bam",
    "contigs": [
        // ...
    ],
    "read_length": 50,
    "reference": "fake_reference.fa",
    "unaligned_reads": 0
}
```

## <a name='Readcountsongraph'></a>Read counts on graph

### <a name='multiparagraph.py'></a>multiparagraph.py

### <a name='paragraph'></a>paragraph

core program used in grmpy to calculate read counts on sequence graphs.
The main idea in this tool is to support read-based calling and disambiguation
of reads between different paths in a graph.

The input for `paragraph` are a graph model JSON file and a BAM file + reference
fasta to retrieve reads from.

The output file contains:

* read counts for every node and edge in the graph (for uniquely-aligning reads)
* locations and read counts for variants (i.e. non-match alignments)
* (optional) individual read alignments for every read

## <a name='Genotyper'></a>Genotyper

### <a name='multigrmpy.py'></a>multigrmpy.py

Wrapper for the graph tool genotyper. Transform multiparagraph output to grmpy (graph genotyper) input format then run grmpy.

Basic Usage:

```bash
python3 bin/multigrmpy.py -m \<manifest> -r \<reference fasta> -o \<output directory>
```

*  **manifest**: list of all paragraph output files.
```
    sample1 test/paragraph/sample1/raw_pg.json      1       50
    sample2 test/paragraph/sample2/raw_pg.json      1       50
```

### <a name='grmpy'></a>grmpy

Core executable for graph genotyper. Process one site at a time. Executed within multigrmpy.py

Paragraph output cannot be directly used for grmpy unless file format conversion has been done manually or through multigrmpy.py.

An example input for grmpy:

```javascript
    "model_name": "Variant1"
    // ... original graph specification
    "samples": {
        "sample1": {
            "read_counts_by_edge": {
                // ...
            }
        },
        "sample2": {
            // ...
        }
    },
    // ...
```

## <a name='Othertools'></a>Other tools

### <a name='vcf2paragraph.py'></a>vcf2paragraph.py

VCF2Paragraph can convert a subset of the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) into a labelled sequence graph.

It supports two different conversion modes:

* **Allele graphs**: each VCF record creates a reference and alternate path. When two alleles overlap on the reference, their paths will not be compatible in the graph. When two alleles are subsequent on the reference, their respective reference / alternate paths can be combined. Sequence labels for the graph are generated automatically using a greedy algorithm which makes sure each edge has a sequence label.
* **Haplotype graphs**: Each VCF sample may specify one haplotype. One sequence / set of contiguous paths are generated for each sample, which must contain only haploid genotypes specifying the allele
to be included in the sequence path. Reference paths may be included using `0` genotypes.

The script supports standard sequence alleles (REF/ALT), as well a subset of symbolic ALT alleles:

*  Sequence alleles: SNPs and small indels / substitutions produce subgraphs with two paths: one for the reference, and one for the alternate sequence.
*  Deletions: explicit deletions or symbolic `<DEL>` alleles segment the reference sequence into three parts. The middle part can be bypassed using an edge. The `END` INFO field can be used to indicate the length of the deletion for long deletions.
* Forward breakends on the same chromosome, e.g. REF: `A`; ALT: `AC[chr1:10000[`, can be used to encode
  long deletions with short inserted sequence (swaps / substitutions).

All variants in one VCF passed to VCF2paragraph.py must be on the same chromosome. Our
other wrapper scripts (e.g. `runGraphTyping.py`) will split up the VCF file into individual records 
before conversion.

Paragraph doesn't support cycles in the graphs, and VCF2paragraph will raise an error if the resulting
graph contains a cycle (e.g. created by a breakend). To break cycles, vcf2paragraph may split long reference nodes and remove long reference stretches (our graph alignments start from BAM files, 
so realignments to long reference stretches are not necessary).

Here is a full list of command line options:

```
usage: vcf2paragraph.py [-h] -r REF [-g {alleles,haplotypes}] [-R]
                        [-l MAX_REF_LEN] [-p READ_LEN] [-T TARGET_REGIONS]
                        [--ref-path] [--crossovers]
                        input output
vcf2paragraph.py: error: the following arguments are required: input, output, -r/--reference-sequence
grm.py-build $ python3 bin/vcf2paragraph.py -h
usage: vcf2paragraph.py [-h] -r REF [-g {alleles,haplotypes}] [-R]
                        [-l MAX_REF_LEN] [-p READ_LEN] [-T TARGET_REGIONS]
                        [--ref-path] [--crossovers]
                        input output

positional arguments:
  input                 Input VCF / BCF file
  output                Output JSON file

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --reference-sequence REF
                        Reference FASTA for checking REF and resolving <DEL>

Common VCF graph options:
  -g {alleles,haplotypes}, --graph-type {alleles,haplotypes}
                        Select the type of graph to generate
  -R, --retrieve-reference-sequence
                        Retrieve reference sequence for REF nodes
  -l MAX_REF_LEN, --max-ref-node-length MAX_REF_LEN
                        Maximum length of reference nodes before they get
                        padded and truncated.
  -p READ_LEN, --read-length READ_LEN
                        Read length -- this can be used to add reference
                        padding for disambiguation.
  -T TARGET_REGIONS, --target-region TARGET_REGIONS
                        Target regions for read retrieval; also, reference
                        nodes will not be clipped inside target regions.

Haploid VCF graph options:
  --ref-path            Add edges for the reference path.
  --crossovers          Add crossover edges (connect directly adjacent nodes
                        even when no sample gives these combinations).
```


### <a name='compare_alignments.py'></a>compare_alignments.py

### <a name='findgrm.py'></a>findgrm.py

### <a name='msa2vcf.py'></a>msa2vcf.py

### <a name='paragraph2dot.py'></a>paragraph2dot.py


