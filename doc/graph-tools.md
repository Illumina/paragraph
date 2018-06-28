
## <a name='Furtherinformationonexecutablesinbin'></a>Further information on executables in bin/

<!-- vscode-markdown-toc -->
* [Further information on executables in bin/](#Furtherinformationonexecutablesinbin)
* [Bam stats estimation](#Bamstatsestimation)
	* [idxdepth](#idxdepth)
* [Read counts on graph](#Readcountsongraph)
	* [paragraph](#paragraph)
* [Genotyper](#Genotyper)
	* [multigrmpy.py](#multigrmpy.py)
	* [grmpy](#grmpy)
* [Other tools](#Othertools)
	* [vcf2paragraph.py](#vcf2paragraph.py)
	* [addVariants.py](#addVariants.py)
	* [compare-alignments.py](#compare-alignments.py)
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
Idxdepth will rapidly provide depth estimates for all chromosomes in a
BAM / CRAM file. When passed a human genome, idxdepth will also guess which
chromosomes are in the autosome and provide an autosome depth estimate and
a guess for the read length.

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
python3 bin/multigrmpy.py -i \<input\> -m \<manifest> -r \<reference fasta> -o \<output directory>
```

* **input**: VCF file with variants or JSON file with graphs

* **manifest**: a list of BAM files or alignments to use.

Example using BAM files and known depth + read length:
```
id      path          read length  depth
sample1 sample1.bam   150          50
sample2 sample2.bam   150          50
```

Idxdepth JSON files may be used instead of read length and depth estimates:
```
id      path          idxdepth
sample1 sample1.bam   sample1.idxdepth.json
sample2 sample2.bam   sample2.idxdepth.json
```

When paragraph alignments have been computed beforehand, these may be passed instead of a BAM file:
```
id      paragraph      idxdepth
sample1 sample1.paragraph.json   sample1.idxdepth.json
sample2 sample2.paragraph.json   sample2.idxdepth.json
```

If the BAM index isn't in a standard location that can be derived from a file name (e.g.
when the BAM is located on remote storage and accessed via a pre-signed URL that may be
different from the one of the index), the `index_path` column can be added to specify the
index location for each BAM file.

* **reference fasta**: Fasta reference to use. This must be the correct reference both for input BAM/CRAM files and the input graphs.

* **output directory**: This directory will contain the following files:
    - *variants.json.gz*: A JSON file with all input variants.
    - *genotypes.json.gz*: A JSON file with breakpoint and genotype information.
    - *genotypes.vcf.gz*: When the input file is in VCF format, then multigrmpy.py will write an annotated version which contains genotypes as computed by grmpy. In addition to the genotypes, this file also contains the following new fields:
        + `INFO/GRMPY_ID`: This value can be used to match the genotypes to the JSON records in *genotypes.json.gz* (where it matches record.graphinfo.ID) and *variants.json.gz* (where it matches record.ID).
        + `FORMAT/GT`: genotype computed by grmpy. Note that the `SAMPLE` column in the VCF file must match a sample id from the manifest shown above.
        + `FORMAT/AD`, `FORMAT/ADF`, `FORMAT/ADR`: depth for each allele (including the reference), total and by strand.
        + `FORMAT/DP`: Total depth used to genotype.
        + `FORMAT/FT`: grmpy filter status for each call, in each sample.
        + `FORMAT/PL`: phred-scaled genotype likelihood.
    - *grmpy.log*: log file from grmpy with progress and diagnostic information.
    - *alignments/\*.json.gz*: when multigrmpy.py is run with the `-A` option, this folder will contain a paragraph alignment for each graph.

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

* **Allele graphs**: each VCF record creates a path for the reference and each alternate allele. When two alleles overlap on the reference, their paths will not be compatible in the graph. When two alleles are subsequent on the reference, their respective reference / alternate paths can be combined. Sequence labels for each allele in the form *vcfID*:*AlleleNumber*, e.g. *var1:0* for the reference allele of a VCF variant with ID *var1*.
* **Haplotype graphs**: Each VCF sample is used to create one haplotype path by specifying haploid genotypes. Reference paths may be included using `0` genotypes, variants can be undetermined on a haplotype (allowing any allele) using `.`. Alt alleles that are not used in any haplotype are excluded from the graph. One sequence label is created for each haplotype, defining a *path family* for all paths consistent with the haplotype.

The script supports standard sequence alleles (REF/ALT), as well a subset of symbolic ALT alleles:

*  Sequence alleles: SNPs and small indels / substitutions produce subgraphs with two paths: one for the reference, and one for the alternate sequence.
*  Deletions: explicit deletions or symbolic `<DEL>` alleles segment the reference sequence into three parts. The middle part can be bypassed using an edge. The `END` INFO field can be used to indicate the length of the deletion for long deletions.
* Forward breakends on the same chromosome, e.g. REF: `A`; ALT: `AC[chr1:10000[`, can be used to encode
  long deletions with short inserted sequence (swaps / substitutions).
* Insertions with symbolic `<INS>` allele: the `SEQ` INFO field must be non-empty and contain the reference sequence. The `END` INFO field may contain the end of the reference sequence that is replaced with the inserted sequence. Note that our convention is to assume that the first reference base in the record is a padding base, and that this base is not present in the value of `INFO/SEQ`.
* When variants have an ID value then it must be unique and will be used to label the variant edges in the graph. A VCF file with non-unique IDs or other problematic INFO/FORMAT fields can be cleaned up using bcftools:
```bash
bcftools annotate --set-id '.' -x 'INFO,FORMAT,^FORMAT/GT' \
        input.vcf.gz -o input-for-vcf2paragraph.py
```

All variants in one VCF passed to VCF2paragraph.py must be on the same chromosome. Our
other wrapper scripts (e.g. `multigrmpy.py`) will split up the VCF file into individual records or superloci (clusters of close-by variants) before conversion.

Paragraph doesn't support cycles in the graphs, and VCF2paragraph will raise an error if the resulting
graph contains a cycle (e.g. created by a breakend). To break cycles, vcf2paragraph may split long reference nodes and remove long reference stretches (our graph alignments start from BAM files, so realignments to long reference stretches are not necessary).

Run `vcf2paragraph.py -h` for full list of command line options.

### <a name='addVariants.py'></a>addVariants.py

addVariants takes a graph (in paragraph JSON format) and a set of paragraph variant calls on that graph, to produce a new graph with the variants included explicitly as nodes and edges. Reads can then be re-aligned to that new graph for accurate genotyping of the variants.

Variants can be in any node (ref or alt) and are relative to the sequence of that node. addVariants will split the node into a left-flank, the variant and a right-flank (either flank maybe empty). For the variant part two nodes will be created one with the original node sequence, one with the alternative sequence (either maybe empty in case of insertions / deletions). The new nodes are connected to the flanks and the whole subgraph replaces the original node in the graph. Existing sequence labels are maintained, but no labels are added to the new edges for the added variant; It has to be halplotyped in a later alignment (paragraph) step.

Run `addVariants.py -h` for full list of command line options.

### <a name='compare-alignments.py'></a>compare-alignments.py

This is a helper script to compare two paragraph alignment JSON files.

### <a name='findgrm.py'></a>findgrm.py

This is a helper script / module which is used by the other scripts to find the installation base directory and set up Python library path for a paragraph installation.

### <a name='msa2vcf.py'></a>msa2vcf.py

This script converts multiple sequence alignment Fasta files (Fasta files where sequences may contain `-` characters to denote gaps) into VCF files which can then be converted into graphs.

### <a name='paragraph2dot.py'></a>paragraph2dot.py

This script is used to convert paragraph graph JSON files to dot files for visualisation.


