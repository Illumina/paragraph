# Paragraph: a suite of graph-based genotyping tools

<!-- vscode-markdown-toc -->
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Run Paragraph from VCF](#RunParagraphFromVCF)
    * [Example](#Example)
    * [Input requirements](#InputRequirements)
* [Run Paragraph on complex variants](#RunParagraphOnComplexVariants)
* [Further Information](#FurtherInformation)
	* [Documentation](#Documentation)
	* [External links](#ExternalLinks)
* [License](#License)

<!-- vscode-markdown-toc-config
	numbering=false
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

## <a name='Introduction'></a>Introduction

Accurate genotyping of known variants is a critical for analysis of whole-genome sequencing data.

Paragraph aims to facilitate these tasks by providing:
- an accurate genotyper for Structural Variations in short-read data
- a suite of graph-based tools to align and genotype complex events.

Please reference Paragraph using:

- Chen, et al (2019) [Paragraph: A graph-based structural variant genotyper for short-read sequence data](https://www.biorxiv.org/content/10.1101/635011v1). *bioRxiv*. doi: https://doi.org/10.1101/635011

Genotyping calls in this paper can be found at [paper-data/download-instructions.txt](paper-data/download-instructions.txt)

## <a name='Installation'></a>Installation

Please check [doc/Installation.md](doc/Installation.md) for system requirements and installation instructions.

## <a name='RunParagraphFromVCF'></a>Run Paragraph from VCF
### <a name='Example'></a>Example
After installation, run `multigrmpy.py` script from the build/bin directory on an example dataset as follows:

```bash
python3 bin/multigrmpy.py -i share/test-data/round-trip-genotyping/candidates.vcf \
                          -m share/test-data/round-trip-genotyping/samples.txt \
                          -r share/test-data/round-trip-genotyping/dummy.fa \
                          -o test \
```

This runs a simple genotyping example for two test samples.
*  **candidates.vcf**: this specifies candidate SV events in a vcf format.
*  **samples.txt**: Manifest that specifies some test BAM files. Tab or comma delimited.
*  **dummy.fa** a short dummy reference which only contains `chr1`

The output folder `test` then contains gzipped json for final genotypes:

```bash
$ tree test
```
```
test
├── grmpy.log            #  main workflow log file
├── genotypes.vcf.gz     #  Output VCF with individual genotypes
├── genotypes.json.gz    #  More detailed output than genotypes.vcf.gz
├── variants.vcf.gz      #  The input VCF with unique ID from Paragraph
└── variants.json.gz     #  The converted graphs from input VCF (no genotypes)
```

If successful, the last 3 lines of genotypes.vcf.gz will the same as in [expected file](share/test-data/round-trip-genotyping/expected-vcf-record.txt).

## <a name='InputRequirements'></a>Input requirements
### VCF format
paraGRAPH will independently genotype each entry of the input VCF. You can use either indel-style representation (full REF and ALT allele sequence in 4th and 5th columns) or symbolic alleles, as long as they meet the format requirement of VCF 4.0+.

It typically takes up to a few seconds to genotype a single event in one sample (single-threaded).
It took us 30 minutes to genotype ~20,000 SVs using 20 CPU cores (with I/O).

Currently we support 4 symbolic alleles:
- `<DEL>` for deletion
    - Must have END key in INFO field.
- `<INS>` for insertion
    - Must have a key in INFO field for insertion sequence (without padding base). The default key is SEQ.
    - For blockwise swap, we strongly recommend using indel-style representation, other than symbolic alleles.
- `<DUP>` for duplication
    - Must have END key in INFO field. paraGRAPH assumes the sequence between POS and END being duplicated for one more time in the alternative allele.
- `<INV>` for inversion
    - Must have END key in INFO field. paraGRAPH assumes the sequence between POS and END being reverse-complemented in the alternative allele.

### Sample Manifest
Must be tab-deliemited.

Required columns:
- ID: Each sample must have a unique ID. The output VCF will include genotypes for all samples in the manifest
- path: Path to the BAM/CRAM file.
- depth: Average depth across the genome. Can be calculated with bin/idxdepth or samtools.
- read length: Average read length (bp) across the genome.

Optional columns:

- depth sd: Specify standard deviation for genome depth. Used for the normal test of breakpoint read depth. Default is sqrt(5*depth).

- depth variance: Square of depth sd.

- sex: Affects chrX and chrY genotyping. 
Allow "male" or "M", "female" or "F", and "unknown" (quotes shouldn't be included in the manifest). 
If not specified, the sample will be treated as unknown.

## <a name='RunParagraphOnComplexVariants'></a>Run Paragraph on complex variants
For more complicated events (e.g. genotype a deletion together with its nearby SNP), you can provide a custimized JSON to paraGRAPH:

Please follow the pattern in [example JSON](share/test-data/paragraph/pg-het-ins/pg-het-ins.json) and make sure all required keys are provided. Here is a visualization of this [sample graph](share/test-data/paragraph/pg-het-ins/pg-het-ins.png).

To obtain graph alignments for this graph (including all reads), run:
```bash
bin/paragraph -b <input BAM> \
              -r <reference fasta> \
              -g <input graph JSON> \
              -o <output JSON path> \
              -E 1
```

To obtain the algnment summary, genotypes of each breakpoint, and the whole graph, run:
```bash
bin/grmpy -m <input manifest> \
          -r <reference fasta> \
          -i <input graph JSON> \
          -o <output JSON path> \
          -E 1
```

If you have multiple events listed in the input JSON, `multigrmpy.py` can help you to run multiple `grmpy` jobs together.

## <a name='FurtherInformation'></a>Further Information

### <a name='Documentation'></a>Documentation

*    More **information about all tools we provide in this package** can be found in 
    [doc/graph-tools.md](doc/graph-tools.md).

*   In [doc/graph-models.md](doc/graph-models.md) we describe the graph and genotyping 
    models we implement.

*    Some developer documentation about our code analysis and testing process can be found in 
    [doc/linting-and-testing.md](doc/linting-and-testing.md).

*    Procedures for read level alignment validation 
    [doc/validation-with-simulated-reads.md](doc/validation-with-simulated-reads.md).

*    How we count reads for variants and paths
    [doc/graph-counting.md](doc/graph-counting.md).

*    Documentation of genotyping model parameters
    [doc/genotyping-parameters.md](doc/genotyping-parameters.md).

*   [Doc/graphs-ashg-2017.pdf](doc/graphs-ashg-2017.pdf) contains the poster about this method we showed at 
    [ASHG 2017](http://www.ashg.org/2017meeting/)

### <a name='ExternalLinks'></a>External links

*   The [Illumina/Polaris](https://github.com/Illumina/Polaris) repository gives the
    short-read sequencing data we used to test our method in population.

## <a name='License'></a>License

The [LICENSE](LICENSE) file contains information about libraries and other tools we use, 
and license information for these.
