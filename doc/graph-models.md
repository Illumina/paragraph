# Graph model summary

Each graph model is based on a reference sequence and alternate paths 
associated with alternative alleles or some characteristics of the region 
we would like to evaluate read support for.

Our graphs may be created specifically for a particular event type 
(deletion, insertion, swap) or from a VCF file. In this document, we describe
the main concepts behind our graph genotyping methods.

<!-- vscode-markdown-toc -->
* [Graph Models](#GraphModels)
	* [Example 1 -- Simple Swap](#Example1--SimpleSwap)
	* [Example 2 -- Short swap](#Example2--Shortswap)
	* [Example 3 -- Short deletion](#Example3--Shortdeletion)
	* [Example 4 -- Long deletion](#Example4--Longdeletion)
	* [Example 5 -- Long Swap](#Example5--LongSwap)
* [Breakpoint Genotyping](#BreakpointGenotyping)
* [Whole Variant Genotyping](#WholeVariantGenotyping)
* [Example](#Example)
	* [Breakpoints](#Breakpoints)
	* [Whole Variant](#WholeVariant)

<!-- vscode-markdown-toc-config
	numbering=false
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->


## <a name='GraphModels'></a>Graph Models

Graphs may be created from VCF or for particular events specified in 
a JSON file. For details on the VCF to graph conversion see 
[graph-tools.md#vcf2paragraph.py](graph-tools.md#vcf2paragraph.py)
and [README.md#Usingvcf2paragraph.pytorunparaGRAPH](../README.md#Usingvcf2paragraph.pytorunparaGRAPH).

Each graph we create contains the following information:

*  **Sequence labels**: Sequence labels are used to tag particular paths through our graph. 
   Each sequence label corresponds to a set of *paths* through the graph.
*  **Nodes**: Each node represents a piece of sequence that is at least one nucleotide long. 
   Nodes may be associated with one or more *sequence labels*. The sequence for each node
   may be specified explicitly, or can be retrieved from a reference FASTA file. Nodes
   may thus be associated with linear reference locations.
*  **Edges**: Edges connect nodes in a directed fashion (`from` node &rarr; `to` node). 
   Our graphs are required to be directed and acyclic, and edges may be associated with one or 
   more sequence *sequence labels*.
*  **Paths**: A path is specified by a list of node names, where we required forward edges between each 
   pair of subsequent nodes to exist (i.e. we must be able to walk along our graph along this path).
   Each path is associated with exactly one sequence label (but one sequence label may correspond
   to multiple paths). 

Graphs may be created for different event types. The following examples 
illustrate a few such events.

### <a name='Example1--SimpleSwap'></a>Example 1 -- Simple Swap

Manually created JSON file to illustrate how to name nodes. See [swap-example-1.json](swap-example-1.json).

Creating these files from VCF is better since it creates edge labels,
paths and disambiguation sequences too.

![Manual Swap](../share/test-data/paragraph/simple/swap-example-1.png)

### <a name='Example2--Shortswap'></a>Example 2 -- Short swap

```bash
vcf2paragraph.py swap-example-2.vcf swap-example-2.json -r hg19.fa -p 100
paragraph2dot.py swap-example-2.json swap-example-2.dot
dot -Tsvg swap-example-2.dot > swap-example-2.svg
```

Short swap with REF and ALT branch:

![Short Swap Example](../share/test-data/paragraph/simple/swap-example-2.png)

### <a name='Example3--Shortdeletion'></a>Example 3 -- Short deletion

```bash
vcf2paragraph.py del-example-3.vcf del-example-3.json -r hg19.fa -p 100
paragraph2dot.py del-example-3.json del-example-3.dot
dot -Tsvg del-example-3.dot > del-example-3.svg
```

This is a short deletion with two straightforward paths: REF and ALT.

![Short Deletion Example](../share/test-data/paragraph/simple/del-example-3.png)

### <a name='Example4--Longdeletion'></a>Example 4 -- Long deletion

```bash
vcf2paragraph.py del-example-4.vcf del-example-3.json -r hg19.fa -p 100
paragraph2dot.py del-example-4.json del-example-3.dot
dot -Tsvg del-example-4.dot > del-example-3.svg
```

This is a longer deletion. We get two branches, one for each breakpoint.

![Long Deletion Example](../share/test-data/paragraph/simple/del-example-4.png)

### <a name='Example5--LongSwap'></a>Example 5 -- Long Swap

```bash
vcf2paragraph.py swap-example-5.vcf swap-example-5.json -r hg19.fa -p 100
paragraph2dot.py swap-example-5.json swap-example-5.dot
dot -Tsvg swap-example-5.dot > swap-example-5.svg
```

This is a longer swap. We get two branches, one for each breakpoint.

![Long Swap Example](../share/test-data/paragraph/simple/swap-example-5.png)

## <a name='BreakpointGenotyping'></a>Breakpoint Genotyping

Breakpoints in a linear reference produce branchings in our graphs. We consider 
two types of breakpoint: forward breakpoints where one node has multiple edges
going out, and reverse breakpoints where one node has multiple edges incoming 
(this excludes source and sink nodes).

![Breakpoint types](breakpoint-types.png)

Each incoming / outgoing edge corresponds to one *allele*. Breakpoint genotyping
uses read counts on these edges to determine the most likely one or two alleles
supported by the sequencing experiment.

Each breakpoint is genotyped independently. After genotyping all breakpoints
for a particular event we compare the genotypes at each breakpoint to check 
the calls are consistent. To do this, we use our sequence labels on the edges
to name each allele -- each sequence label is assumed to support one particular
allele/haplotype for each event.

Note that one node may correspond to two breakpoints if it has more than one 
input edge and more than one output edge. Genotyping of these two breakpoints
will be performed independently.

Consider a breakpoint in one sample. We have:

&nbsp;&nbsp;&nbsp;&nbsp;**[G]<sub>a/b</sub>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A possible underlying genotype with haplotypes named *a*, *b*.

&nbsp;&nbsp;&nbsp;&nbsp;**R**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;set of reads mapped to edges connected with this breakpoint and with this direction

&nbsp;&nbsp;&nbsp;&nbsp;**R<sub>a,b</sub>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;reads in **R** mapped to the edge of haplotype *a* or *b*.

&nbsp;&nbsp;&nbsp;&nbsp;**R<sub>≠a,b</sub>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;reads in **R** mapped to edges that belong to any other haplotypes rather than *a* or *b*.

&nbsp;&nbsp;&nbsp;&nbsp;**P(G<sub>a/b</sub>)**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Prior of genotype G<sub>a/b</sub>.
  
The likelihood of observing such reads with the underlying genotype of **G<sub>a/b</sub>**

&nbsp;&nbsp;&nbsp;&nbsp;**P(R | G<sub>a/b</sub>) =  &nbsp;P(R<sub>a,b</sub> | G<sub>a/b</sub>)&nbsp; * &nbsp;P(R<sub>≠a,b</sub> | G<sub>a/b</sub>)**

is the product of the likelihood of reads mapped to haplotype *a*,*b*, and 
the likelihood of reads not mapped haplotype *a* or *b*.

The likelihood of observing genotype **G<sub>a/b</sub>** given a set of reads **R** is then

&nbsp;&nbsp;&nbsp;&nbsp;**P(G<sub>a/b</sub> | R) &prop; &nbsp;P(G<sub>a/b</sub>)&nbsp; * &nbsp;P(R | G<sub>a/b</sub>)&nbsp;**

where &nbsp;P(G<sub>a/b</sub>)&nbsp; is the prior genotype probability for **G<sub>a/b</sub>** (e.g. a
population frequency, or uniform prior over all possible genotypes).
 
Assuming read length ***l***, average depth ***d***, and a pre-defined minimum overlap of read bases ***m***, we estimate Poisson distribution mean **λ** for **R** as

&nbsp;&nbsp;&nbsp;&nbsp;**λ = round( *d* * (*l* - *m*) / *l* )**

Let **ε<sub>a/b</sub>** be the genotype error rate for **[G]<sub>a/b</sub>**. We get the Poisson distribution mean for **R<sub>≠a,b</sub>**:

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>≠a,b</sub> = round( λ * ε<sub>a/b</sub> )**

With the following notation:

&nbsp;&nbsp;&nbsp;&nbsp;**PMF<sub>≠a,b</sub>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Poisson probability mass function with parameter **λ<sub>≠a,b</sub>**

&nbsp;&nbsp;&nbsp;&nbsp;**N<sub>≠a,b</sub>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number of reads mapped to neither haplotype *a* or *b*

The genotype likelihood part of **R<sub>≠a,b</sub>** (reads mapped to neither *a* or *b*) is

&nbsp;&nbsp;&nbsp;&nbsp;**P(R<sub>≠a,b</sub> | G<sub>a/b</sub>) = PMF<sub>≠a,b</sub>( N<sub>≠a,b</sub> )**

If *a*, *b* are the same haplotype:

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>a,b</sub> = round( λ * (1 - ε<sub>a/b</sub>) )**

With similar notation above, the genotype likelihood part of **R<sub>a,b</sub>** is

&nbsp;&nbsp;&nbsp;&nbsp;**P(R<sub>a,b</sub> | G<sub>a/b</sub>) = PMF<sub>a,b</sub>( N<sub>a,b</sub> )**

If *a*,*b* are different haplotypes:

With expected read fraction in heterozygotes of haplotype *a*, *b* as *μ<sub>a</sub>*, *μ<sub>b</sub>*

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>a</sub> = λ * μ<sub>a</sub>**

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>b</sub> = λ * μ<sub>b</sub>**

The resulting genotype likelihood part of **R<sub>a,b</sub>** (reads mapped to *a* or *b*) is

&nbsp;&nbsp;&nbsp;&nbsp;**P(R<sub>a,b</sub> | G<sub>a/b</sub>) = PMF<sub>a</sub>( N<sub>a</sub> ) * PMF<sub>b</sub>( N<sub>b</sub> ) * PMF<sub>≠a,b</sub>( N<sub>≠a,b</sub> )**.

## <a name='WholeVariantGenotyping'></a>Whole Variant Genotyping
We use the above breakpoint genotyping model to determine the most likely 
genotype for each breakpoint.

To report a final genotype for the variant, we follow these conventions: 

* If all breakpoints have the same most likely genotype, report this genotype as the overall variant genotype. The resulting filter value is *PASS*
* If some breakpoints have homozygous reference genotype, and the other breakpoints have 
another consistent non-reference genotype, report this non-reference genotype and label
the filter as *CONFLICTS*.
* If none of the above, the genotype is reported as *MISSING* (**./.**) and the filter is set to *CONFLICTS*.

## <a name='Example'></a>Example

A graph model for a sequence swap consists of four nodes representing left flank, deleted sequence, inserted sequence, and right flank. It can be represented schematically like this:

+![breakpoint-genotyper.png](sequence-swap-graph.png)

This graph has two breakpoints: breakpoint #1 defined by two "from" edges and breakpoint #2 defined by two "to" edges (see the diagram above).

### <a name='Breakpoints'></a>Breakpoints

Here we calculate **G<sub>REF/INS</sub>**, the likelihood of *REF/INS* genotype of breakpoint #1.
We have haplotype *a*, *b* representing edge *LF_REF* and *LF_INS*, respectively.

We set **N<sub>a</sub>** to the number of reads overlapping left flank / reference edge (*LF_REF*), and **N<sub>b</sub>** to the number of reads mapped to left flank / insertion edge (*LF_INS*). Assume:

&nbsp;&nbsp;&nbsp;&nbsp;**N<sub>a</sub> = 19**

&nbsp;&nbsp;&nbsp;&nbsp;**N<sub>b</sub> = 1**

Since there are no other edges, the number of reads mapped to neither edges is zero, which means

&nbsp;&nbsp;&nbsp;&nbsp;**N<sub>≠a,b</sub> = 0**

Assume the sample has depth ***d*** = 20, read length ***l*** = 100. Using the default parameter of minimum base overlap

&nbsp;&nbsp;&nbsp;&nbsp;***m* = 4**

&nbsp;&nbsp;&nbsp;&nbsp;**λ = round( d * (*l* - *m*) / *l* ) = 19**

Using the default ε as 0.05 for all genotypes:

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>≠a,b</sub> = round(λ * ε) = 1**

&nbsp;&nbsp;&nbsp;&nbsp;**P(R<sub>≠a,b</sub> | G<sub>a/b</sub>) = PMF<sub>≠a,b</sub>(N<sub>≠a,b</sub>) = PMF<sub>≠a,b</sub>(0) = 0.368**

Using the default haplotype fraction **μ** = 0.5 for both **μ<sub>a</sub>** and **μ<sub>b</sub>**

&nbsp;&nbsp;&nbsp;&nbsp;**λ<sub>a</sub> = λ<sub>a</sub> = round(λ * (1 - ε) * μ) = 10**

&nbsp;&nbsp;&nbsp;&nbsp;**P(R<sub>a,b</sub> | G<sub>a/b</sub>) = PMF<sub>a</sub>(N<sub>a</sub>) * PMF<sub>b</sub>(N<sub>b</sub>)**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**= 3.73e<sup>-3</sup> * 4.54e<sup>-4</sup> = 1.69e<sup>-6<sup>**

Using the default fixed prior for all genotypes as 1 devided by the total number of genotypes. Here we have 3 genotypes (REF/REF, REF/INS and INS/INS), so

&nbsp;&nbsp;&nbsp;&nbsp;**P(G<sub>a/b</sub>) = 1/3 ≈ 0.33**

The resulting likelihood of genotype REF/INS is:

&nbsp;&nbsp;&nbsp;&nbsp;**P(R | G<sub>a/b</sub>) = 0.33 * 0.368 * 1.69e<sup>-6</sup> = 2.08e<sup>-7</sup>**

Similarly, we can get the likelihood of genotype REF/REF as **1.08e<sup>-2</sup>**, and INS/INS as **4.11e<sup>-26</sup>**.

Therefore, the maximum likelihood genotype of this breakpoint is REF/REF.

### <a name='WholeVariant'></a>Whole Variant

As we have calculated above, genotype of BP1 is REF/REF. 
If BP2 is also genotyped as REF/REF, then the variant genotype will be REF/REF. 
The filter flag will be set as *PASS*.

If BP2 is not REF/REF and BP1 is the reference genotype, then the variant genotype will be the genotype of BP2. However, the filter flag will be set as *CONFLICTS* to indicate that not all breakpoints
support this genotype. This may happen when for example when we type an event where the 
start of the of the event is correct but the end is not.
