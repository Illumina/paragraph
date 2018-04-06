# Alignment validation with simulated reads

The following describes steps required to assess accuracy of paragraph alignment using simulated data reads.

## Generating reference sequences

In order to simulate data reads, an input reference is required. `bin/graph-to-fasta` will produce a separate contig for each path in the supplied graphs:

`/path/to/paragraph/bin/graph-to-fasta -r /path/to/full/genome.fa -g /path/to/graphs/*.json >graphs.fa`

## Running EAGLE simulation to produce reads

When building paragraph `configure --with-eagle=/path/to/EAGLE` can be used to enable paragraph/EAGLE integration. If this has not been done `libexec/simulate-reads.sh --with-eagle /path/to/EAGLE` can still be used to perform simulation. For example:

`/path/to/paragraph/libexec/simulate-reads.sh --graph-genome graphs.fa --reference-genome=/path/to/whole/genome.fa --with-eagle /path/to/EAGLE`

By default `simulate-reads.sh` will attempt to use all cores available on the node. Use `--jobs` to override.

Once the `simulate-reads.sh` completes, the current folder will contain `simulated.bam` file with reads generated in such way that paragraph will load them when processing corresponding graph.json files.

## Running validation

It is important to run validation on a single thread. Otherwise each threads will produce a validation report which might not be a desired outcome.

The following example runs validation for a combination of kmer and graph aligners:

`/path/to/paragraph/bin/paragraph -r /path/to/whole/genome.fa -g /path/to/graphs/*.json -b simulated.bam -o aligned.json --exact no --graph-se yes --kmer no --validate --threads 1`

Once paragraph completes, the log contains the validation summary table:

```
[2018-02-28 18:03:22.190] [paragraph] [info] [Done with alignment step 376 total aligned (exact: 0 / kmers: 0 / sw: 376) ; 2 were filtered]
[2018-02-28 18:03:22.190] [paragraph] [info] [VALIDATION]       MAPQ    EmpMAPQ Wrong   Total
[2018-02-28 18:03:22.190] [paragraph] [info] [VALIDATION]       unalgnd 0       0       0
[2018-02-28 18:03:22.190] [paragraph] [info] [VALIDATION]       repeat  0       0       0
[2018-02-28 18:03:22.190] [paragraph] [info] [VALIDATION]       60      60      0       376
```

The table columns have the following meaning:
* MAPQ     MAPQ assigned by validation
  * unalgnd unaligned reads
  * repeat  reads that had more than one equivalent best candiate alignment
  * 60      Uniquely aligned reads
* EmpMAPQ  MAPQ computed using the following formula: -10 * log10(Wrong/Total)
* Wrong    Number of alignments that don't support the path from which they have been simulated
* Total    Total number of aligned reads

# References

EAGLE on github [https://github.com/sequencing/EAGLE](https://github.com/sequencing/EAGLE)
