PARAGRAPH=~/winhome/dev/paralog/paragraph/src/python/
HG19=/mnt/c/Users/fschlesing/graphs/genome/hg19.fa
HG38=/mnt/c/Users/fschlesing/graphs/genome/hg38.fa

# vcf2paragraph
ls simple/*vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -r $HG38"
ls simple/*.json | parallel "python3 $PARAGRAPH/bin/paragraph2dot.py {} >(dot -Tpng > {.}.png)"
 
ls pg-het-ins/*vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -r $HG19"
ls haplo-complex/*.vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -r $HG38"

python3 $PARAGRAPH/bin/vcf2paragraph.py long-del/chr4-21369091-21376907.vcf >(python3 -mjson.tool --sort-keys >long-del/chr4-21369091-21376907.json) -r $HG19

ls insertions/*.vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -r {.}.ref.fa -g alleles --alt-splitting --read-len 5 --max-ref-node-length 10 --alt-paths --retrieve-reference-sequence"
ls insertions/*.vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.noas.json) -r {.}.ref.fa -g alleles --read-len 5 --max-ref-node-length 10 --alt-paths --retrieve-reference-sequence"

ls pg-complex/*.vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -R -r $HG38 -g alleles"

ls ../genotyping_test_2/chr*.vcf | parallel "python3 $PARAGRAPH/bin/vcf2paragraph.py {} >(python3 -mjson.tool --sort-keys >{.}.json) -r ../genotyping_test_2/swaps.fa -R -g alleles"

# Alignment
$PARAGRAPH/bin/paragraph -b long-del/chr4-21369091-21376907.bam -g long-del/chr4-21369091-21376907.json -o >(python3 -mjson.tool --sort-keys > long-del/chr4-21369091-21376907.paragraph.json) -r $HG19
$PARAGRAPH/bin/paragraph -r ${HG38} -g phasing/long-phasing.json -b phasing/na12878-many-haps.bam -E 1 -o >(python3 -mjson.tool --sort-keys > phasing/expected.json) --bad-align-uniq-kmer-len 64

# Variants
python3 $PARAGRAPH/bin/addVariants.py -v variants/ref.json - | python3 -mjson.tool --sort-keys > variants/ref-vars.json
