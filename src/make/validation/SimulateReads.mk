#
# uses EAGLE simulator to generate reads from graph path sequences
#
# Author: Roman Petrovski <rpetrovski@illumina.com>

# The shell pipefail option is required for accurate workflow control
SHELL = /bin/bash -o pipefail

ifeq (,$(EAGLE_SHARE))
$(error "EAGLE_SHARE must point to EAGLE installation share subfolder")
endif

ifeq (,$(EXEC_PATH))
$(error "EXEC_PATH is not defined")
endif

ifeq (,$(GRAPH_FA))
$(error "GRAPH_FA must point to fasta containing path sequences for simulation")
endif

ifeq (,$(REFERENCE_FA))
$(error "REFERENCE_FA must point to whole genome fasta")
endif


$(REFERENCE_FA).fai:
	echo $@ is required to exist. && exit 2

EAGLE/Makefile:
	configureEAGLE.pl   --run-info=$(EAGLE_SHARE)/RunInfo/RunInfo_PairedReadsBarcode8x32Tiles.xml    --reference-genome=$(GRAPH_FA)   --variant-list=$(EAGLE_SHARE)/Variants/None.vcf   --coverage-depth=30   --genome-mutator-options="--organism-ploidy=2"

EAGLE/eagle.bam: EAGLE/Makefile
	$(MAKE) -C EAGLE bam

simulated.bam: EAGLE/eagle.bam $(REFERENCE_FA).fai
	$(EXEC_PATH)/swapReference.sh $^ |samtools view -Sb - >$@.tmp && mv $@.tmp $@

simulated.bam.bai: simulated.bam
	samtools index $<
