#!/bin/bash

GENE_REF=$1
GENE_MODEL=$2
FASTQ_IN=$3
OUT_DIR=$4

REF_BASE=$(basename $GENE_REF)
FILE_BASE=$(echo $FASTQ_IN | egrep -e "lib[0-9]+(_.*XX)*" -o)
OUT_BASE="${FILE_BASE}_to_${REF_BASE}"

# Align with bowtie2
TH_OUT_DIR="${OUT_BASE}_tophat_out"
BAM_OUT="${OUT_BASE}.bam"
JUNC_BED_OUT="${OUT_BASE}.junctions.bed"

tophat -G $GENE_MODEL -o $TH_OUT_DIR $GENE_REF $FASTQ_IN
mv ${TH_OUT_DIR}/accepted_hits.bam $BAM_OUT
mv ${TH_OUT_DIR}/junctions.bed $JUNC_BED_OUT
rm -rf $TH_OUT_DIR

# Count reads
COUNTS_OUT="${OUT_BASE}.counts.txt"
featureCounts -f -a $GENE_MODEL -o $COUNTS_OUT $BAM_OUT

RESULTS_OUT="${OUT_BASE}.results.txt"
python scripts/summarizeResults.py $GENE_MODEL $JUNC_BED_OUT $COUNTS_OUT $RESULTS_OUT
