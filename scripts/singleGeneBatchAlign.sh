#!/bin/bash

GENE_REF=$1
GENE_MODEL=$2
IN_DIR=$3
OUT_DIR=$4

while read file; do
    sbatch <(cat scripts/singleGeneAlign.sh | \
        awk -v param=$GENE_REF '{gsub("\$1", param)}1' |
        awk -v param=$GENE_MODEL '{gsub("\$2", param)}1' |
        awk -v param=$file '{gsub("\$3", param)}1' |
        awk -v param=$OUT_DIR '{gsub("\$4", param)}1')

done < <(find $IN_DIR -name "*.fastq")
