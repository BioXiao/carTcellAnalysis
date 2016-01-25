#!/bin/bash

SEQUENCE_DIR="data/sequence"
INDEX_DIR="data/indexes/rapmap"

RAPMAP_EXEC="docker run --rm -v ${PWD}/data:/home/data jaeddy/rapmap:0.1.0-pre rapmap"

### BUILD SALMON INDEXES

# Treat each part of the gene construct as a separate transcript
GENE_PARTS_FASTA="${SEQUENCE_DIR}/carGeneParts.fasta"
GENE_PARTS_INDEX="${INDEX_DIR}/carGeneParts"

if [ ! -e "$GENE_PARTS_INDEX" ]; then
    echo "Building gene parts index..."
    $RAPMAP_EXEC quasiindex \
        -t $GENE_PARTS_FASTA \
        -i $GENE_PARTS_INDEX \
        -k 19
fi

# ### TEST SALMON WITH SIMULATED READS
#
# SIM_DATA_DIR="data/simFastqs"
# SIM_FASTQ_NO_MULTI="${SIM_DATA_DIR}/simCarParts_noMulti.fastq"
#
# OUT_DIR="data/results"
# if [ ! -e "$OUT_DIR" ]; then
#     mkdir -p $OUT_DIR
# fi
#
# # Test 1
#
# # Map simulated reads to gene construct parts
# GENE_PART_OUT_DIR="${OUT_DIR}/rapMapTest"
#
# $RAPMAP_EXEC quasimap \
#     -i $GENE_PARTS_INDEX \
#     -r $SIM_FASTQ_NO_MULTI \
#     -o $GENE_PART_OUT_DIR
