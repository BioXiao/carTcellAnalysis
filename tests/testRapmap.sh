#!/bin/bash

SEQUENCE_DIR="data/sequence"
INDEX_DIR="data/indexes/rapmap"

RAPMAP_EXEC="docker run --rm -v ${PWD}/data:/home/data jaeddy/rapmap:0.1.0-pre rapmap"
SAMTOOLS_EXEC="samtools"

# Set k-mer size
K=31

### BUILD RAPMAP INDEXES

# Build quasiindex for full CAR transcript
CAR_FASTA="${SEQUENCE_DIR}/carTranscript.fasta"
CAR_INDEX="${INDEX_DIR}/carTranscript_k${K}"

if [ ! -e "$CAR_INDEX" ]; then
    echo "Building CAR transcript quasiindex..."
    $RAPMAP_EXEC quasiindex \
        -t $CAR_FASTA \
        -i $CAR_INDEX \
        -k $K
fi

# Build pseudoindex for full CAR transcript
CAR_PINDEX="${INDEX_DIR}/carTranscript_k${K}_p"

if [ ! -e "$CAR_PINDEX" ]; then
    echo "Building CAR transcript pseudoindex..."
    $RAPMAP_EXEC pseudoindex \
        -t $CAR_FASTA \
        -i $CAR_PINDEX \
        -k 19
fi

# # Build quasiindex for GRCh38 + CAR transcript
# GRCH38_CAR_FASTA="${SEQUENCE_DIR}/GRCh38_CAR_transcripts.fa"
# GRCH38_CAR_INDEX="${INDEX_DIR}/GRCh38_CAR"
#
# if [ ! -e "$GRCH38_CAR_INDEX" ]; then
#     echo "Building GRCh38 + CAR quasiindex..."
#     $RAPMAP_EXEC quasiindex \
#         -t $GRCH38_CAR_FASTA \
#         -i $GRCH38_CAR_INDEX \
#         -k 19
# fi


### TEST SALMON WITH SIMULATED READS

SIM_DATA_DIR="data/simFastqs"
SIM_FASTQ="${SIM_DATA_DIR}/simCarParts.fastq"

OUT_DIR="data/results/rapMap"
if [ ! -e "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

# Test 1

# Map simulated reads to CAR transcript
SIM_OUT="${OUT_DIR}/simTest_k${K}"

if [ ! -e "${SIM_OUT}.bam" ]; then
    echo "Quasimap ${SIM_FASTQ}"
    $RAPMAP_EXEC quasimap \
        -i $CAR_INDEX \
        -r $SIM_FASTQ \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $SIM_OUT && \
        $SAMTOOLS_EXEC index ${SIM_OUT}.bam
fi

### TEST SALMON WITH REAL LIBRARY

DATA_DIR="data/P89"
LIB7582_FASTQ="${DATA_DIR}/lib7582_C6VC6ANXX_trimmed.fastq"

# Test 2

# Map simulated reads to CAR transcript
LIB7582_OUT="${OUT_DIR}/lib7582_k${K}"

if [ ! -e "${LIB7582_OUT}.bam" ]; then
    echo "Quasimap ${LIB7582_FASTQ} to ${CAR_INDEX}"
    $RAPMAP_EXEC quasimap \
        -i $CAR_INDEX \
        -r $LIB7582_FASTQ \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $LIB7582_OUT && \
        $SAMTOOLS_EXEC index ${LIB7582_OUT}.bam
fi

# Test 3

# Pseudoalign simulated reads to CAR transcript
LIB7582_POUT="${OUT_DIR}/lib7582_k${K}_p"

if [ ! -e "${LIB7582_POUT}.bam" ]; then
    echo "Pseudomap ${LIB7582_FASTQ} to ${CAR_PINDEX}"
    $RAPMAP_EXEC pseudomap \
        -i $CAR_PINDEX \
        -r $LIB7582_FASTQ \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $LIB7582_POUT && \
        $SAMTOOLS_EXEC index ${LIB7582_POUT}.bam
fi

# Test 4

# # Map simulated reads to CAR transcript
# LIB7582_XOUT="${OUT_DIR}/lib7582_x"
#
# if [ ! -e "${LIB7582_XOUT}.bam" ]; then
#     echo "Quasimap ${LIB7582_FASTQ} to ${GRCH38_CAR_INDEX}"
#     $RAPMAP_EXEC quasimap \
#         -i $GRCH38_CAR_INDEX \
#         -r $LIB7582_FASTQ \
#         | $SAMTOOLS_EXEC view -bS - \
#         | $SAMTOOLS_EXEC sort - $LIB7582_XOUT && \
#         $SAMTOOLS_EXEC index ${LIB7582_XOUT}.bam
# fi
