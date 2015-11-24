#!/bin/bash

# This script is used to generate simulated reads from the gene construct
# sequence for testing mapping & quantification procedures.

# Specify raw sequence text file input
RAW_GENE_SEQ=$1

### FORMAT REFERENCE SEQUENCE

# Create folder to store formatted reference sequence data
SEQUENCE_DIR="data/sequence"
if [ ! -e "$SEQUENCE_DIR" ]; then
    mkdir -p $SEQUENCE_DIR
fi

# Convert raw sequence text file to FASTA with records for each part of the
# gene construct
GENE_NAME="CAR"
GENE_PARTS_FASTA="${SEQUENCE_DIR}/carGeneParts.fasta"
python scripts/formatFasta.py $RAW_GENE_SEQ $GENE_NAME $GENE_PARTS_FASTA

### SIMULATE FASTQ READS

# Create folder for outputs
SIM_DATA_DIR="data/simFastqs"
if [ ! -e "$SIM_DATA_DIR" ]; then
    mkdir -p $SIM_DATA_DIR
fi

# Number of simulated reads to generate
NUM_READS=20

# Generate reads that map only to individual parts in the gene construct
MODE=1
SIM_FASTQ_NO_MULTI="${SIM_DATA_DIR}/simCarParts_noMulti.fastq"
python scripts/simulateReads.py $GENE_PARTS_FASTA $SIM_FASTQ_NO_MULTI \
    $NUM_READS $MODE

# Generate only reads that span multiple adjacent parts in the gene construct
MODE=2
SIM_FASTQ_MULTI_ONLY="${SIM_DATA_DIR}/simCarParts_multiOnly.fastq"
python scripts/simulateReads.py $GENE_PARTS_FASTA $SIM_FASTQ_MULTI_ONLY \
    $NUM_READS $MODE

# Generate a combination of intra-part and part-spanning reads
MODE=3
SIM_FASTQ_ALL="${SIM_DATA_DIR}/simCarParts.fastq"
python scripts/simulateReads.py $GENE_PARTS_FASTA $SIM_FASTQ_ALL \
    $NUM_READS $MODE
