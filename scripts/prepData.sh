#!/bin/bash

# This script is used to generate simulated reads from the gene construct
# sequence for testing mapping & quantification procedures.

# Specify raw sequence text file input
RAW_GENE_SEQ=$1
SALMON_EXEC="tools/SalmonBeta-0.5.0_OSX-10.10/bin/salmon"

### FORMAT REFERENCE SEQUENCE

# Create folder to store formatted reference sequence data
SEQUENCE_DIR="data/sequence"
if [ ! -e "$SEQUENCE_DIR" ]; then
    mkdir -p $SEQUENCE_DIR
fi

# Convert raw sequence text file to FASTA with records for each part of the
# gene construct
GENE_NAME="CAR"
MERGE="True"
TRANSCRIPT="True"
GENE_XCRIPT_FASTA="${SEQUENCE_DIR}/carTranscript.fasta"

if [ ! -e "$GENE_XCRIPT_FASTA" ]; then
    python scripts/formatFasta.py $RAW_GENE_SEQ $GENE_NAME $GENE_XCRIPT_FASTA \
        $MERGE $TRANSCRIPT
fi


### ADD SEQUENCE TO TRANSCRIPTOME FASTA

# Concatenate FASTA record for gene construct with whole-transcriptome
# reference FASTA
CUSTOM_XCRIPTOME_FASTA="${SEQUENCE_DIR}/GRCh38_CAR_transcripts.fa"

if [ ! -e "$CUSTOM_XCRIPTOME_FASTA" ]; then
    cat $GENE_XCRIPT_FASTA "${SEQUENCE_DIR}/GRCh38_transcripts.fa" > \
        $CUSTOM_XCRIPTOME_FASTA
fi


### BUILD SALMON INDEXES

# Create an index with only one transcript - the sequence of the gene construct
GENE_XCRIPT_FASTA="${SEQUENCE_DIR}/carTranscript.fasta"
GENE_XCRIPT_INDEX="${INDEX_DIR}/carTranscript"

if [ ! -e "$GENE_XCRIPT_INDEX" ]; then
    echo "Building transcript index..."
    $SALMON_EXEC index \
        -t $GENE_XCRIPT_FASTA \
        -i $GENE_XCRIPT_INDEX \
        --type quasi \
        -k 19
fi

# Build index for the modified reference transcriptome, which includes the
# gene construct sequence as an additional transcript
GENE_XCRIPTOME_FASTA="${SEQUENCE_DIR}/GRCh38_CAR_transcripts.fa"
GENE_XCRIPTOME_INDEX="${INDEX_DIR}/GRCh38_CAR"

if [ ! -e "$GENE_XCRIPTOME_INDEX" ]; then
    echo "Building modified reference transcriptome index..."
    $SALMON_EXEC index \
        -t $GENE_XCRIPTOME_FASTA \
        -i $GENE_XCRIPTOME_INDEX \
        --type quasi \
        -k 19
fi
