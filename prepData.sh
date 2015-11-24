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
MERGE="True"
TRANSCRIPT="True"
GENE_XCRIPT_FASTA="${SEQUENCE_DIR}/carTranscript.fasta"

python scripts/formatFasta.py $RAW_GENE_SEQ $GENE_NAME $GENE_XCRIPT_FASTA \
    $MERGE $TRANSCRIPT


### ADD SEQUENCE TO TRANSCRIPTOME FASTA

# Concatenate FASTA record for gene construct with whole-transcriptome
# reference FASTA
CUSTOM_XCRIPTOME_FASTA="${SEQUENCE_DIR}/GRCh38_CAR_transcripts.fa"
cat $GENE_XCRIPT_FASTA "${SEQUENCE_DIR}/GRCh38_transcripts.fa" > \
    $CUSTOM_XCRIPTOME_FASTA
