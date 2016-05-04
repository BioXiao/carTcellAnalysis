#!usr/bin/env bash

# This script is used to prepare reference sequence and annotation data used
# for

# Specify whether to force generation of formatted reference data, even if
# files already exist
FORCE=0
while getopts "f" ARG; do
	case "$ARG" in
	    f ) FORCE=1;;
	esac
done
shift $(($OPTIND - 1))


# Specify raw sequence text file input; file should already exist in a folder
# labeled 'sequence' under the 'data' directory
SEQUENCE_DIR="data/sequence"
RAW_CAR_SEQ="${SEQUENCE_DIR}/car_transcript_raw.txt"

# Specify human genome FASTA and gene model GTF input; file should exist in a
# folder labeled 'annotation' under the 'data' directory
HG38_GENOME_FASTA="${SEQUENCE_DIR}/hg38/hg38_genome.fa"

ANNOTATION_DIR="data/annotation"
HG38_XCRIPT_GTF="${ANNOTATION_DIR}/hg38/hg38_genes.gtf"

### FORMAT REFERENCE SEQUENCE

# Convert raw sequence text file to FASTA with records for each part of the
# gene construct
#
# RAW_GENE_SEQ: unformatted sequences, copied from Word doc into text file
# NAME: name of gene/transcript/construct (used for labelling)
# CAR_XCRIPT_FASTA: output FASTA file for merged transcript sequence
# MERGE: if 'True', merge individual segments into a single FASTA record
# TRANSCRIPT: if 'True', output is transcript (don't add artificial 'intron'
#   buffer between segments)

NAME="CAR"
CAR_XCRIPT_FASTA="${SEQUENCE_DIR}/car_transcript.fa"
MERGE="True"
TRANSCRIPT="True"

if [ ! -e "$GENE_XCRIPT_FASTA" ] || [ ${FORCE} == 1 ]; then
    python scripts/format_fasta.py \
        $RAW_CAR_SEQ \
        $NAME \
        $CAR_XCRIPT_FASTA \
        $MERGE \
        $TRANSCRIPT
fi


### CONVERT HUMAN GENE MODEL GTF TO FASTA

# HG38_XCRIPT_FASTA: output transcriptome FASTA
# HG38_GENOME_FASTA: reference human genome (hg38) from iGenomes
# HG38_XCRIPT_GTF: reference gene models for hg38 from iGenomes

HG38_XCRIPT_FASTA="${ANNOTATION_DIR}/hg38/hg38_transcripts.fa"

if [ ! -e "$HG38_XCRIPT_FASTA" ] || [ ${FORCE} == 1 ]; then
    gffread \
        -w $HG38_XCRIPT_FASTA
        -g $HG38_GENOME_FASTA
        $HG38_XCRIPT_GTF


### ADD SEQUENCE TO TRANSCRIPTOME FASTA

# Concatenate FASTA record for gene construct with whole-transcriptome
# reference FASTA
#
# CAR_XCRIPT_FASTA: CAR transcript FASTA file
# HG38_XCRIPT_FASTA: hg38 transcriptome FASTA file
# CUSTOM_XCRIPTOME_FASTA: concatenated transcriptome FASTA

CUSTOM_XCRIPTOME_FASTA="${SEQUENCE_DIR}/hg38/hg38_car_transcripts.fa"

if [ ! -e "$CUSTOM_XCRIPTOME_FASTA" ] || [ ${FORCE} == 1 ]; then
    cat $CAR_XCRIPT_FASTA \
        $HG38_XCRIPT_FASTA \
        > $CUSTOM_XCRIPTOME_FASTA
fi


### CREATE REFERENCE GENE MODEL GTF

# Convert raw CAR sequence into FASTA again, but this time create separate
# records for each segment.

# RAW_CAR_SEQ: unformatted sequences, copied from Word doc into text file
# NAME: name of gene/transcript/construct (used for labelling)
# CAR_PARTS_FASTA: FASTA file with a separate record for each segment of the
#   CAR construct
# MERGE: when 'False', creates a separate labeled record for each segment of
#   the CAR construct
# TRANSCRIPT: default is 'False', but has no effect when `MERGE` is 'False'

CAR_PARTS_FASTA="${SEQUENCE_DIR}/car_parts.fa"
MERGE="False"
TRANSCRIPT="False"

if [ ! -e "$CAR_PARTS_FASTA" ]; || [ ${FORCE} == 1 ] then
    python scripts/format_fasta.py \
        $RAW_CAR_SEQ \
        $NAME \
        $CAR_PARTS_FASTA \
        $MERGE \
        $TRANSCRIPT
fi

# From CAR parts FASTA, create gene model GTF; transcript ID is used as
# chromosome name; each segment is labeled as a unique transcript (note:
# segments are trimmed by 1bp from each direction such that the start/end
# is easier to see with track visualization like IGV)

# CAR_PARTS_FASTA: FASTA file with a separate record for each segment of the
#   CAR construct
# CAR_PARTS_GTF: name of output gene model GTF file
# SPACER: number of bp to trim from both ends of each part to create a visual
#   gap when viewing in browser

CAR_PARTS_GTF="${ANNOTATION_DIR}/car_parts.gtf"
SPACER=1

if [ ! -e "$CAR_PARTS_GTF" ] || [ ${FORCE} == 1 ]; then
    python scripts/parts_fasta_to_gtf.py \
        $CAR_PARTS_FASTA \
        $CAR_PARTS_GTF \
        $SPACER
fi

# Create pseudo-GTF records for human transcripts that overlap with CAR
# segments; transcript ID is used as chromosome name; common gene name is used
# as gene ID; corresponding CAR segment is used as transcript ID

# CAR_PARTS_OVERLAP_GTF: gene model records for transcripts that overlap with
#   CAR construct sequence
# CAR_PARTS_OVERLAP_GTF: name of output gene model GTF file
# SPACER: number of bp to trim from both ends of each part to create a visual
#   gap when viewing in browser

CAR_PARTS_OVERLAP_GTF="${ANNOTATION_DIR}/car_parts_overlap.gtf"
SPACER=0

if [ ! -e "$CAR_PARTS_OVERLAP_GTF" ] || [ ${FORCE} == 1 ]; then
    python scripts/parts_fasta_to_gtf.py \
        $CAR_PARTS_OVERLAP_FASTA \
        $CAR_PARTS_OVERLAP_GTF \
        $SPACER
fi

# Finally, create combined GTF file

CAR_PLUS_OVERLAP_GTF="${ANNOTATION_DIR}/car_plus_overlap.gtf"

if [ ! -e "$CAR_PLUS_REF_GTF" || [ ${FORCE} == 1 ]]; then
    cat \
        $CAR_PARTS_GTF \
        $CAR_PARTS_OVERLAP_GTF \
        > $CAR_PLUS_OVERLAP_GTF
fi
