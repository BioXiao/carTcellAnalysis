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

echo "[1/8]"
if [ ! -e "$CAR_XCRIPT_FASTA" ] || [ ${FORCE} == 1 ]; then
    echo "Formatting raw sequence to produce merged FASTA..."
    python scripts/format_fasta.py \
        $RAW_CAR_SEQ \
        $NAME \
        $CAR_XCRIPT_FASTA \
        $MERGE \
        $TRANSCRIPT
else
    echo "Merged FASTA already exists; skipping."
fi


### CONVERT HUMAN GENE MODEL GTF TO FASTA

# HG38_XCRIPT_FASTA: output transcriptome FASTA
# HG38_GENOME_FASTA: reference human genome (hg38) from iGenomes
# HG38_XCRIPT_GTF: reference gene models for hg38 from iGenomes

HG38_XCRIPT_FASTA="${SEQUENCE_DIR}/hg38/hg38_transcripts.fa"

echo "[2/8]"
if [ ! -e "$HG38_XCRIPT_FASTA" ] || [ ${FORCE} == 1 ]; then
    echo "Converting reference gene model GTF to transcriptome FASTA..."
    gffread \
        -w $HG38_XCRIPT_FASTA \
        -g $HG38_GENOME_FASTA \
        $HG38_XCRIPT_GTF
else
    echo "Reference transcriptome FASTA already exists; skipping."
fi


### ADD SEQUENCE TO TRANSCRIPTOME FASTA

# Concatenate FASTA record for gene construct with whole-transcriptome
# reference FASTA
#
# CAR_XCRIPT_FASTA: CAR transcript FASTA file
# HG38_XCRIPT_FASTA: hg38 transcriptome FASTA file
# CUSTOM_XCRIPTOME_FASTA: concatenated transcriptome FASTA

CUSTOM_XCRIPTOME_FASTA="${SEQUENCE_DIR}/hg38/hg38_car_transcripts.fa"

echo "[3/8]"
if [ ! -e "$CUSTOM_XCRIPTOME_FASTA" ] || [ ${FORCE} == 1 ]; then
    echo "Adding construct transcript sequence to reference transcriptome..."
    cat $CAR_XCRIPT_FASTA \
        $HG38_XCRIPT_FASTA \
        > $CUSTOM_XCRIPTOME_FASTA
else
    echo "Combined reference+construct transcriptome already exists; skipping."
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

echo "[4/8]"
if [ ! -e "$CAR_PARTS_FASTA" ] || [ ${FORCE} == 1 ]; then
    echo "Formatting raw sequence to produce FASTA records for all segments..."
    python scripts/format_fasta.py \
        $RAW_CAR_SEQ \
        $NAME \
        $CAR_PARTS_FASTA \
        $MERGE \
        $TRANSCRIPT
else
    echo "Construct parts/segments FASTA already exists; skipping."
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

echo "[5/8]"
if [ ! -e "$CAR_PARTS_GTF" ] || [ ${FORCE} == 1 ]; then
    echo "Converting construct parts FASTA to pseudo-GTF..."
    python scripts/parts_fasta_to_gtf.py \
        $CAR_PARTS_FASTA \
        $CAR_PARTS_GTF \
        $SPACER
else
    echo "Construct parts pseudo-GTF already exists; skipping."
fi

# Extract transcripts from human reference corresponding to genes that
# overlap with CAR construct (i.e., map to CAR segments).

# HG38_XCRIPT_FASTA: hg38 transcriptome FASTA file
# NAME: name of gene/transcript/construct (used for labelling)
# GENE_MAP_YAML: YAML file with CAR segments mapped to common gene names
# CAR_PARTS_OVERLAP_FASTA: FASTA records for transcripts that overlap with
#   CAR construct sequence

GENE_MAP_YAML="${ANNOTATION_DIR}/gene_map.yaml"
CAR_PARTS_OVERLAP_FASTA="${SEQUENCE_DIR}/car_parts_overlap.fa"

echo "[6/8]"
if [ ! -e "$CAR_PARTS_OVERLAP_FASTA" ] || [ ${FORCE} == 1 ]; then
    echo "Extracting transcript FASTA records for genes overlapping with construct..."
    python scripts/get_overlap_transcripts.py \
        $HG38_XCRIPT_FASTA \
        $NAME \
        $GENE_MAP_YAML \
        $CAR_PARTS_OVERLAP_FASTA
else
    echo "Overlapping transcript FASTA already exists; skipping."
fi

# Create pseudo-GTF records for human transcripts that overlap with CAR
# segments; transcript ID is used as chromosome name; common gene name is used
# as gene ID; corresponding CAR segment is used as transcript ID

# CAR_PARTS_OVERLAP_GTF: gene model records for transcripts that overlap with
#   CAR construct sequence
# SPACER: number of bp to trim from both ends of each part to create a visual
#   gap when viewing in browser

CAR_PARTS_OVERLAP_GTF="${ANNOTATION_DIR}/car_parts_overlap.gtf"
SPACER=0

echo "[7/8]"
if [ ! -e "$CAR_PARTS_OVERLAP_GTF" ] || [ ${FORCE} == 1 ]; then
    echo "Converting construct overlap FASTA to pseudo-GTF..."
    python scripts/parts_fasta_to_gtf.py \
        $CAR_PARTS_OVERLAP_FASTA \
        $CAR_PARTS_OVERLAP_GTF \
        $SPACER
else
    echo "Construct overlap pseudo-GTF already exists; skipping."
fi

# Finally, create combined GTF file

CAR_PLUS_OVERLAP_GTF="${ANNOTATION_DIR}/car_plus_overlap.gtf"

echo "[8/8]"
if [ ! -e "$CAR_PLUS_OVERLAP_GTF" ] || [ ${FORCE} == 1 ]; then
    echo "Combining pseudo-GTF records into a single file..."
    cat \
        $CAR_PARTS_GTF \
        $CAR_PARTS_OVERLAP_GTF \
        > $CAR_PLUS_OVERLAP_GTF
else
    echo "Combined pseudo-GTF file already exists; skipping."
fi
