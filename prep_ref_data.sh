#!usr/bin/env bash

# This script is used to prepare reference sequence and annotation data used
# for

# Specify raw sequence text file input; file should already exist in a folder
# labeled 'sequence' under the 'data' directory
SEQUENCE_DIR="data/sequence"
RAW_CAR_SEQ="${SEQUENCE_DIR}/car_transcript_raw.txt"

# Specify human genome FASTA and gene model GTF input; file should exist in a
# folder labeled 'annotation' under the 'data' directory
HG38_GENOME_FASTA="${SEQUENCE_DIR}/hg38_genome.fa"

ANNOTATION_DIR="data/annotation"
HG38_XCRIPT_GTF="${ANNOTATION_DIR}/hg38_genes.gtf"

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

if [ ! -e "$GENE_XCRIPT_FASTA" ]; then
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

HG38_XCRIPT_FASTA="${ANNOTATION_DIR}/hg38_transcripts.fa"

if [ ! -e "$HG38_XCRIPT_FASTA" ]; then
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

CUSTOM_XCRIPTOME_FASTA="${SEQUENCE_DIR}/hg38_car_transcripts.fa"

if [ ! -e "$CUSTOM_XCRIPTOME_FASTA" ]; then
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

if [ ! -e "$CAR_PARTS_FASTA" ]; then
    python scripts/format_fasta.py \
        $RAW_CAR_SEQ \
        $NAME \
        $CAR_PARTS_FASTA \
        $MERGE \
        $TRANSCRIPT
fi

# From CAR parts FASTA, create gene model GTF; transcript ID is used as
# chromosome name; each segment is labeled as a unique transcript (note:
# segments are trimmed by 2bp from each direction such that the start/end
# is easier to see with track visualization like IGV)

# CAR_PARTS_FASTA: FASTA file with a separate record for each segment of the
#   CAR construct
# NAME: name of gene/transcript/construct (used for labelling)
# CAR_PARTS_GTF: name of output gene model GTF file

CAR_PARTS_GTF="${ANNOTATION_DIR}/car_parts.gtf"

if [ ! -e "$CAR_PARTS_GTF" ]; then
    python scripts/gene_fasta_to_gtf.py \
        $CAR_PARTS_FASTA \
        $NAME \
        $CAR_PARTS_GTF
fi

# Create pseudo-GTF records for human transcripts that overlap with CAR
# segments; transcript ID is used as chromosome name; common gene name is used
# as gene ID; corresponding CAR segment is used as transcript ID

# CAR_PARTS_OVERLAP_GTF: gene model records for transcripts that overlap with
#   CAR construct sequence

CAR_PARTS_OVERLAP_GTF="${ANNOTATION_DIR}/car_parts_overlap.gtf"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_000758' 'custom' 'transcript' '1' '784' '.' '+' '.' \
    'gene_id "CSF2"; transcript_id "GMCSFRss_r1";' \
    > $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_001243077' 'custom' 'transcript' '1' '4588' '.' '+' '.' \
    'gene_id "CD28"; transcript_id "CD28tm_r1";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_001243078' 'custom' 'transcript' '1' '4522' '.' '+' '.' \
    'gene_id "CD28"; transcript_id "CD28tm_r2";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_006139' 'custom' 'transcript' '1' '4878' '.' '+' '.' \
    'gene_id "CD28"; transcript_id "CD28tm_r3";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_001561' 'custom' 'transcript' '1' '5984' '.' '+' '.' \
    'gene_id "TNFRSF9"; transcript_id "IgG4hinge_r1";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_000734' 'custom' 'transcript' '1' '1676' '.' '+' '.' \
    'gene_id "CD247"; transcript_id "CD3Zeta_r1";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_198053' 'custom' 'transcript' '1' '1679' '.' '+' '.' \
    'gene_id "CD247"; transcript_id "CD3Zeta_r2";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_201283' 'custom' 'transcript' '1' '1570' '.' '+' '.' \
    'gene_id "EGFR"; transcript_id "EGFRt_r1";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_201282' 'custom' 'transcript' '1' '2236' '.' '+' '.' \
    'gene_id "EGFR"; transcript_id "EGFRt_r2";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_201284' 'custom' 'transcript' '1' '2861' '.' '+' '.' \
    'gene_id "EGFR"; transcript_id "EGFRt_r3";' \
    >> $CAR_PARTS_OVERLAP_GTF
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    'NM_005228' 'custom' 'transcript' '1' '5592' '.' '+' '.' \
    'gene_id "EGFR"; transcript_id "EGFRt_r4";' \
    >> $CAR_PARTS_OVERLAP_GTF

# Finally, create combined GTF file

CAR_PLUS_REF_GTF="${ANNOTATION_DIR}/car_plus_ref.gtf"

if [ ! -e "$CAR_PLUS_REF_GTF"]; then
    cat \
        $CAR_PARTS_GTF \
        $CAR_PARTS_OVERLAP_GTF \
        > $CAR_PLUS_REF_GTF
fi
