#!/bin/bash

SEQUENCE_DIR="data/sequence"
INDEX_DIR="data/indexes/salmon"
SALMON_EXEC="/Users/jaeddy/code/github/resources/SalmonBeta-0.5.0_OSX-10.10/bin/salmon"

### BUILD SALMON INDEXES

# Treat each part of the gene construct as a separate transcript
GENE_PARTS_FASTA="${SEQUENCE_DIR}/carGeneParts.fasta"
GENE_PARTS_INDEX="${INDEX_DIR}/carGeneParts"

if [ ! -e "$GENE_PARTS_INDEX" ]; then
    echo "Building gene parts index..."
    $SALMON_EXEC index \
        -t $GENE_PARTS_FASTA \
        -i $GENE_PARTS_INDEX \
        --type quasi \
        -k 19
fi

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

# Create an index with only the full gene construct as well as individual parts
# included as transcripts
GENE_XCRIPT_PLUS_FASTA="${SEQUENCE_DIR}/carTranscriptPlus.fasta"
cat $GENE_XCRIPT_FASTA $GENE_PARTS_FASTA > $GENE_XCRIPT_PLUS_FASTA
GENE_XCRIPT_PLUS_INDEX="${INDEX_DIR}/carTranscriptPlus"

if [ ! -e "$GENE_XCRIPT_PLUS_INDEX" ]; then
    echo "Building transcript plus parts index..."
    $SALMON_EXEC index \
        -t $GENE_XCRIPT_PLUS_FASTA \
        -i $GENE_XCRIPT_PLUS_INDEX \
        --type quasi \
        -k 19
fi

# Build the index for the reference transcriptome
XCRIPTOME_FASTA="${SEQUENCE_DIR}/GRCh38_transcripts.fa"
XCRIPTOME_INDEX="${INDEX_DIR}/GRCh38"

if [ ! -e "$XCRIPTOME_INDEX" ]; then
    echo "Building reference transcriptome index..."
    $SALMON_EXEC index \
        -t $XCRIPTOME_FASTA \
        -i $XCRIPTOME_INDEX \
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


### TEST SALMON WITH SIMULATED READS

SIM_DATA_DIR="data/simFastqs"
SIM_FASTQ_ALL="${SIM_DATA_DIR}/simCarParts.fastq"
SIM_FASTQ_NO_MULTI="${SIM_DATA_DIR}/simCarParts_noMulti.fastq"
SIM_FASTQ_MULTI_ONLY="${SIM_DATA_DIR}/simCarParts_multiOnly.fastq"

OUT_DIR="data/results"
if [ ! -e "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

# Test 1

# Map simulated reads to gene construct parts
GENE_PART_OUT_DIR_1="${OUT_DIR}/genePartTest1"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $SIM_FASTQ_NO_MULTI \
    -o $GENE_PART_OUT_DIR_1

# Count total mapped reads
total1_1=$(grep -v "#" ${GENE_PART_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

GENE_PART_OUT_DIR_2="${OUT_DIR}/genePartTest2"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $SIM_FASTQ_MULTI_ONLY \
    -o $GENE_PART_OUT_DIR_2

# Count total mapped reads
total1_2=$(grep -v "#" ${GENE_PART_OUT_DIR_2}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')



# Test 2

# Map simulated reads to gene transcript
GENE_XCRIPT_OUT_DIR="${OUT_DIR}/transcriptTest"

$SALMON_EXEC quant \
    -i $GENE_XCRIPT_INDEX \
    -l U \
    -r $SIM_FASTQ_ALL \
    -o $GENE_XCRIPT_OUT_DIR

# Count total mapped reads
total2=$(grep -v "#" ${GENE_XCRIPT_OUT_DIR}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')


# Test 3

# Map simulated reads to full gene transcript as well as individual parts
GENE_XCRIPT_PLUS_OUT_DIR_1="${OUT_DIR}/transcriptPlusTest1"

$SALMON_EXEC quant \
    -i $GENE_XCRIPT_PLUS_INDEX \
    -l U \
    -r $SIM_FASTQ_NO_MULTI \
    -o $GENE_XCRIPT_PLUS_OUT_DIR_1

# Count total mapped reads
total3_1=$(grep -e "CAR-1" ${GENE_XCRIPT_PLUS_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

GENE_XCRIPT_PLUS_OUT_DIR_2="${OUT_DIR}/transcriptPlusTest2"

$SALMON_EXEC quant \
    -i $GENE_XCRIPT_PLUS_INDEX \
    -l U \
    -r $SIM_FASTQ_MULTI_ONLY \
    -o $GENE_XCRIPT_PLUS_OUT_DIR_2

# Count total mapped reads
total3_2=$(grep -e "CAR-1" ${GENE_XCRIPT_PLUS_OUT_DIR_2}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

GENE_XCRIPT_PLUS_OUT_DIR_3="${OUT_DIR}/transcriptPlusTest3"

$SALMON_EXEC quant \
    -i $GENE_XCRIPT_PLUS_INDEX \
    -l U \
    -r $SIM_FASTQ_ALL \
    -o $GENE_XCRIPT_PLUS_OUT_DIR_3

# Count total mapped reads
total3_3=$(grep -e "CAR-1" ${GENE_XCRIPT_PLUS_OUT_DIR_3}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')


# Test 4

XCRIPTOME_OUT_DIR_1="${OUT_DIR}/transcriptomeTest1"

$SALMON_EXEC quant \
    -i $XCRIPTOME_INDEX \
    -l U \
    -r $SIM_FASTQ_ALL \
    -o $XCRIPTOME_OUT_DIR_1

# Count total mapped reads
total4_1=$(grep -v "#" ${XCRIPTOME_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')


XCRIPTOME_OUT_DIR_2="${OUT_DIR}/transcriptomeTest2"

$SALMON_EXEC quant \
    -i $XCRIPTOME_INDEX \
    -l U \
    -r $SIM_FASTQ_MULTI_ONLY \
    -o $XCRIPTOME_OUT_DIR_2

# Count total mapped reads
total4_2=$(grep -v "#" ${XCRIPTOME_OUT_DIR_2}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Report test results

cat << EOM

################################################################
TEST RESULTS
################################################################

>>> Test 1:
Intra-part reads should map more successfully to separate parts than reads that
span multiple parts. Because reads are split into k-mers for mapping, some (if
not most) of the inter-part reads will still map.

EOM

printf "%s intra-part reads mapped\n" $total1_1
printf "%s inter-part reads mapped\n" $total1_2


cat << EOM

>>> Test 2:
All reads, whether intra-part or inter-part, should map to the full gene
construct sequence (treated as a transcript by Salmon).

EOM

printf "%s combined reads mapped\n" $total2

cat << EOM

>>> Test 3:
When both the full gene and individual parts are included in the transcript
index, the quantification model should be able to distinguish the following
scenarios:

1. If all reads only map to individual parts, this provides very little
evidence for expression of the full transcript relative to the various parts.
2. If all reads span multiple parts, this provides very strong evidence that
the full gene is expressed, while the levels of parts should be very low.
3. For a mixture of intra- and inter-part reads, the model should still predict
relatively high expression of the full transcript, with some low-to-moderate
expression of individual parts.

EOM

printf "%s intra-part reads mapped to full gene/transcript\n" $total3_1
printf "%s inter-part reads mapped to full gene/transcript\n" $total3_2
printf "%s combined reads mapped to full gene/transcript\n" $total3_3

cat << EOM

>>> Test 4:
Only some reads (a subset of those mapping to individual parts) should map to
the reference transcriptome. Very few (or none) of the part-spanning reads
should map.

EOM

printf "%s combined reads mapped\n" $total4_1
printf "%s inter-part reads mapped\n" $total4_2
