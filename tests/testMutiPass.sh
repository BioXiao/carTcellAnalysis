#!/bin/bash

SEQUENCE_DIR="data/sequence"
INDEX_DIR="data/indexes/salmon"
SALMON_EXEC="tools/SalmonBeta-0.5.0_OSX-10.10/bin/salmon"

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


### TEST SALMON WITH P89-7 SAMPLES

DATA_DIR="data/P89"
LIB7536_FASTQ="${DATA_DIR}/lib7536_C6VC6ANXX_trimmed.fastq"
LIB7536_FASTQ="${DATA_DIR}/lib7536_C6VC6ANXX_trimmed.fastq"
LIB7549_FASTQ="${DATA_DIR}/lib7549_C6VC6ANXX_trimmed.fastq"
LIB7582_FASTQ="${DATA_DIR}/lib7582_C6VC6ANXX_trimmed.fastq"


OUT_DIR="data/tests"
if [ ! -e "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

# Test 1

# Map simulated reads to gene construct parts
LIB7536_OUT_DIR_1="${OUT_DIR}/lib7536Test1"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $LIB7536_FASTQ \
    -o $LIB7536_OUT_DIR_1 \
    --extraSensitive \
    -k 13

# Count total mapped reads
total1_1a=$(grep -v "#" ${LIB7536_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Count reads mapped to exogenous segments
total1_1b=$(grep -e "CD19scFv" -e "T2A" ${LIB7536_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Test 2

# Map simulated reads to gene construct parts
LIB7546_OUT_DIR_1="${OUT_DIR}/lib7546Test1"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $LIB7546_FASTQ \
    -o $LIB7546_OUT_DIR_1 \
    --extraSensitive \
    -k 13

# Count total mapped reads
total2_1a=$(grep -v "#" ${LIB7546_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Count reads mapped to exogenous segments
total2_1b=$(grep -e "CD19scFv" -e "T2A" ${LIB7546_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Test 3

# Map simulated reads to gene construct parts
LIB7549_OUT_DIR_1="${OUT_DIR}/lib7549Test1"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $LIB7549_FASTQ \
    -o $LIB7549_OUT_DIR_1 \
    --extraSensitive \
    -k 13

# Count total mapped reads
total3_1a=$(grep -v "#" ${LIB7549_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Count reads mapped to exogenous segments
total3_1b=$(grep -e "CD19scFv" -e "T2A" ${LIB7549_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Test 3

# Map simulated reads to gene construct parts
LIB7582_OUT_DIR_1="${OUT_DIR}/lib7582Test1"

$SALMON_EXEC quant \
    -i $GENE_PARTS_INDEX \
    -l U \
    -r $LIB7582_FASTQ \
    -o $LIB7582_OUT_DIR_1 \
    --extraSensitive \
    -k 13

# Count total mapped reads
total4_1a=$(grep -v "#" ${LIB7582_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')

# Count reads mapped to exogenous segments
total4_1b=$(grep -e "CD19scFv" -e "T2A" ${LIB7582_OUT_DIR_1}/quant.sf \
    | awk '{ sum+=$4} END {print sum}')


### REPORT TEST RESULTS

SUMMARY_FILE="data/results/testSummary_P89-7.txt"

cat > $SUMMARY_FILE << EOM

################################################################
TEST RESULTS
################################################################

>>> Test 1:
lib7536:

EOM

printf "%s total reads mapped\n" $total1_1a >> $SUMMARY_FILE
printf "%s reads mapped to exogenous segments\n" $total1_1b >> $SUMMARY_FILE


cat >> $SUMMARY_FILE << EOM

>>> Test 2:
lib7536:

EOM

printf "%s total reads mapped\n" $total2_1a >> $SUMMARY_FILE
printf "%s reads mapped to exogenous segments\n" $total2_1b >> $SUMMARY_FILE


cat >> $SUMMARY_FILE << EOM

>>> Test 3:
lib7549:

EOM

printf "%s total reads mapped\n" $total3_1a >> $SUMMARY_FILE
printf "%s reads mapped to exogenous segments\n" $total3_1b >> $SUMMARY_FILE


cat >> $SUMMARY_FILE << EOM

>>> Test 4:
lib7582:

EOM

printf "%s total reads mapped\n" $total4_1a >> $SUMMARY_FILE
printf "%s reads mapped to exogenous segments\n" $total4_1b >> $SUMMARY_FILE


cat $SUMMARY_FILE
