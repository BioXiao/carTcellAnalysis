#!/bin/bash

SEQUENCE_DIR="data/sequence"
INDEX_DIR="data/indexes/rapmap"

RAPMAP_EXEC="docker run --rm -v ${PWD}/data:/home/data jaeddy/rapmap:0.1.0-pre rapmap"
SAMTOOLS_EXEC="samtools"

# Set k-mer size
K=19

### BUILD RAPMAP INDEXES

# Build quasiindex for full CAR transcript
CAR_FASTA="${SEQUENCE_DIR}/carTranscript.fasta"
CAR_INDEX="${INDEX_DIR}/carTranscript"

# if [ ! -e "$CAR_INDEX" ]; then
#     echo "Building CAR transcript quasiindex..."
#     $RAPMAP_EXEC quasiindex \
#         -t $CAR_FASTA \
#         -i $CAR_INDEX \
#         -k $K
# fi
#
# # Build pseudoindex for full CAR transcript
# CAR_PINDEX="${INDEX_DIR}/carTranscript_k${K}_p"
#
# if [ ! -e "$CAR_PINDEX" ]; then
#     echo "Building CAR transcript pseudoindex..."
#     $RAPMAP_EXEC pseudoindex \
#         -t $CAR_FASTA \
#         -i $CAR_PINDEX \
#         -k 19
# fi

# Build quasiindex for GRCh38 + CAR transcript
GRCH38_CAR_FASTA="${SEQUENCE_DIR}/GRCh38_CAR_transcripts.fa"
GRCH38_CAR_INDEX="${INDEX_DIR}/GRCh38_CAR"

if [ ! -e "$GRCH38_CAR_INDEX" ]; then
    echo "Building GRCh38 + CAR quasiindex..."
    $RAPMAP_EXEC quasiindex \
        -t $GRCH38_CAR_FASTA \
        -i $GRCH38_CAR_INDEX \
        -k 19
fi

# Build quasiindex for GRCh38 + CAR transcript
HG38_CAR_FASTA="${SEQUENCE_DIR}/hg38_CAR_transcripts.fa"
HG38_CAR_INDEX="${INDEX_DIR}/hg38_CAR"

if [ ! -e "$HG38_CAR_INDEX" ]; then
    echo "Building GRCh38 + CAR quasiindex..."
    $RAPMAP_EXEC quasiindex \
        -t $HG38_CAR_FASTA \
        -i $HG38_CAR_INDEX \
        -k 19
fi


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
    echo "Quasimap ${SIM_FASTQ} to ${CAR_INDEX}"
    $RAPMAP_EXEC quasimap \
        -i $CAR_INDEX \
        -r $SIM_FASTQ \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $SIM_OUT && \
        $SAMTOOLS_EXEC index ${SIM_OUT}.bam
fi

### TEST SALMON WITH REAL LIBRARY

DATA_DIR="data/P89_bulk"

# Test on a few libs

XCRIPT_LIST=$(grep ">" ${SEQUENCE_DIR}/combined_CAR_endo.fa \
              | awk -F" " '{gsub(">", ""); print $1}')

while read fastq; do
    in_file=${fastq##*/}
    lib_id=${in_file%.*}
    out_prefix="data/results/rapmap/P89_bulk_hg38/${lib_id}_hg38_CAR"

    printf "\n>> Quasimap %s to %s...\n" "${in_file}" "${HG38_CAR_INDEX}"
    # if [ ! -e "${out_prefix}.bam" ]; then
        time $RAPMAP_EXEC quasimap -t 4 \
        -i $HG38_CAR_INDEX \
        -r $fastq \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $out_prefix && \
        $SAMTOOLS_EXEC index ${out_prefix}.bam
    # fi
done < <(find ${DATA_DIR} -name "*fastq*")

while read fastq; do
    in_file=${fastq##*/}
    lib_id=${in_file%.*}
    out_prefix="data/results/rapmap/P89_bulk_CAR/${lib_id}_CAR"

    printf "\n>> Quasimap %s to %s...\n" "${in_file}" "${CAR_INDEX}"
    # if [ ! -e "${out_prefix}.bam" ]; then
        time $RAPMAP_EXEC quasimap -t 4 \
        -i $CAR_INDEX \
        -r $fastq \
        | $SAMTOOLS_EXEC view -bS - \
        | $SAMTOOLS_EXEC sort - $out_prefix && \
        $SAMTOOLS_EXEC index ${out_prefix}.bam
    # fi
done < <(find ${DATA_DIR} -name "*fastq*")
