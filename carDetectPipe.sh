#!/bin/bash

# Define directory paths
SEQUENCE_DIR="data/sequence"
RAPMAP_INDEX_DIR="data/indexes/rapmap"
SALMON_INDEX_DIR="data/indexes/salmon"

# Define executable paths
RAPMAP_EXEC="docker run --rm -v ${PWD}/data:/home/data jaeddy/rapmap:0.1.0-pre rapmap"
SAMTOOLS_EXEC="samtools"
SALMON_EXEC="tools/SalmonBeta-0.5.0_OSX-10.10/bin/salmon"

# Build indexes for the modified reference transcriptome, which includes the
# gene construct sequence as an additional transcript
HG19_CAR_FASTA="${SEQUENCE_DIR}/hg19_CAR_transcripts.fa"

# Build RapMap quasiindex for hg19 + CAR transcript
HG19_CAR_RAPMAP_INDEX="${INDEX_DIR}/hg19_CAR"
if [ ! -e "$HG19_CAR_RAPMAP_INDEX" ]; then
    echo "Building hg19 + CAR RapMap quasiindex..."
    $RAPMAP_EXEC quasiindex \
        -t $HG19_CAR_FASTA \
        -i $HG19_CAR_RAPMAP_INDEX \
        -k 19
fi

# Build RapMap quasiindex for hg19 + CAR transcript
HG19_CAR_SALMON_INDEX="${INDEX_DIR}/hg19_CAR"
if [ ! -e "$HG19_CAR_SALMON_INDEX" ]; then
    echo "Building hg19 + CAR Salmon quasiindex..."
    $SALMON_EXEC index \
        -t $HG19_CAR_FASTA \
        -i $HG19_CAR_SALMON_INDEX \
        --type quasi \
        -k 19
fi

# Run RapMap and Salmon on all input FASTQs
SAMPLE_LIST=#

# Run RapMap on all samples
while read input; do
    in_file=${input##*/}
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

# Run Salmon on all samples
while read fastq; do
    $SALMON_EXEC quant \
        -i $GENE_XCRIPT_PLUS_INDEX \
        -l U \
        -r $SIM_FASTQ_NO_MULTI \
        -o $GENE_XCRIPT_PLUS_OUT_DIR_1
done < #
