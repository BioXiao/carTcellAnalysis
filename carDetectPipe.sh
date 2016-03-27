#!/bin/bash

SAMPLE_LIST=$1

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
HG38_CAR_FASTA="${SEQUENCE_DIR}/hg38_CAR_transcripts.fa"

# Build RapMap quasiindex for hg38 + CAR transcript
HG38_CAR_RAPMAP_INDEX="${RAPMAP_INDEX_DIR}/hg38_CAR"
if [ ! -e "$HG38_CAR_RAPMAP_INDEX" ]; then
    echo "Building hg38 + CAR RapMap quasiindex..."
    $RAPMAP_EXEC quasiindex \
        -t $HG38_CAR_FASTA \
        -i $HG38_CAR_RAPMAP_INDEX \
        -k 19
fi

# Build RapMap quasiindex for hg38 + CAR transcript
HG38_CAR_SALMON_INDEX="${SALMON_INDEX_DIR}/hg38_CAR"
if [ ! -e "$HG38_CAR_SALMON_INDEX" ]; then
    echo "Building hg38 + CAR Salmon quasiindex..."
    $SALMON_EXEC index \
        -t $HG38_CAR_FASTA \
        -i $HG38_CAR_SALMON_INDEX \
        --type quasi \
        -k 19
fi

RAPMAP_OUT_DIR="data/results/rapmap"
SALMON_OUT_DIR="data/results/salmon"

XCRIPT_LIST=$(grep ">" ${SEQUENCE_DIR}/combined_CAR_endo.fa \
              | awk -F" " '{gsub(">", ""); print $1}')

# Run RapMap and Salmon on all samples
while read line; do
    start_time=$(date +%s)
    sleep 1

    file=$(echo $line | awk '{print $1}')
    lib_id=$(echo $file | egrep -o -e "lib[0-9]+.*XX")
    fc=$(echo $line | awk '{print $2}')
    proj=$(echo $line | awk '{print $3}')

    out_subdir="${fc}/${proj}"
    if [ ! -e "${RAPMAP_OUT_DIR}/${out_subdir}" ]; then
        mkdir -p "${RAPMAP_OUT_DIR}/${out_subdir}"
    fi

    if [ ! -e "${SALMON_OUT_DIR}/${out_subdir}" ]; then
        mkdir -p "${SALMON_OUT_DIR}/${out_subdir}"
    fi

    # Run RapMap and format output
    tmp_dir="data/tmp"
    if [ ! -e "$tmp_dir" ]; then
        mkdir "$tmp_dir"
    fi
    tmp_prefix="${tmp_dir}/${lib_id}"
    out_prefix="${RAPMAP_OUT_DIR}/${out_subdir}/${lib_id}"
    printf "\n>> Quasimap %s to %s...\n" "${lib_id}" "${HG38_CAR_RAPMAP_INDEX}"
    if [ ! -e "${out_prefix}.bam" ]; then
        if [ ! -e "${tmp_dir}/$(basename $file)" ]; then
            time rsync -h --progress $file "${tmp_dir}/$(basename $file)"
        fi
        time $RAPMAP_EXEC quasimap -t 4 \
            -i $HG38_CAR_RAPMAP_INDEX \
            -r "${tmp_dir}/$(basename $file)" > ${tmp_prefix}.sam

        echo "Formatting output..."
        time ($SAMTOOLS_EXEC view -bS -@ 4 ${tmp_prefix}.sam \
            > ${tmp_prefix}.bam && \
            $SAMTOOLS_EXEC sort -@ 4 ${tmp_prefix}.bam ${tmp_prefix}_sorted && \
            $SAMTOOLS_EXEC index -@ 4 ${tmp_prefix}_sorted.bam && \
            $SAMTOOLS_EXEC view -b -@ 4 ${tmp_prefix}_sorted.bam $XCRIPT_LIST \
            > ${out_prefix}.bam && \
            $SAMTOOLS_EXEC index -@ 4 ${out_prefix}.bam)
    else
        echo "Output already exists."
    fi
    rm -rf $tmp_dir

    # Run Salmon
    out_prefix="${SALMON_OUT_DIR}/${out_subdir}/${lib_id}"
    printf "\n>> Quant %s with %s...\n" "${lib_id}" "${HG38_CAR_SALMON_INDEX}"
    if [ ! -e "${out_prefix}/quant.sf" ]; then
        time $SALMON_EXEC quant \
            -i $HG38_CAR_SALMON_INDEX \
            -l U \
            -r $file \
            -o $out_prefix
    else
        echo "Output already exists."
    fi

    end_time=$(date +%s)
    time_diff=$(echo "$end_time - $start_time" | bc)
    printf "\n>> jobs for ${lib_id} completed in %s seconds\n" $time_diff

done < $SAMPLE_LIST
