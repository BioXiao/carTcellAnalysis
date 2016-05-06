#!/bin/bash

SAMPLE_LIST=$1

# Define directory paths
SEQUENCE_DIR="data/sequence"
RAPMAP_INDEX_DIR="data/indexes/rapmap"
SALMON_INDEX_DIR="data/indexes/salmon"

# Define executable paths
RAPMAP_IMG="docker run --rm -v ${PWD}/data:/home/data jaeddy/rapmap:0.2.1"
SAMTOOLS_EXEC="samtools"
SALMON_EXEC="tools/SalmonBeta-0.5.0_OSX-10.10/bin/salmon"

# Build indexes for the modified reference transcriptome, which includes the
# gene construct sequence as an additional transcript
HG38_CAR_FASTA="${SEQUENCE_DIR}/hg38/hg38_car_transcripts.fa"

# Build RapMap quasiindex for hg38 + CAR transcript
HG38_CAR_RAPMAP_INDEX="${RAPMAP_INDEX_DIR}/hg38_car"
if [ ! -e "$HG38_CAR_RAPMAP_INDEX" ]; then
    echo "Building hg38 + CAR RapMap quasiindex..."
    $RAPMAP_IMG rapmap quasiindex \
        -t $HG38_CAR_FASTA \
        -i $HG38_CAR_RAPMAP_INDEX \
        -k 19
else
    echo "hg38 + CAR RapMap quasiindex already exists; skipping."
fi

# Build RapMap quasiindex for hg38 + CAR transcript
HG38_CAR_SALMON_INDEX="${SALMON_INDEX_DIR}/hg38_car"
if [ ! -e "$HG38_CAR_SALMON_INDEX" ]; then
    echo "Building hg38 + CAR Salmon quasiindex..."
    $SALMON_EXEC index \
        -t $HG38_CAR_FASTA \
        -i $HG38_CAR_SALMON_INDEX \
        --type quasi \
        -k 19
else
    echo "hg38 + CAR Salmon quasiindex already exists; skipping."
fi

RAPMAP_OUT_DIR="data/rapmap"
SALMON_OUT_DIR="data/salmon"

XCRIPT_LIST=$(grep ">" ${SEQUENCE_DIR}/car_plus_overlap.fa \
              | awk -F" " '{gsub(">", ""); print $1}')

# Run RapMap and Salmon on all samples
while read line; do
    start_time=$(date +%s)
    sleep 1

    file=$(echo $line | awk '{print $1}')
    lib_id=$(echo $file | egrep -o -e "lib[0-9]+")
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
    tmp_dir="data/input"
    if [ ! -e "$tmp_dir" ]; then
        mkdir "$tmp_dir"
    fi
    tmp_prefix="${lib_id}"
    out_prefix="${RAPMAP_OUT_DIR}/${out_subdir}/${lib_id}"
    printf "\n>> Quasimap %s to %s...\n" "${lib_id}" "${HG38_CAR_RAPMAP_INDEX}"
    if [ ! -e "${out_prefix}.bam" ]; then
        printf "\n...getting input FASTQ...\n"
        if [ ! -e "${tmp_dir}/$(basename $file)" ]; then
            time rsync -h --progress $file "${tmp_dir}/$(basename $file)"
        fi
        read -r -d '' rapmap_cmd << EOF
printf "\n...running RapMap...\n" && \
    time rapmap quasimap -t 4 \
    -i $HG38_CAR_RAPMAP_INDEX \
    -r "${tmp_dir}/$(basename $file)" \
    | samtools view -bS -@ 4 - > ${tmp_prefix}.bam && \
    printf "\n...sorting BAM output...\n" && \
    time samtools sort -@ 4 -m 2G -o ${tmp_prefix}_sorted.bam ${tmp_prefix}.bam && \
    samtools index ${tmp_prefix}_sorted.bam && \
    printf "\n...filtering sorted BAM...\n" && \
    time samtools view -b -@ 4 ${tmp_prefix}_sorted.bam $XCRIPT_LIST \
    > ${out_prefix}.bam && \
    samtools index ${out_prefix}.bam
EOF
        $RAPMAP_IMG /bin/bash -c "$(echo $rapmap_cmd)"
    else
        echo "Output already exists."
    fi
    rm -rf $tmp_dir

    # Run Salmon
    out_prefix="${SALMON_OUT_DIR}/${out_subdir}/${lib_id}"
    printf "\n>> Quant %s with %s...\n" "${lib_id}" "${HG38_CAR_SALMON_INDEX}"
    if [ ! -e "${out_prefix}/quant.sf" ]; then
        printf "\n...running Salmon...\n"
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
