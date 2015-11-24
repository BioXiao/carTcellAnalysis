

salmon index \
    -t sequence/GRCh38_CAR_transcripts.fa \
    -i GRCh38_CAR --type quasi -k 19

salmon quant \
    -i indexes/salmon/GRCh38_CAR \
    -l U \
    -r /Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-5Processed_150916/TrimmedFastqs/lib8576_C6REVANXX_trimmed.fastq \
    -o transcripts_quant6
