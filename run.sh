#!/bin/bash

# INPUT_DIR="data/P89/P89-5"
# OUTPUT_BASE="data/results/P89-5_test"
#
# python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE


INPUT_DIR="/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-5Processed_150916/TrimmedFastqs/"
OUTPUT_BASE="data/results/P89-5_150910"

python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE


INPUT_DIR="/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-7Processed_150916/TrimmedFastqs/"
OUTPUT_BASE="data/results/P89-7_150910"

python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE


INPUT_DIR="/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-9Processed_150916/TrimmedFastqs/"
OUTPUT_BASE="data/results/P89-9_150910"

python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE


INPUT_DIR="/Volumes/genomics/Illumina/150729_D00565_0092_AC6VC6ANXX/Project_P89-7Processed_150916/TrimmedFastqs/"
OUTPUT_BASE="data/results/P89-7_150729"

python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE


INPUT_DIR="/Volumes/genomics/Illumina/150729_D00565_0092_AC6VC6ANXX/Project_P89-9Processed_150916/TrimmedFastqs/"
OUTPUT_BASE="data/results/P89-9_150729"

python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE
