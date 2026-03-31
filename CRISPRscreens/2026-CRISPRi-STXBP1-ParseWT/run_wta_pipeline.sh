#!/bin/bash

set -euo pipefail

GENOME_DIR="/home/users/jo912684/data/genomes/GRCh38/"
KIT="WT"
CHEM="v3"
SAMPLE_LOADING_TABLE="BM084_Parse_WT_SLT_v2.1_260323.xlsm"

shopt -s nullglob

# Loop over all WT sublibraries present in the directory
for SUB in $(printf '%s\n' BM084-WT-Sublibrary-*_S*_L001_R1_001.fastq.gz | sed -E 's/.*Sublibrary-([0-9][0-9])_.*/\1/' | sort -u)
do
    echo "Processing WT Sublibrary $SUB"

    # Collect all R1 and R2 files for this sublibrary
    FQ1=$(printf '%s ' BM084-WT-Sublibrary-${SUB}_S*_L*_R1_001.fastq.gz)
    FQ2=$(printf '%s ' BM084-WT-Sublibrary-${SUB}_S*_L*_R2_001.fastq.gz)

    # Define output directory
    OUTDIR="output_WT_sublibrary_${SUB}"

    # Run split-pipe
    split-pipe \
        --mode all \
        --chemistry "${CHEM}" \
        --kit "${KIT}" \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --genome_dir "${GENOME_DIR}" \
        --output_dir "${OUTDIR}" \
        --samp_slt "${SAMPLE_LOADING_TABLE}"

    echo "Finished WT Sublibrary $SUB"
    echo "----------------------------------------"
done