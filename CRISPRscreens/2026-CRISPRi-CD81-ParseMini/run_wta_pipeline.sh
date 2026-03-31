#!/bin/bash

GENOME_DIR="/home/users/jo912684/data/genomes/GRCh38/"
KIT="WT_mini"
CHEM="v3"

# Loop over all WT sublibraries present in the directory
for SUB in $(ls BM076_WT_Sublibrary_*_R1_001.fastq.gz | sed 's/.*Sublibrary_\([0-9][0-9]\).*/\1/' | sort -u)
do
    echo "Processing Sublibrary $SUB"

    # Collect all R1 and R2 files for this sublibrary
    FQ1=$(ls BM076_WT_Sublibrary_${SUB}_*_R1_001.fastq.gz | tr '\n' ' ')
    FQ2=$(ls BM076_WT_Sublibrary_${SUB}_*_R2_001.fastq.gz | tr '\n' ' ')

    # Define output directory
    OUTDIR="output_sublibrary_${SUB}"

    # Run split-pipe
    split-pipe \
        --mode all \
        --chemistry ${CHEM} \
        --kit ${KIT} \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --genome_dir ${GENOME_DIR} \
        --output_dir ${OUTDIR} \
        --samp_slt BM076_sample_loading_table.xlsm

    echo "Finished Sublibrary $SUB"
    echo "----------------------------------------"
done
