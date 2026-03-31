#!/bin/bash

GUIDES="guides.csv"
KIT="WT"
CHEM="v3"

# Loop over all CRISPR sublibraries present
for SUB in $(ls BM073_CRISPR_Sublibrary_*_R1_001.fastq.gz | sed 's/.*Sublibrary_\([0-9][0-9]\).*/\1/' | sort -u)
do
    echo "Processing CRISPR Sublibrary $SUB"

    # Collect all R1 and R2 files for this CRISPR sublibrary
    FQ1=$(ls BM073_CRISPR_Sublibrary_${SUB}_*_R1_001.fastq.gz | tr '\n' ' ')
    FQ2=$(ls BM073_CRISPR_Sublibrary_${SUB}_*_R2_001.fastq.gz | tr '\n' ' ')

    # Parent directory from WT pipeline
    PARENT_DIR="output_sublibrary_${SUB}"

    # Output directory for CRISPR run
    OUTDIR="output_crispr_sublibrary_${SUB}"

    # Run split-pipe CRISPR mode
    split-pipe \
        --mode all \
        --crispr \
        --crsp_guides ${GUIDES} \
        --parent_dir ${PARENT_DIR} \
        --output_dir ${OUTDIR} \
        --fq1 ${FQ1} \
        --fq2 ${FQ2} \
        --chemistry ${CHEM} \
        --kit ${KIT} \
        --chem_score_skip

    echo "Finished CRISPR Sublibrary $SUB"
    echo "----------------------------------------"
done