# Process raw fastqs for each sublibrary of the WTA (transcriptome) and CRISPR (sgRNAs) libraries.

## ------------------------------------------------------------
## STEP 0: Activate conda environment.
## ------------------------------------------------------------

```bash
cd /home/users/jo912684
source miniconda3/etc/profile.d/conda.sh
conda activate spipe
cd data/Seqmatic_122025_Mega_MPRA
```

## ------------------------------------------------------------
## STEP 1: Loop over all available WTA sublibraries.
## ------------------------------------------------------------
```bash
#!/bin/bash

GENOME_DIR="/home/users/jo912684/data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/genomes/GRCh38/"
KIT="WT_mega"
CHEM="v3"

# Loop over all WT sublibraries present in the directory
for SUB in $(ls BM075_WT_Sublibrary_*_R1_001.fastq.gz | sed 's/.*Sublibrary_\([0-9][0-9]\).*/\1/' | sort -u)
do
    echo "Processing Sublibrary $SUB"

    # Collect all R1 and R2 files for this sublibrary
    FQ1=$(ls BM075_WT_Sublibrary_${SUB}_*_R1_001.fastq.gz | tr '\n' ' ')
    FQ2=$(ls BM075_WT_Sublibrary_${SUB}_*_R2_001.fastq.gz | tr '\n' ' ')

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
        --samp_slt BM075_Sample_Loading.xlsm

    echo "Finished Sublibrary $SUB"
    echo "----------------------------------------"
done
```

## ------------------------------------------------------------
## STEP 2: Merge WTA sublibraries into a combined matrix.
## ------------------------------------------------------------

Note: The regular "--mode comb" command crashed since it's trying to merge count matrices of 1.7 million (filtered) or ~14 million cells (unfiltered). Therefore, I crafted scripts to merge files of the sublibrary in a memory efficient manner, which are located in the 'utils' folder.


```bash
python3 utils/merge_dge_filtered.py \
  --sublibs \
    output_sublibrary_01 \
    output_sublibrary_02 \
    output_sublibrary_03 \
    output_sublibrary_04 \
    output_sublibrary_05 \
    output_sublibrary_06 \
    output_sublibrary_07 \
    output_sublibrary_08 \
    output_sublibrary_09 \
    output_sublibrary_10 \
    output_sublibrary_11 \
    output_sublibrary_12 \
    output_sublibrary_13 \
    output_sublibrary_14 \
    output_sublibrary_15 \
    output_sublibrary_16 \
    --outdir output_combined_filtered \
  --gzip_out

  python3 utils/merge_dge_unfiltered.py \
  --sublibs \
    output_sublibrary_01 \
    output_sublibrary_02 \
    output_sublibrary_03 \
    output_sublibrary_04 \
    output_sublibrary_05 \
    output_sublibrary_06 \
    output_sublibrary_07 \
    output_sublibrary_08 \
    output_sublibrary_09 \
    output_sublibrary_10 \
    output_sublibrary_11 \
    output_sublibrary_12 \
    output_sublibrary_13 \
    output_sublibrary_14 \
    output_sublibrary_15 \
    output_sublibrary_16 \
    --outdir output_combined_filtered \
  --gzip_out
```

Now need to add the sublibrary information to the barcodes to match Parse's pipeline.

```bash
awk -F',' 'BEGIN{OFS=","}
NR==1 {
    print $0",bc_wells_s"
    next
}
{
    match($1, /([0-9]+)$/, m)
    s_num = m[1] + 0
    print $0, $2 "__s" s_num
}' output_combined_filtered/cell_metadata.csv > output_combined_filtered/cell_metadata_sub.csv

awk -F',' 'BEGIN{OFS=","}
NR==1 {
    print $0",bc_wells_s"
    next
}
{
    match($1, /([0-9]+)$/, m)
    s_num = m[1] + 0
    print $0, $2 "__s" s_num
}' output_combined_unfiltered/cell_metadata.csv > output_combined_unfiltered/cell_metadata_sub.csv
```

## ------------------------------------------------------------
## STEP 3: Check sgRNA list for potential errors
## ------------------------------------------------------------

We encountered an error due to several sgRNAs having a Hamming distance => 1. Therefore, the following script identifies such cases and separate conflicting sgRNAs into subsets.

```bash
python utils/partition_guides_mm1.py
# Make sure the script is calling the proper guides CSV file (in this case: guides_MPRA_dedup.csv)
```

This identified 2 subsets:
- 4639 sgRNAs
- 5 sgRNAs

Moving forward, I dropped the 5 conflicting sgRNAs.

### Loop over CRISPR sublibraries using the new subset of sgRNAs.

```bash
#!/bin/bash

GUIDES="guides_subset_1.csv"
KIT="WT_mega"
CHEM="v3"

# Loop over all CRISPR sublibraries present
for SUB in $(ls BM075_CRISPR_Sublibrary_*_R1_001.fastq.gz | sed 's/.*Sublibrary_\([0-9][0-9]\).*/\1/' | sort -u)
do
    echo "Processing CRISPR Sublibrary $SUB"

    # Collect all R1 and R2 files for this CRISPR sublibrary
    FQ1=$(ls BM075_CRISPR_Sublibrary_${SUB}_*_R1_001.fastq.gz | tr '\n' ' ')
    FQ2=$(ls BM075_CRISPR_Sublibrary_${SUB}_*_R2_001.fastq.gz | tr '\n' ' ')

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
```

## ------------------------------------------------------------
## STEP 4: Merge CRISPR sublibraries into a combined matrix.
## ------------------------------------------------------------

Since the CRISPR matrices are much smaller (sgRNAs by total cells), this process is not memory-intensive and can be ran using Parse's commands.

The script can call any of the WTA sublibraries via 'parent_dir' to extract the chemistry, version, and kit information.

```bash
split-pipe --mode comb \
--parent_dir output_sublibrary_01 \
--output_dir output_combined_crispr \
--sublibraries output_crispr_sublibrary_01 output_crispr_sublibrary_02 output_crispr_sublibrary_03 output_crispr_sublibrary_04 output_crispr_sublibrary_05 output_crispr_sublibrary_06 output_crispr_sublibrary_07 output_crispr_sublibrary_08 output_crispr_sublibrary_09 output_crispr_sublibrary_10 output_crispr_sublibrary_11 output_crispr_sublibrary_12 output_crispr_sublibrary_13 output_crispr_sublibrary_14 output_crispr_sublibrary_15 output_crispr_sublibrary_16

```
