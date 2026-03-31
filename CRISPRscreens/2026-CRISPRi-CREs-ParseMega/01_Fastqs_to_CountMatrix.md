# Process raw fastqs for each sublibrary of the WTA (transcriptome) and CRISPR (sgRNAs) libraries.

## ------------------------------------------------------------
## STEP 0: Activate conda environment.
## ------------------------------------------------------------

```bash
cd /home/users/jo912684
source miniconda3/etc/profile.d/conda.sh
conda activate spipe
cd data/Seqmatic_022026_Parse_Mega_CREs
```

## ------------------------------------------------------------
## STEP 1: Loop over all available WTA sublibraries.
## ------------------------------------------------------------
```bash
#!/bin/bash
set -euo pipefail
shopt -s nullglob

GENOME_DIR="/home/users/jo912684/data/genomes/GRCh38/"
KIT="WT_mega"
CHEM="v3"
SLT="02a-BM081-Parse_Biosciences_Evercode_WT_Mega_INTEGRA_Sample_Loading_Table_v2.1-260204_-_Copy.xlsm"

# find sublibrary numbers like 01, 02, ... from filenames
subs=()
for f in BM081_WT_Sublibrary??_*_R1_001.fastq.gz; do
  [[ $f =~ Sublibrary([0-9]{2})_ ]] && subs+=( "${BASH_REMATCH[1]}" )
done
# unique sort
mapfile -t subs < <(printf "%s\n" "${subs[@]}" | sort -u)

for SUB in "${subs[@]}"; do
  echo "Processing Sublibrary $SUB"

  fq1=( BM081_WT_Sublibrary${SUB}_*_R1_001.fastq.gz )
  fq2=( BM081_WT_Sublibrary${SUB}_*_R2_001.fastq.gz )

  if ((${#fq1[@]} == 0 || ${#fq2[@]} == 0)); then
    echo "WARNING: no FASTQs found for Sublibrary $SUB (fq1=${#fq1[@]}, fq2=${#fq2[@]}). Skipping."
    continue
  fi

  OUTDIR="output_sublibrary_${SUB}"

  split-pipe \
    --mode all \
    --chemistry "${CHEM}" \
    --kit "${KIT}" \
    --fq1 "${fq1[@]}" \
    --fq2 "${fq2[@]}" \
    --genome_dir "${GENOME_DIR}" \
    --output_dir "${OUTDIR}" \
    --samp_slt "${SLT}"

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
    --outdir output_combined_filtered

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
    --outdir output_combined_filtered
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
## STEP 3: Check sgRNA list to make sure there are no issues
## ------------------------------------------------------------

In a previous run, we encountered an error due to several sgRNAs having a Hamming distance => 1. Therefore, the following script identifies such cases and separates conflicting sgRNAs into subsets.

```bash
python utils/partition_guides_mm1.py
# Make sure the script is calling the proper guides CSV file (in this case: guides_CRE.csv)
```

This identified 2 subsets:
- 6665 sgRNAs
- 34 sgRNAs

I dropped out the 34 sgRNAs moving forward to avoid issues.

## ------------------------------------------------------------
## STEP 4: Loop over CRISPR libraries
## ------------------------------------------------------------

```bash
#!/bin/bash
set -euo pipefail
shopt -s nullglob

GUIDES="guides_CRE_subset_1.csv"
KIT="WT_mega"
CHEM="v3"

CRISPR_PREFIX="BM081_CRISPR_Sublibrary"
WT_PARENT_PREFIX="output_sublibrary_"

# Discover sublibraries like 01, 02, ... from existing CRISPR FASTQs
subs=()
for f in ${CRISPR_PREFIX}??_*_R1_001.fastq.gz; do
  [[ $f =~ Sublibrary([0-9]{2})_ ]] && subs+=( "${BASH_REMATCH[1]}" )
done
mapfile -t subs < <(printf "%s\n" "${subs[@]}" | sort -u)

if ((${#subs[@]} == 0)); then
  echo "ERROR: No CRISPR R1 FASTQs found matching ${CRISPR_PREFIX}??_*_R1_001.fastq.gz"
  exit 1
fi

for SUB in "${subs[@]}"; do
  echo "Processing CRISPR Sublibrary ${SUB}"

  fq1=( ${CRISPR_PREFIX}${SUB}_*_R1_001.fastq.gz )
  fq2=( ${CRISPR_PREFIX}${SUB}_*_R2_001.fastq.gz )

  if ((${#fq1[@]} == 0 || ${#fq2[@]} == 0)); then
    echo "WARNING: Missing FASTQs for SUB=${SUB} (fq1=${#fq1[@]}, fq2=${#fq2[@]}). Skipping."
    continue
  fi

  PARENT_DIR="${WT_PARENT_PREFIX}${SUB}"
  if [[ ! -d "${PARENT_DIR}" ]]; then
    echo "WARNING: WT parent dir not found: ${PARENT_DIR}. Skipping SUB=${SUB}."
    continue
  fi

  if [[ ! -s "${GUIDES}" ]]; then
    echo "ERROR: Guides file not found or empty: ${GUIDES}"
    exit 1
  fi

  OUTDIR="output_crispr_sublibrary_${SUB}"

  split-pipe \
    --mode all \
    --crispr \
    --crsp_guides "${GUIDES}" \
    --parent_dir "${PARENT_DIR}" \
    --output_dir "${OUTDIR}" \
    --fq1 "${fq1[@]}" \
    --fq2 "${fq2[@]}" \
    --chemistry "${CHEM}" \
    --kit "${KIT}" \
    --chem_score_skip

  echo "Finished CRISPR Sublibrary ${SUB}"
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