# Process raw fastqs for each sublibrary of the WTA (transcriptome) and CRISPR (sgRNAs) libraries.

## ------------------------------------------------------------
## STEP 0: Activate conda environment.
## ------------------------------------------------------------

```bash
cd /home/users/jo912684
source miniconda3/etc/profile.d/conda.sh
conda activate spipe
cd data/Seqmatic_032026_Parse_Mega_druggable
```

## ------------------------------------------------------------
## STEP 1: Loop over all available WTA sublibraries via SLURM
## ------------------------------------------------------------

Trying different options to speed up the analysis, currently runtime in the HPC is ~1 week.

 Save the following script as BM082R_WT_splitpipe_array.sh
```bash
#!/bin/bash
#SBATCH --job-name=BM082R_WT_splitpipe
#SBATCH --output=/home/users/jo912684/data/Seqmatic_032026_Parse_Mega_druggable/slurm_logs/BM082R_WT_splitpipe_%A_%a.out
#SBATCH --error=/home/users/jo912684/data/Seqmatic_032026_Parse_Mega_druggable/slurm_logs/BM082R_WT_splitpipe_%A_%a.err
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=5-00:00:00
#SBATCH --array=1-16%4

set -euo pipefail
shopt -s nullglob
umask 002

# Activate conda environment
source /home/users/jo912684/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# -----------------------------
# User-configurable settings
# -----------------------------
WORKDIR="/home/users/jo912684/data/Seqmatic_032026_Parse_Mega_druggable"
GENOME_DIR="/home/users/jo912684/data/genomes/GRCh38/"
KIT="WT_mega"
CHEM="v3"
SLT="${WORKDIR}/03a-BM082-R-Parse_WT_Mega_SLT_v2.1-260226_Copy.xlsm"

# This script is for WT only
PREFIX="BM082R-WT-Sublibrary"

mkdir -p "${WORKDIR}"
mkdir -p "${WORKDIR}/slurm_logs"

cd "${WORKDIR}"

# Use Slurm CPU allocation
NTHREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Map array task 1..16 -> sublibrary 01..16
SUB=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
OUTDIR="${WORKDIR}/output_wt_sublibrary_${SUB}"
DONEFILE="${OUTDIR}/.splitpipe_done"

echo "============================================================"
echo "Job started: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Host: $(hostname)"
echo "Working dir: $(pwd)"
echo "Base workdir: ${WORKDIR}"
echo "Sublibrary: ${SUB}"
echo "Threads: ${NTHREADS}"
echo "CPU affinity: $(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')"
echo "============================================================"

# Basic input validation
if [[ ! -d "${WORKDIR}" ]]; then
  echo "ERROR: WORKDIR not found: ${WORKDIR}"
  exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "ERROR: GENOME_DIR not found: ${GENOME_DIR}"
  exit 1
fi

if [[ ! -f "${SLT}" ]]; then
  echo "ERROR: Sample loading table not found: ${SLT}"
  exit 1
fi

fq1=( ${PREFIX}${SUB}_*_R1_001.fastq.gz )
fq2=( ${PREFIX}${SUB}_*_R2_001.fastq.gz )

if ((${#fq1[@]} == 0)); then
  echo "ERROR: No R1 FASTQs found for WT Sublibrary ${SUB}"
  exit 1
fi

if ((${#fq2[@]} == 0)); then
  echo "ERROR: No R2 FASTQs found for WT Sublibrary ${SUB}"
  exit 1
fi

if ((${#fq1[@]} != ${#fq2[@]})); then
  echo "ERROR: Mismatch between R1 and R2 FASTQ counts for WT Sublibrary ${SUB} (R1=${#fq1[@]}, R2=${#fq2[@]})"
  exit 1
fi

# Skip if already completed successfully
if [[ -f "${DONEFILE}" ]]; then
  echo "Sublibrary ${SUB} already completed. Skipping."
  exit 0
fi

# Protect against accidental reruns into a partial output directory
if [[ -d "${OUTDIR}" && ! -f "${DONEFILE}" ]]; then
  echo "ERROR: Output directory ${OUTDIR} already exists but is not marked complete."
  echo "Inspect it before rerunning. Remove it manually if you want to restart this sublibrary."
  exit 1
fi

echo "FASTQ R1 files (${#fq1[@]}):"
printf '  %s\n' "${fq1[@]}"

echo "FASTQ R2 files (${#fq2[@]}):"
printf '  %s\n' "${fq2[@]}"

echo "Starting split-pipe for WT Sublibrary ${SUB} at $(date)"

split-pipe \
  --mode all \
  --chemistry "${CHEM}" \
  --kit "${KIT}" \
  --fq1 "${fq1[@]}" \
  --fq2 "${fq2[@]}" \
  --genome_dir "${GENOME_DIR}" \
  --output_dir "${OUTDIR}" \
  --samp_sltab "${SLT}" \
  --nthreads "${NTHREADS}"

touch "${DONEFILE}"

echo "Finished WT Sublibrary ${SUB} at $(date)"
echo "Output: ${OUTDIR}"
```

Submit the SLURM job
```bash
cp "03a-BM082-R-Parse Biosciences Evercode WT Mega INTEGRA Sample Loading Table v2.1-260226 - Copy.xlsm" 03a-BM082-R-Parse_WT_Mega_SLT_v2.1-260226_Copy.xlsm
mkdir -p /home/users/jo912684/data/Seqmatic_032026_Parse_Mega_druggable/slurm_logs
sbatch BM082R_WT_splitpipe_array.sh
```

## ------------------------------------------------------------
## STEP 2: Merge WTA sublibraries into a combined matrix.
## ------------------------------------------------------------

Note: The regular "--mode comb" command crashed since it's trying to merge count matrices of 1.7 million (filtered) or ~14 million cells (unfiltered). Therefore, I crafted scripts to merge files of the sublibrary in a memory efficient manner, which are located in the 'utils' folder.

```bash
python3 utils/merge_dge_filtered.py \
  --sublibs \
    output_wt_sublibrary_01 \
    output_wt_sublibrary_02 \
    output_wt_sublibrary_03 \
    output_wt_sublibrary_04 \
    output_wt_sublibrary_05 \
    output_wt_sublibrary_06 \
    output_wt_sublibrary_07 \
    output_wt_sublibrary_08 \
    output_wt_sublibrary_09 \
    output_wt_sublibrary_10 \
    output_wt_sublibrary_11 \
    output_wt_sublibrary_12 \
    output_wt_sublibrary_13 \
    output_wt_sublibrary_14 \
    output_wt_sublibrary_15 \
    output_wt_sublibrary_16 \
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
- 7492 sgRNAs
- 2 sgRNAs

I dropped out the 2 sgRNAs moving forward to avoid issues.

There were the 2 dropped sgRNAs:
AKR1C3_P2_3,GGCTTTATATATCTTGTGGAAAGGACGAAACACCG,GTGCATAGGTGCCAAATCCC,GTTTAAGAGCTATGCTGGA,AKR1C3_P2
TUBB6_P1_4,GGCTTTATATATCTTGTGGAAAGGACGAAACACCG,GTGCACGATCTCCCTCATGG,GTTTAAGAGCTATGCTGGA,TUBB6_P1

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

CRISPR_PREFIX="BM082R-CRISPR-Sublibrary"
WT_PARENT_PREFIX="output_wt_sublibrary_"

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
--parent_dir output_wt_sublibrary_01 \
--output_dir output_combined_crispr \
--sublibraries output_crispr_sublibrary_01 output_crispr_sublibrary_02 output_crispr_sublibrary_03 output_crispr_sublibrary_04 output_crispr_sublibrary_05 output_crispr_sublibrary_06 output_crispr_sublibrary_07 output_crispr_sublibrary_08 output_crispr_sublibrary_09 output_crispr_sublibrary_10 output_crispr_sublibrary_11 output_crispr_sublibrary_12 output_crispr_sublibrary_13 output_crispr_sublibrary_14 output_crispr_sublibrary_15 output_crispr_sublibrary_16
```