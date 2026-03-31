# Use control data to determine the number of cells needed to detect X fold expression change for a specific gene of interest

## ------------------------------------------------------------
## STEP 0: Load proper environment
## ------------------------------------------------------------

```{r}
# Load conda environment
source /home/users/jo912684/miniconda/etc/profile.d/conda.sh
conda activate seurat5_env
```

## ------------------------------------------------------------
## STEP 1: Run script to perform calculations
## ------------------------------------------------------------

```{r}
Rscript utils/power_curve.R --seurat output_combined_filtered/scobj_filt.rds --genes APP,FMR1,MAPT,MECP2,SCN2A,SOD1,CDKL5,STXBP1,UBE3A,SYNGAP1,C9orf72 --target_col target --ctrl_label NTC --outdir output_combined_filtered --sim_reps 500 --sim_test wilcox
```