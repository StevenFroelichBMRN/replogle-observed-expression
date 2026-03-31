# Load H5 file into Seurat for further processing.

Need to use Seurat v5 to use the BPCells method to load large data matrices

## ------------------------------------------------------------
## STEP 0: Start and R session with Seurat v5
## ------------------------------------------------------------

## Option 1: activate conda environment
```bash
source miniconda3/etc/profile.d/conda.sh
conda activate seurat5_env
R # start R
```

## Option 2: run the rstudio node (https://hpcprod.bmrn.com/rstudio) and make sure to select 'R 4.5.0'

## Load libraries
```{r}
library(Seurat) # make sure v5 is loaded
library(BPCells)
```

## ------------------------------------------------------------
## STEP 1: Convert from H5 to Seurat object.
## ------------------------------------------------------------

```{r}
# Read H5 file from the WTA matrix.
wta.data <- open_matrix_10x_hdf5(path = "data/Seqmatic_032026_Parse_Mega_druggable/output_combined_filtered/count_matrix.h5")

# Write the matrix to a directory
write_matrix_dir(mat = wta.data, dir = 'data/Seqmatic_032026_Parse_Mega_druggable/output_combined_filtered/bpcells')

# Now that we have the matrix on disk, we can load it
wta.mat <- BPCells::open_matrix_dir(dir = "data/Seqmatic_032026_Parse_Mega_druggable/output_combined_filtered/bpcells/")
metadata <- read.csv("data/Seqmatic_032026_Parse_Mega_druggable/output_combined_filtered/cell_metadata_sub.csv", head=T)
rownames(metadata) <- metadata$bc_wells_s

# Create Seurat Object
scobj <- CreateSeuratObject(counts = wta.mat, meta.data = metadata)

# Save file
saveRDS(scobj, file="analyses/2026_PerturbSeq_Mega_druggable/scobj.rds")
```