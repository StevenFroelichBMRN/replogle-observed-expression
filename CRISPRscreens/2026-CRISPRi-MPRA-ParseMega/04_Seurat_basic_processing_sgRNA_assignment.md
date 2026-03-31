# Process the Seurat object

## ------------------------------------------------------------
## STEP 0: Load libraries
## ------------------------------------------------------------

```{r}
library(Seurat) # make sure v5 is loaded
library(BPCells)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(edgeR)
library(Matrix)
library(tidyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
```

## ------------------------------------------------------------
## STEP 1: Assign sgRNAs to cells
## ------------------------------------------------------------

```{r}
# Grab the metadata
md <- scobj@meta.data %>% rownames_to_column("cell")

# Load sgRNA assignments from Parse's CRISPR pipeline
ga <- read_csv("output_combined_crispr/MPRA/guide_RNAs_filtered/guide_assignment.csv")

# Collapse duplicates per barcode
ga_clean <- ga %>%
  group_by(bc_wells) %>%
  summarize(
    gRNA = if (n_distinct(guide) == 1) first(guide) else "ambiguous",
    .groups = "drop"
  )

# Join on barcode: bc_wells_s (meta) == bc_wells (CSV)
md_joined <- md %>%
  left_join(ga_clean, by = c("bc_wells_s" = "bc_wells")) %>%
  mutate(
    # unmatched -> "negative"
    gRNA = coalesce(gRNA, "negative"),
    # target rules:
    target = case_when(
      gRNA %in% c("negative", "ambiguous") ~ gRNA,
      str_detect(gRNA, "NTC")              ~ "NTC",
      str_detect(gRNA, "_\\d+$")           ~ str_replace(gRNA, "_\\d+$", ""),
      TRUE                                 ~ gRNA
    )
  )

# Add new columns back to Seurat object (aligned by rownames)
rownames(md_joined) <- md_joined$cell
scobj <- Seurat::AddMetaData(scobj, metadata = md_joined[, c("gRNA", "target")])
```

## ------------------------------------------------------------
## STEP 2: Quick summary of sgRNA assignment
## ------------------------------------------------------------

```{r}
# Quick summary of sgRNA assignment
mean(table(scobj$gRNA[scobj$gRNA != "negative" & scobj$gRNA != "ambiguos" & scobj$gRNA != "NTC"]))
# 194.0198 cells per sgRNA

# Cells with 'negative' assignment: 866,405
# Cells with NTC controls: 21,904
# Cells with 'ambiguous': 275,815
```

## ------------------------------------------------------------
## STEP 3: QC
## ------------------------------------------------------------

### QC filterign will include:
- Remove cells with UMI counts below 1500 or above 4 SD from mean
- More than 15% of MT reads
- 'Negative' sgRNA assignment (cells that did not receive sgRNAs)
- 'Ambiguous' sgRNA assignment (cells that received >1 sgRNAs)

Following same guidelines as Chardon et al 2024 (https://www.nature.com/articles/s41467-024-52490-4#Sec7)

```{r}
# Calculate percentage of mitochondrial genes:
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")

# UMI filtering thresholds:
upper_threshold= mean(scobj$nFeature_RNA) + 4*sd(scobj$nFeature_RNA) # 6399.926 UMIs
lower_threshold= 1500

# Filtering:
scobj_filt <- subset(scobj, subset = nFeature_RNA > lower_threshold & 
nFeature_RNA < upper_threshold & 
percent.mt < 15 & 
gRNA != 'negative' & 
gRNA != 'ambiguous')
# Total cells retained= 537,915

# Plots post-filtering:
scobj_filt$orig.ident <- scobj_filt$sample
Idents(scobj_filt) <- scobj_filt$orig.ident

p_QC1 <- FetchData(scobj_filt, vars = c("nCount_RNA", "nFeature_RNA")) %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(alpha = 0.3, size = 0.3) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "nCount_RNA (UMIs)", y = "nFeature_RNA (genes)") +
  theme_classic()

ggsave("output_combined_filtered/plots/QC_Feature_Counts.pdf", p_QC1, height=5, width = 5)
```
## ------------------------------------------------------------
## STEP 4: Basic metrics post-filtering
## ------------------------------------------------------------

```{r}
# Average number of cells for each sgRNA = 116.005
mean(table(scobj_filt$gRNA))

# Number of unique sgRNAs detected = 4637
length(unique(scobj_filt$gRNA))

# Number of cells with NTC sgRNAs = 18,930
length(scobj_filt$target[scobj_filt$target == "NTC"])
```

## ------------------------------------------------------------
## STEP 5: Checkpoint to save RDS object
## ------------------------------------------------------------

```{r}
saveRDS(scobj_filt, file="output_combined_filtered/scobj_filt.rds")
```