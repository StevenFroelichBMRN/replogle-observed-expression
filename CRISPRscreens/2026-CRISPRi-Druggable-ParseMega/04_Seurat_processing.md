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
scobj <- readRDS("analyses/2026_PerturbSeq_Mega_druggable/scobj.rds")
md <- scobj@meta.data %>% rownames_to_column("cell")

# Load sgRNA assignments from Parse's CRISPR pipeline
ga <- read_csv("data/Seqmatic_032026_Parse_Mega_druggable/output_combined_crispr/DRUG/guide_RNAs_filtered/guide_assignment.csv")

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
mean(table(scobj$gRNA[scobj$gRNA != "negative" & scobj$gRNA != "ambiguos"]))
# 75.33618 cells per sgRNA

# NTC targets
ntc_targets <- unique(grep("Non_Targeting_Human_CRi_", scobj$target, value = T))
# Cells with 'NTC' target: 7750
# Cells with 'negative' assignment: 388,782
# Cells with 'ambiguous': 204,304
```

## ------------------------------------------------------------
## STEP 3: QC
## ------------------------------------------------------------

### QC filterign will include:
- Remove cells with UMI counts below 1000
- More than 15% of MT reads
- 'Negative' sgRNA assignment (cells that did not receive sgRNAs)
- 'Ambiguous' sgRNA assignment (cells that received >1 sgRNAs)

Following same guidelines as Chardon et al 2024 (https://www.nature.com/articles/s41467-024-52490-4#Sec7)

```{r}
# Calculate percentage of mitochondrial genes:
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")

# UMI filtering thresholds:
lower_threshold= 1000

# Filtering:
scobj_filt <- subset(scobj, subset = nFeature_RNA > lower_threshold & 
percent.mt < 15 & 
gRNA != 'negative' & 
gRNA != 'ambiguous')
# Total cells retained= 337,668

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

ggsave("analyses/2026_PerturbSeq_Mega_druggable/QC_Feature_Counts.pdf", p_QC1, height=5, width = 5)
```

## ------------------------------------------------------------
## STEP 4: Basic metrics post-filtering
## ------------------------------------------------------------

```{r}
# Average number of cells for each sgRNA = 46.54921
mean(table(scobj_filt$gRNA))

# Number of unique sgRNAs detected = 7254
length(unique(scobj_filt$gRNA))

# Number of cells with NTC sgRNAs = 7,634
length(scobj_filt$target[scobj_filt$target == "NTC"])
```

## ------------------------------------------------------------
## STEP 5: Checkpoint to save RDS object
## ------------------------------------------------------------

```{r}
saveRDS(scobj_filt, file="analyses/2026_PerturbSeq_Mega_druggable/scobj_filt.rds")
```