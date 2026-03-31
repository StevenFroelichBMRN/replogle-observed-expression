# Perturb-seq run using just 1 sgRNA targeting CD81

# Location of working directory
```bash
cd /home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81
```

# Increase memory of the session
```{r}
library(future)
future::plan(future::sequential)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
```

# Load libraries
```{r}
library(Seurat)
library(Matrix)
library(ggplot2)
library(viridis)
library(tidyverse)
library(scCustomize)
library(ggsignif)
library(cowplot)
```

# Load data from WT samples
```{r}
# Load matrix
wt_expr <- readMM("data/Seqmatic_102025_Parse_10K_CD81_WT/output_combined/WT/DGE_filtered/count_matrix.mtx.gz")
wt_expr <- t(wt_expr)

# Load features (genes)
wt_features <- read.csv("data/Seqmatic_102025_Parse_10K_CD81_WT/output_combined/WT/DGE_filtered/all_genes.csv.gz", head=T)
wt_features$gene_name <- make.unique(wt_features$gene_name, sep="-dup") # make unique names
rownames(wt_expr) <- wt_features$gene_name

# Load cell barcodes
wt_barcodes <- read.csv("data/Seqmatic_102025_Parse_10K_CD81_WT/output_combined/WT/DGE_filtered/cell_metadata.csv.gz", head=T)
colnames(wt_expr) <- wt_barcodes$bc_wells

# Make Seurat object
wt_scobj <- CreateSeuratObject(counts= wt_expr, meta.data = wt_barcodes)
```

# Load data from 
```{r}
# Load matrix
cd81_expr <- readMM("data/Seqmatic_012926_Parse_Mini_CD81/output_combined/CD81/DGE_filtered/count_matrix.mtx")
cd81_expr <- t(cd81_expr)

# Load features (genes)
cd81_features <- read.csv("data/Seqmatic_012926_Parse_Mini_CD81/output_combined/CD81/DGE_filtered/all_genes.csv", head=T)
cd81_features$gene_name <- make.unique(cd81_features$gene_name, sep="-dup") # make unique names
rownames(cd81_expr) <- cd81_features$gene_name

# Load cell barcodes
cd81_barcodes <- read.csv("data/Seqmatic_012926_Parse_Mini_CD81/output_combined/CD81/DGE_filtered/cell_metadata.csv", head=T)
colnames(cd81_expr) <- cd81_barcodes$bc_wells

# Make Seurat object
cd81_scobj <- CreateSeuratObject(counts= cd81_expr, meta.data = cd81_barcodes)

# Add sgRNA assignment
md <- cd81_scobj@meta.data %>% rownames_to_column("cell")
ga <- read_csv("data/Seqmatic_012926_Parse_Mini_CD81/output_combined_crispr/CD81/guide_RNAs_filtered/guide_assignment.csv")
md_joined <- md %>% left_join(ga, by = "bc_wells")
rownames(md_joined) <- md_joined$cell
cd81_scobj@meta.data <- md_joined

# Any unassigned cell should have "none" in the "guide" column, indicating they didn't receive a sgRNA
cd81_scobj$guide <- as.character(cd81_scobj$guide)
cd81_scobj$guide[is.na(cd81_scobj$guide) | cd81_scobj$guide == ""] <- "none"
```

# Basic QC
```{r}
# Calculate percentage of mitochondrial genes:
wt_scobj[["percent.mt"]] <- PercentageFeatureSet(wt_scobj, pattern = "^MT-")
cd81_scobj[["percent.mt"]] <- PercentageFeatureSet(cd81_scobj, pattern = "^MT-")

# UMI filtering threshold:
lower_threshold= 1500

# Filtering:
wt_scobj <- subset(wt_scobj, subset = nFeature_RNA > lower_threshold &  percent.mt < 15)
cd81_scobj <- subset(cd81_scobj, subset = nFeature_RNA > lower_threshold &  percent.mt < 15)

# Counts per object:
# WT= 2,052 cells
# CD81= 14,715 cells
```


# Integrate both datasets
```{r}
# Add a "guide" column to the WT object
wt_scobj$guide <- "none"

# Define process to run sequentially to avoid memory issues
future::plan("sequential")

# Add a condition variable
wt_scobj$condition   <- "WT"
cd81_scobj$condition <- "CD81"

# Make cell names unique across objects
wt_scobj   <- RenameCells(wt_scobj,   add.cell.id = "WT")
cd81_scobj <- RenameCells(cd81_scobj, add.cell.id = "CD81")

# Standard Seurat integration workflow
obj_list <- list(WT = wt_scobj, CD81 = cd81_scobj)
obj_list <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x)
  x
})
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features)
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT")

# Save integrated object
saveRDS(integrated, file = "analyses/2026_PerturbSeq_Mini_CD81/integrated.rds")
```


# Standard dimensionality reduction
```{r}
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.5, verbose = FALSE)

# Enforce "guide" to have 2 levels:
integrated$guide <- as.character(integrated$guide)
integrated$guide[is.na(integrated$guide)] <- "none"
integrated$guide[!(integrated$guide %in% c("K21", "none"))] <- NA
integrated$guide[is.na(integrated$guide)] <- "none"
integrated$guide <- factor(integrated$guide, levels = c("K21", "none"))
```

# Evaluate knockdown of CD81 in cells that received the K21 sgRNA
```{r}
# Set assay to RNA and integrate counts for the two datasets
DefaultAssay(integrated) <- "RNA"
integrated <- JoinLayers(integrated, assay = "RNA")
integrated <- NormalizeData(integrated, assay = "RNA")
integrated <- ScaleData(integrated, assay = "RNA", features = rownames(integrated))

# Extract CD81 expression
cd81_expr <- FetchData(object = integrated, vars = "CD81")

# Extract sgRNA/guide status from metadata
sgRNA_status <- integrated@meta.data$guide
sgRNA_status <- as.character(sgRNA_status)
sgRNA_status[is.na(sgRNA_status) | sgRNA_status == ""] <- "none"

# Enforce only two levels: K21 and none (anything else -> none)
sgRNA_status[!(sgRNA_status %in% c("K21", "none"))] <- "none"
sgRNA_status <- factor(sgRNA_status, levels = c("none", "K21"))

# Combine into a dataframe (similar structure to your example)
df <- data.frame(
  CD81 = cd81_expr[["CD81"]],
  sgRNA_status = sgRNA_status
)

# Make boxplot comparing CD81 expression across groups
p_KD <- ggplot(df, aes(x = sgRNA_status, y = CD81+1, fill = sgRNA_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = 0.3) +
  labs(
    x = "",
    y = paste0("CD81 expression")
  ) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("none" = "#56B4E9", "K21" = "#D55E00")) +
  geom_signif(comparisons = list(c("none", "K21")), test = "wilcox.test") +
  theme(legend.position = "none")
ggsave("/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/plots/Plot_CD81_KD_test.pdf", p_KD, height = 5, width = 4)

# Statistical test
wilcox.test(CD81 ~ sgRNA_status, data = df) # W = 30586140, p-value < 2.2e-16

# Fold-change reduction
mean_none <- mean(df$CD81[df$sgRNA_status == "none"], na.rm = TRUE)
mean_K21  <- mean(df$CD81[df$sgRNA_status == "K21"],  na.rm = TRUE)

fold_change <- mean_K21 / mean_none # 0.03768917
(1 - fold_change) * 100 # 96.23 % redudction
```

# How many cells are neurons?
```{r}
integrated$combined <- "combined"
p_markers <- DotPlot_scCustom(seurat_object = integrated, group.by = "combined",
                 features = rev(c("GFAP","AQP4",
                              "ACHE",
                              "TH",
                              "CLDN5",
                              "PVALB", "SLC32A1", "DLX5", "DLX2", "GAD1", "GAD2", "DLX1", "SST", "SLC6A1", "ADARB2",
                              "SLC17A6", "NRP1", "SLC17A7", "NEUROD6",
                              "SYT1", "MAP2", "SYP", "DLG4", "MAPT","RBFOX3", "TUBB3", "DCX", "NCAM1",
                              "CD74",
                              "ISL1", "OLIG2",
                              "PAX6", "NES", "SOX2", "HES5",
                              "VCAN", "PDGFRA",
                              "MBP",
                              "POU5F1")), 
                 colors_use = viridis_plasma_dark_high,
                 flip_axes = T)
ggsave("/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/plots/QC_MarkersCheck.pdf", p_markers, width=4, height = 6)

```

# What other cell types are there?
```{r}
gene.markers <- FindAllMarkers(integrated, only.pos = T)
write.csv(gene.markers, file = "/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/gene.markers.csv")

gene.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

pdf("/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/plots/Heatmap_topmarkers.pdf", width=16, height=12)
DoHeatmap(integrated, features = top10$gene) + NoLegend()
dev.off()

pdf("/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/plots/UMAP_clusters_red0.5.pdf", width=5, height=5)
DimPlot_scCustom(seurat_object = integrated, repel = TRUE)
dev.off()

p_markers_clust <- DotPlot_scCustom(seurat_object = integrated, group.by = "seurat_clusters",
                 features = rev(c("GFAP","AQP4",
                              "ACHE",
                              "TH",
                              "CLDN5",
                              "PVALB", "SLC32A1", "DLX5", "DLX2", "GAD1", "GAD2", "DLX1", "SST", "SLC6A1", "ADARB2",
                              "SLC17A6", "NRP1", "SLC17A7", "NEUROD6",
                              "SYT1", "MAP2", "SYP", "DLG4", "MAPT","RBFOX3", "TUBB3", "DCX", "NCAM1",
                              "CD74",
                              "ISL1", "OLIG2",
                              "PAX6", "NES", "SOX2", "HES5",
                              "VCAN", "PDGFRA",
                              "MBP",
                              "POU5F1")), 
                 colors_use = viridis_plasma_dark_high,
                 flip_axes = T)
ggsave("/home/users/jo912684/analyses/2026_PerturbSeq_Mini_CD81/plots/QC_MarkersCheck_by_cluster.pdf", p_markers_clust, width=10, height = 6)


```