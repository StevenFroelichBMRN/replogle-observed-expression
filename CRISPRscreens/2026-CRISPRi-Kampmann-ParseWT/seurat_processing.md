# Perturb-seq run using the sgRNA from the published study Tian et al 2020 from the Kampmann UCSF lab

# Location of working directory
```bash
cd /home/users/jo912684/analyses/2026_PerturbSeq_WT_Kampmann
```

# Increase memory of the session
```{r}
library(future)
future::plan(future::sequential)
options(future.globals.maxSize = 10 * 1024^3)
```

# Load libraries
```{r}
library(Seurat)
library(Matrix)
library(anndataR)
library(ggplot2)
library(viridis)
library(tidyverse)
library(data.table)
library(harmony)
library(ggsignif)
library(cowplot)
library(zellkonverter)
library(SingleCellExperiment)
library(readr)
library(purrr)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
```

# Load WTA data
```{r}
# Load matrix
expr_mat <- readMM("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined/K1/DGE_filtered/count_matrix.mtx")
expr_mat <- t(expr_mat)

# Load features (genes)
features <- read.csv("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined/K1/DGE_filtered/all_genes.csv", head=T)
features$gene_name <- make.unique(features$gene_name, sep="-dup") # make unique names
rownames(expr_mat) <- features$gene_name

# Load cell barcodes
barcodes <- read.csv("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined/K1/DGE_filtered/cell_metadata.csv", head=T)
colnames(expr_mat) <- barcodes$bc_wells
```

# Load CRISPR data
```{r}
# Load matrix
crispr_matrix <- readMM("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined_crispr/K1/guide_RNAs_unfiltered/count_matrix.mtx")
crispr_matrix <- t(crispr_matrix)

# Load features (sgRNAs)
crispr_features <- read.csv("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined_crispr/K1/guide_RNAs_unfiltered/all_guides.csv", head=T)
crispr_features$gene_name <- make.unique(crispr_features$gene_name, sep="-dup") # make unique names
rownames(crispr_matrix) <- crispr_features$gene_name

# Load cell barcodes
crispr_barcodes <- read.csv("data/Seqmatic_112025_Parse_100K_Kampmann/newvolume/analysis/combined_crispr/K1/guide_RNAs_unfiltered/cell_metadata.csv", head=T)
colnames(crispr_matrix) <- crispr_barcodes$bc_wells
```

# Assign sgRNA to each cell (needed to do this manually due to very low sgRNA call - not ideal)
```{r}
# Sanity: column names must match exactly (order and values)
if (!identical(colnames(crispr_matrix), colnames(expr_mat))) {
  stop("Column names of crispr_matrix and expression_matrix do not match. 
       Align colnames before assigning.")
}

min_count_per_cell <- 1    # minimum sgRNA count required to consider a guide present
tie_method <- "none"    # choose "none", "first", "random"

## ---- Compute per-cell top sgRNA efficiently ----
# dgTMatrix stores triplets i (row), j (col), x (value).
# We'll convert to a data.table of non-zero entries, then pick the max x per j (column).

# Extract triplets
trip <- data.table(
  i = crispr_matrix@i + 1L,            # 1-based row index
  j = crispr_matrix@j + 1L,            # 1-based column index
  x = crispr_matrix@x                  # counts
)

# Filter by minimum count if desired
if (min_count_per_cell > 1) {
  trip <- trip[x >= min_count_per_cell]
}

# If there are no non-zero entries after filtering, return all NA
if (nrow(trip) == 0L) {
  assigned_sgRNA <- rep(NA_character_, ncol(crispr_matrix))
  names(assigned_sgRNA) <- colnames(crispr_matrix)
} else {
  # For each column j, find the maximum x; keep all rows that equal the column max (to resolve ties later)
  # Step 1: compute max x per j
  col_max <- trip[, .(max_x = max(x)), by = j]

  # Step 2: join back to keep only maxima rows (potentially multiple per j)
  trip_max <- trip[col_max, on = "j"][x == max_x]

  # Resolve ties based on tie_method
  if (tie_method == "first") {
    # Keep the first row within each j (data.table's fast head-by-group)
    trip_max_unique <- trip_max[, .SD[1L], by = j]
  } else if (tie_method == "random") {
    set.seed(1) # for reproducibility; change or remove as needed
    trip_max_unique <- trip_max[, .SD[sample.int(.N, 1L)], by = j]
  } else if (tie_method == "none") {
    # Mark all tied maxima as NA by collapsing to a single NA record per j
    # Identify j's with >1 rows (ties)
    tie_info <- trip_max[, .N, by = j]
    ties <- tie_info[N > 1L, j]
    noties <- tie_info[N == 1L, j]

    # Unique pick for non-ties
    trip_max_unique <- trip_max[j %in% noties, .SD[1L], by = j]

    # Build NA records for ties
    if (length(ties)) {
      na_rows <- data.table(j = ties, i = NA_integer_, x = NA_real_)
      trip_max_unique <- rbind(trip_max_unique, na_rows, fill = TRUE)
    }
  } else {
    stop("Unknown tie_method. Use 'first', 'random', or 'none'.")
  }

  # Map row index i -> sgRNA name
  sgRNA_names <- rownames(crispr_matrix)
  # Initialize vector with NA
  assigned_sgRNA <- rep(NA_character_, ncol(crispr_matrix))
  names(assigned_sgRNA) <- colnames(crispr_matrix)

  # Fill assignments for columns that had at least one entry
  valid_rows <- !is.na(trip_max_unique$i)
  assigned_sgRNA[trip_max_unique$j[valid_rows]] <- sgRNA_names[trip_max_unique$i[valid_rows]]
}

## ---- Build a tidy result and a vector aligned to expression_matrix ----
assignment_df <- data.frame(
  cell_barcode    = colnames(crispr_matrix),   # same as expression_matrix colnames
  assigned_sgRNA  = assigned_sgRNA,
  stringsAsFactors = FALSE
)

# Optional: also return the counts of the assigned sgRNA per cell (useful for QC)
# We can fetch the count by looking it up in the sparse matrix for each (i,j).
# For speed, build a lookup on triplets that correspond to the chosen (i,j).
# Create a keyed data.table for fast joins
trip_keyed <- copy(trip)
setkey(trip_keyed, j, i)

# Build (j,i) pairs for assigned cells
assigned_pairs <- data.table(
  j = seq_len(ncol(crispr_matrix)),
  i = match(assigned_sgRNA, rownames(crispr_matrix))
)

# Lookup counts (will be NA for cells with no assignment)
assigned_counts <- trip_keyed[assigned_pairs, x]

assignment_df$assigned_count <- assigned_counts

# Add the targeted gene information
guides_info <- read.csv("data/Seqmatic_112025_Parse_100K_Kampmann/guides.csv", head=F)
sgRNA_to_gene <- setNames(guides_info$V5, guides_info$V1)

# Map targeted gene for the assigned sgRNA (NA if no assignment or sgRNA not in the lookup)
assignment_df$targeted_gene <- sgRNA_to_gene[assignment_df$assigned_sgRNA]

# A convenience flag for non-targeting guides
assignment_df$is_non_targeting <- !is.na(assignment_df$targeted_gene) &
                                  grepl("^Non_Targeting", assignment_df$targeted_gene)

```

# Create Seurat object and QC
```{r, dev = "RStudioGD"}
scobj <- CreateSeuratObject(counts= expr_mat, meta.data = barcodes)

# Add the sgRNA assignment
scobj@meta.data$sgRNA <- assignment_df$assigned_sgRNA
scobj@meta.data$sgRNA_count <- assignment_df$assigned_count
scobj@meta.data$targeted_gene <- assignment_df$targeted_gene
scobj@meta.data$is_non_targeting <- assignment_df$is_non_targeting

# Calculate percentage of mitochondrial genes:
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")

# Filtering data:
scobj <- subset(scobj, subset = nFeature_RNA > 200 & percent.mt < 15)
# Final count: 127,858 cells

# Add 'sample' code
scobj$sample <- "K1"

# Save object
saveRDS(scobj, file="analyses/2026_PerturbSeq_WT_Kampmann/scobj.rds")
```

# Compare our data to the original Tian 2020 dataset
```{r}

# -----------------------------
# Load and preprocess both datasets
# -----------------------------

# Load the raw Tian 2020 data and normalize
tian <- read_h5ad("data/public_scRNAseq_datasets/Kampmann_lab/TianKampmann2021_CRISPRi.h5ad", as = "Seurat")
rna <- tian[["RNA"]] # H5AD is loaded with counts as "X", so need to correct it
LayerData(rna, layer = "counts") <- LayerData(rna, layer = "X")
DefaultLayer(rna) <- "counts"
tian[["RNA"]] <- rna
Layers(tian[["RNA"]])
DefaultLayer(tian[["RNA"]])

tian <- NormalizeData(tian)
tian <- FindVariableFeatures(tian, selection.method = "vst", nfeatures = 2000)

# Normalize our single-cell data for a fair comparison.
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)


# -----------------------------
# Set output directories
# -----------------------------

out_base <- "analyses/2026_PerturbSeq_WT_Kampmann/compare_to_Tian2020_signature"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

dir_sig   <- file.path(out_base, "01_signature_space")
dir_cor   <- file.path(out_base, "02_deg_correlation")
dir_volc  <- file.path(out_base, "03_volcano")
dir_de    <- file.path(out_base, "04_de_tables")

dir.create(dir_sig,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_cor,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_volc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_de,   showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Parameters
# -----------------------------

set.seed(1)

# DE settings (within each dataset only, perturbed vs that dataset's controls)
de_test     <- "wilcox"
min_pct     <- 0.10
logfc_thr   <- 0.10
only_pos    <- FALSE

# How many genes to keep per (dataset, target) signature when building the shared signature space
# (top by absolute logFC among DE results)
top_genes_per_sig <- 500

# PCA for signature-space distances
sig_pcs_for_distance <- 20

# Permutation test (shuffle Tian target labels across Tian signatures)
n_perm <- 1000

# -----------------------------
# Harmonize metadata
# -----------------------------
DefaultAssay(scobj) <- "RNA"
DefaultAssay(tian)  <- "RNA"

sc_target_col   <- "targeted_gene"
sc_ctrl_label   <- "Non_Targeting_Human_CRi"
tian_target_col <- "perturbation"
tian_ctrl_label <- "control"

scobj$dataset <- "BMRN"
tian$dataset  <- "Tian2020"

scobj$target_gene <- as.character(scobj@meta.data[[sc_target_col]])
tian$target_gene  <- as.character(tian@meta.data[[tian_target_col]])

scobj$is_control <- scobj$target_gene == sc_ctrl_label
tian$is_control  <- tian$target_gene  == tian_ctrl_label

scobj$target_gene[is.na(scobj$target_gene) | scobj$target_gene==""] <- NA
tian$target_gene[is.na(tian$target_gene) | tian$target_gene==""] <- NA

# Align gene space to intersection (so logFCs are comparable by gene)
common_genes <- intersect(rownames(scobj), rownames(tian))
scobj2 <- subset(scobj, features = common_genes)
tian2  <- subset(tian,  features = common_genes)

# Normalize each dataset (kept simple; signatures are computed within-dataset)
scobj2 <- NormalizeData(scobj2, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
tian2  <- NormalizeData(tian2,  normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

# Set identities for DE (target vs CONTROL) within each dataset
Idents(scobj2) <- ifelse(scobj2$is_control, "CONTROL", scobj2$target_gene)
Idents(tian2)  <- ifelse(tian2$is_control,  "CONTROL", tian2$target_gene)

# Determine common targets (excluding controls)
targets_bmrn <- scobj2@meta.data %>%
  filter(!is.na(target_gene), !is_control) %>%
  pull(target_gene) %>% unique()

targets_tian <- tian2@meta.data %>%
  filter(!is.na(target_gene), !is_control) %>%
  pull(target_gene) %>% unique()

targets_both <- intersect(targets_bmrn, targets_tian)
targets_both <- sort(targets_both)

# 173 common targetted genes

# -----------------------------
# Run DE per target within each dataset, save tables + volcano plots
# Also build per-(dataset,target) signature vectors (logFC) for later PCA
# -----------------------------
# We will store DE results in lists keyed by target
de_bmrn_list <- list()
de_tian_list <- list()

# And store signatures as named numeric vectors (gene -> logFC)
sig_bmrn <- list()
sig_tian <- list()

for (tg in targets_both) {
  message("DE: ", tg)

  # Basic cell count checks (avoid unstable DE)
  n_b <- sum(Idents(scobj2) == tg, na.rm = TRUE)
  n_bc <- sum(Idents(scobj2) == "CONTROL", na.rm = TRUE)
  n_t <- sum(Idents(tian2) == tg, na.rm = TRUE)
  n_tc <- sum(Idents(tian2) == "CONTROL", na.rm = TRUE)

  if (n_b < 30 || n_bc < 30 || n_t < 30 || n_tc < 30) {
    message("  Skipping ", tg, " (too few cells). BMRN:", n_b, " ctrl:", n_bc, "  Tian:", n_t, " ctrl:", n_tc)
    next
  }

  # --- BMRN DE ---
  de_b <- FindMarkers(
    scobj2,
    ident.1 = tg, ident.2 = "CONTROL",
    test.use = de_test,
    min.pct = min_pct,
    logfc.threshold = logfc_thr,
    only.pos = only_pos,
    verbose = FALSE
  )
  de_b$gene <- rownames(de_b)

  # --- Tian DE ---
  de_t <- FindMarkers(
    tian2,
    ident.1 = tg, ident.2 = "CONTROL",
    test.use = de_test,
    min.pct = min_pct,
    logfc.threshold = logfc_thr,
    only.pos = only_pos,
    verbose = FALSE
  )
  de_t$gene <- rownames(de_t)

  # Decide logFC column names robustly (Seurat can output avg_log2FC or avg_logFC)
  b_fc_col <- if ("avg_log2FC" %in% colnames(de_b)) "avg_log2FC" else "avg_logFC"
  t_fc_col <- if ("avg_log2FC" %in% colnames(de_t)) "avg_log2FC" else "avg_logFC"

  # Save DE tables
  tg_de_dir <- file.path(dir_de, tg)
  dir.create(tg_de_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(as.data.table(de_b), file.path(tg_de_dir, "DE_BMRN.tsv"), sep = "\t")
  fwrite(as.data.table(de_t), file.path(tg_de_dir, "DE_Tian2020.tsv"), sep = "\t")

  de_bmrn_list[[tg]] <- de_b
  de_tian_list[[tg]] <- de_t

  # Volcano plots (side-by-side)
  b_pcol <- if ("p_val_adj" %in% colnames(de_b)) "p_val_adj" else "p_val"
  t_pcol <- if ("p_val_adj" %in% colnames(de_t)) "p_val_adj" else "p_val"

  vb <- de_b %>% mutate(dataset = "BMRN",
                        logFC = .data[[b_fc_col]],
                        neglog10p = -log10(pmax(.data[[b_pcol]], 1e-300)))
  vt <- de_t %>% mutate(dataset = "Tian2020",
                        logFC = .data[[t_fc_col]],
                        neglog10p = -log10(pmax(.data[[t_pcol]], 1e-300)))

  vboth <- bind_rows(vb, vt)

  p_volc <- ggplot(vboth, aes(x = logFC, y = neglog10p)) +
    geom_point(alpha = 0.5, size = 0.8) +
    facet_wrap(~dataset, nrow = 1, scales = "free_y") +
    theme_classic() +
    ggtitle(paste0("Volcano: ", tg, " (perturbed vs within-dataset control)")) +
    xlab("avg logFC") + ylab("-log10(p)")

  tg_volc_dir <- file.path(dir_volc, tg)
  dir.create(tg_volc_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(tg_volc_dir, "volcano.png"), p_volc, width = 10, height = 4.5, dpi = 200)

  # Build signatures: take TOP genes by |logFC| (from DE results)
  # (This keeps the signature-space matrix compact and focused on perturbation signal)
  de_b2 <- de_b %>% mutate(absFC = abs(.data[[b_fc_col]])) %>% arrange(desc(absFC))
  de_t2 <- de_t %>% mutate(absFC = abs(.data[[t_fc_col]])) %>% arrange(desc(absFC))

  b_keep <- head(de_b2$gene, top_genes_per_sig)
  t_keep <- head(de_t2$gene, top_genes_per_sig)

  # Use union of top genes from both datasets for this target (more fair for cross-study)
  keep_union <- union(b_keep, t_keep)

  b_vec <- de_b %>% filter(gene %in% keep_union) %>% select(gene, !!b_fc_col)
  t_vec <- de_t %>% filter(gene %in% keep_union) %>% select(gene, !!t_fc_col)

  b_named <- b_vec[[b_fc_col]]; names(b_named) <- b_vec$gene
  t_named <- t_vec[[t_fc_col]]; names(t_named) <- t_vec$gene

  # Store
  sig_bmrn[[tg]] <- b_named
  sig_tian[[tg]] <- t_named
}

# Filter targets that actually produced DE results in both
valid_targets <- intersect(names(sig_bmrn), names(sig_tian))
valid_targets <- sort(valid_targets)

# 167 valid targets retained (showed DE in both datasets)

# -----------------------------
# Build a shared signature matrix: genes x (dataset_target)
# Values = logFC; missing genes filled with 0
# -----------------------------
all_sig_genes <- unique(c(unlist(lapply(sig_bmrn[valid_targets], names)),
                          unlist(lapply(sig_tian[valid_targets], names))))
all_sig_genes <- sort(all_sig_genes)

col_names <- c(paste0("BMRN__", valid_targets), paste0("Tian2020__", valid_targets))

sig_mat <- matrix(0, nrow = length(all_sig_genes), ncol = length(col_names),
                  dimnames = list(all_sig_genes, col_names))

# Fill BMRN columns
for (tg in valid_targets) {
  v <- sig_bmrn[[tg]]
  sig_mat[names(v), paste0("BMRN__", tg)] <- as.numeric(v)
}

# Fill Tian columns
for (tg in valid_targets) {
  v <- sig_tian[[tg]]
  sig_mat[names(v), paste0("Tian2020__", tg)] <- as.numeric(v)
}

# -----------------------------
# Signature-space PCA: points are (dataset × target) signatures
# -----------------------------
# Center/scale across signatures (columns are samples)
pca <- prcomp(t(sig_mat), center = TRUE, scale. = TRUE)

pca_df <- data.frame(pca$x[, 1:2, drop = FALSE])
pca_df$signature <- rownames(pca_df)
pca_df <- pca_df %>%
  tidyr::separate(signature, into = c("dataset", "target_gene"), sep = "__", remove = FALSE)

p_sig <- ggplot(pca_df, aes(x = PC1, y = PC2, color = target_gene, shape = dataset)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic() +
  ggtitle("Perturbation signature space (PCA on DE logFC vectors)\nOne point per (dataset × target)") +
  labs(color = "Target gene", shape = "Dataset")

ggsave(file.path(dir_sig, "PCA_signatures.png"), p_sig, width = 11, height = 7, dpi = 200)

# Also save a dataset-only version (often cleaner visually)
p_sig_ds <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic() +
  ggtitle("Perturbation signature space (PCA)\nColor = dataset") +
  labs(color = "Dataset")

ggsave(file.path(dir_sig, "PCA_signatures_dataset_only.png"), p_sig_ds, width = 8, height = 6, dpi = 200)

# Save PCA coordinates
fwrite(pca_df, file.path(dir_sig, "signature_pca_coordinates.tsv"), sep = "\t")

# -----------------------------
# Quantitative similarity: distance between matched targets in PCA space,
# compared to a null by shuffling Tian2020 target labels across Tian signatures.
# -----------------------------
kpcs <- min(sig_pcs_for_distance, ncol(pca$x))
pcs_use <- paste0("PC", seq_len(kpcs))

pc_all <- as.data.frame(pca$x[, seq_len(kpcs), drop = FALSE])
pc_all$signature <- rownames(pc_all)
pc_all <- pc_all %>% tidyr::separate(signature, into = c("dataset", "target_gene"), sep = "__", remove = FALSE)

bmrn_pc <- pc_all %>% filter(dataset == "BMRN", target_gene %in% valid_targets)
tian_pc <- pc_all %>% filter(dataset == "Tian2020", target_gene %in% valid_targets)

# observed distances
obs <- lapply(valid_targets, function(tg) {
  a <- bmrn_pc %>% filter(target_gene == tg) %>% select(all_of(pcs_use)) %>% as.numeric()
  b <- tian_pc %>% filter(target_gene == tg) %>% select(all_of(pcs_use)) %>% as.numeric()
  data.frame(target_gene = tg, dist = sqrt(sum((a - b)^2)))
}) %>% bind_rows()

# permutation: shuffle Tian target labels across Tian signatures
null_mat <- matrix(NA_real_, nrow = length(valid_targets), ncol = n_perm,
                   dimnames = list(valid_targets, NULL))

tian_pc_mat <- as.matrix(tian_pc[, pcs_use, drop = FALSE])
rownames(tian_pc_mat) <- tian_pc$target_gene

bmrn_pc_mat <- as.matrix(bmrn_pc[, pcs_use, drop = FALSE])
rownames(bmrn_pc_mat) <- bmrn_pc$target_gene

message("Permutation test on signature distances: n_perm=", n_perm)
for (b in seq_len(n_perm)) {
  perm_labels <- sample(valid_targets, replace = FALSE)  # permuted mapping for Tian rows
  # permuted Tian matrix: rownames become perm_labels (same order as valid_targets)
  perm_tian <- tian_pc_mat
  rownames(perm_tian) <- perm_labels

  for (tg in valid_targets) {
    a <- bmrn_pc_mat[tg, ]
    bb <- perm_tian[tg, ]  # after renaming rows, tg now points to a random Tian signature
    null_mat[tg, b] <- sqrt(sum((a - bb)^2))
  }
}

res_dist <- obs %>%
  rowwise() %>%
  mutate(
    null = list(null_mat[target_gene, ]),
    p_perm = (sum(unlist(null) <= dist) + 1) / (length(unlist(null)) + 1),
    z = (dist - mean(unlist(null))) / sd(unlist(null))
  ) %>%
  ungroup() %>%
  arrange(p_perm, dist)

fwrite(res_dist, file.path(dir_sig, "signature_distance_per_target.tsv"), sep = "\t")

p_dist <- ggplot(res_dist, aes(x = dist)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  ggtitle("Observed distances between matched targets in signature PCA space\n(smaller = more similar)")

ggsave(file.path(dir_sig, "signature_distance_hist.png"), p_dist, width = 8, height = 5, dpi = 200)

# Ensure consistent ordering
targets_order <- valid_targets
obs2 <- obs %>% arrange(match(target_gene, targets_order))

# Compute per-gene empirical p-values (one-sided: observed <= null)
per_gene_p <- sapply(targets_order, function(tg) {
  nullv <- null_mat[tg, ]
  nullv <- nullv[!is.na(nullv)]
  if (length(nullv) == 0) return(NA_real_)
  (sum(nullv <= obs2$dist[obs2$target_gene == tg]) + 1) / (length(nullv) + 1)
})
names(per_gene_p) <- targets_order

# Compute z-scores per gene (obs - mean(null)) / sd(null)
per_gene_z <- sapply(targets_order, function(tg) {
  nullv <- null_mat[tg, ]
  nullv <- nullv[!is.na(nullv)]
  if (length(nullv) <= 1) return(NA_real_)
  (obs2$dist[obs2$target_gene == tg] - mean(nullv)) / sd(nullv)
})
names(per_gene_z) <- targets_order

# BH adjustment across genes
per_gene_p_adj <- p.adjust(per_gene_p, method = "BH")

# Compose results table
res_table <- obs2 %>%
  mutate(p_perm = per_gene_p[match(target_gene, names(per_gene_p))],
         p_adj  = per_gene_p_adj[match(target_gene, names(per_gene_p_adj))],
         zscore = per_gene_z[match(target_gene, names(per_gene_z))])

# Write per-target table
fwrite(res_table, file.path(dir_sig, "signature_distance_per_target_with_pvals.tsv"), sep = "\t")

# -----------------------------
# Global test: is mean(observed distances) smaller than expected?
# For each permutation, compute the mean of permuted distances across targets,
# then compare observed mean to that null distribution.
# -----------------------------
mean_nulls <- apply(null_mat, 2, function(col) mean(col, na.rm = TRUE))
obs_mean <- mean(res_table$dist, na.rm = TRUE)
global_p <- (sum(mean_nulls <= obs_mean, na.rm = TRUE) + 1) / (length(mean_nulls[!is.na(mean_nulls)]) + 1)

global_summary <- data.frame(
  obs_mean_distance = obs_mean,
  null_mean_mean = mean(mean_nulls, na.rm = TRUE),
  null_mean_sd   = sd(mean_nulls, na.rm = TRUE),
  global_p_empirical = global_p,
  n_perm = length(mean_nulls)
)
fwrite(global_summary, file.path(dir_sig, "signature_distance_global_test.tsv"), sep = "\t")

# Plot global null distribution of mean distances with observed mean
p_global <- ggplot(data.frame(mean_nulls = mean_nulls), aes(x = mean_nulls)) +
  geom_histogram(bins = 40, fill = "grey80", color = "grey30") +
  geom_vline(xintercept = obs_mean, color = "blue", size = 1) +
  theme_bw() +
  ggtitle(sprintf("Global permuted-mean distances (n_perm=%d)\nobserved mean = %.4f, p_emp = %.4g",
                  length(mean_nulls), obs_mean, global_p)) +
  xlab("Mean distance across targets (permutation)") +
  ylab("Count")

ggsave(file.path(dir_sig, "global_mean_distance_perm_hist.png"), p_global, width = 7, height = 5, dpi = 200)

# -----------------------------
# Per-gene null histograms: show null distribution + observed distance vertical line
# Saves one PNG per target in subfolder; skip targets with too few null samples
# -----------------------------
per_gene_plot_dir <- file.path(dir_sig, "per_target_null_plots")
dir.create(per_gene_plot_dir, showWarnings = FALSE, recursive = TRUE)

for (tg in targets_order) {
  nullv <- null_mat[tg, ]
  nullv <- nullv[!is.na(nullv)]
  if (length(nullv) < 20) next
  obs_val <- res_table$dist[res_table$target_gene == tg]

  dfnull <- data.frame(null_dist = nullv)
  p <- ggplot(dfnull, aes(x = null_dist)) +
    geom_histogram(bins = 30, fill = "grey80", color = "grey40") +
    geom_vline(xintercept = obs_val, color = "red", size = 1) +
    theme_bw() +
    ggtitle(sprintf("%s — observed=%.4f, p_perm=%.4g, p_adj=%.4g",
                    tg,
                    obs_val,
                    res_table$p_perm[res_table$target_gene==tg],
                    res_table$p_adj[res_table$target_gene==tg])) +
    xlab("Permuted distances") + ylab("Count")

  ggsave(file.path(per_gene_plot_dir, paste0(tg, "_null_hist.png")), p, width = 6, height = 4, dpi = 200)
}

# -----------------------------
# Make a volcano-like ranking plot: distance vs -log10(p_adj)
# -----------------------------
res_table <- res_table %>%
  mutate(rank_score = -log10(p_adj + 1e-300) * ( -1 * sign(dist - median(dist, na.rm = TRUE)) )) # smaller dist = positive score

p_rank <- ggplot(res_table, aes(x = dist, y = -log10(p_adj + 1e-300))) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  ggtitle("Target ranking: distance (x) vs -log10(p_adj) (y)") +
  xlab("Observed distance (signature PCA)") +
  ylab("-log10(p_adj)")

ggsave(file.path(dir_sig, "distance_vs_neglog10padj.png"), p_rank, width = 7, height = 5, dpi = 200)

# -----------------------------
# DEG logFC correlation per target (BMRN vs Tian2020)
# Use overlap of DE gene sets (after thresholds), plus a global histogram
# -----------------------------
cor_tbl <- list()

for (tg in valid_targets) {
  de_b <- de_bmrn_list[[tg]]
  de_t <- de_tian_list[[tg]]

  # robustly pick FC columns
b_fc_col <- if ("avg_log2FC" %in% colnames(de_b)) "avg_log2FC" else "avg_logFC"
t_fc_col <- if ("avg_log2FC" %in% colnames(de_t)) "avg_log2FC" else "avg_logFC"

mrg <- inner_join(
  de_b %>% dplyr::select(gene, logFC_BMRN    = all_of(b_fc_col)),
  de_t %>% dplyr::select(gene, logFC_Tian2020 = all_of(t_fc_col)),
  by = "gene"
)

if (nrow(mrg) < 20) next

# force numeric (guards against weird classes)
mrg$logFC_BMRN     <- as.numeric(mrg$logFC_BMRN)
mrg$logFC_Tian2020 <- as.numeric(mrg$logFC_Tian2020)

# remove any NAs/infs (cor() will otherwise return NA or error in edge cases)
mrg <- mrg[is.finite(mrg$logFC_BMRN) & is.finite(mrg$logFC_Tian2020), , drop = FALSE]
if (nrow(mrg) < 20) next

pear  <- suppressWarnings(cor(mrg$logFC_BMRN, mrg$logFC_Tian2020, method = "pearson"))
spear <- suppressWarnings(cor(mrg$logFC_BMRN, mrg$logFC_Tian2020, method = "spearman"))

  cor_tbl[[tg]] <- data.frame(target_gene = tg, n_overlap = nrow(mrg), pearson = pear, spearman = spear)

  out_tg <- file.path(dir_cor, tg)
  dir.create(out_tg, showWarnings = FALSE, recursive = TRUE)

  fwrite(as.data.table(mrg), file.path(out_tg, "logFC_overlap.tsv"), sep = "\t")

  p_cor <- ggplot(mrg, aes(x = logFC_BMRN, y = logFC_Tian2020)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    ggtitle(paste0("DE logFC correlation: ", tg,
                   "\nPearson=", round(pear, 3),
                   "  Spearman=", round(spear, 3),
                   "  (n=", nrow(mrg), ")")) +
    xlab("BMRN: avg logFC (perturbed vs BMRN control)") +
    ylab("Tian2020: avg logFC (perturbed vs Tian control)")

  ggsave(file.path(out_tg, "DE_logFC_correlation.png"), p_cor, width = 7, height = 6, dpi = 200)
}

cor_df <- bind_rows(cor_tbl)
if (nrow(cor_df) > 0) {
  fwrite(cor_df, file.path(dir_cor, "DE_correlation_summary.tsv"), sep = "\t")

  p_hist_cor <- ggplot(cor_df, aes(x = pearson)) +
    geom_histogram(bins = 30) +
    theme_classic() +
    ggtitle("Distribution of DEG logFC correlations across targets (Pearson)")

  ggsave(file.path(dir_cor, "DE_correlation_hist_pearson.png"),
         p_hist_cor, width = 8, height = 5, dpi = 200)
}

```

# Evaluate if there is a knockdown when the sgRNA are present.
```{r}
meta <- scobj@meta.data
DefaultAssay(scobj) <- "RNA"
expr <- GetAssayData(scobj, layer = "data")  # v5

nt_label <- "Non_Targeting_Human_CRi"
meta$is_nt <- (!is.na(meta$targeted_gene) & meta$targeted_gene == nt_label) |
              (!is.na(meta$is_non_targeting) & meta$is_non_targeting)

keep <- (!is.na(meta$sgRNA) & !is.na(meta$targeted_gene) & !meta$is_nt) | meta$is_nt
meta <- meta[keep, , drop = FALSE]
expr <- expr[, rownames(meta), drop = FALSE]

plot_dir <- "analyses/2026_PerturbSeq_WT_Kampmann/knockdown_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# genes to test: targeted genes excluding non-targeting label
genes <- sort(unique(na.omit(meta$targeted_gene[!meta$is_nt])))

safe_wilcox <- function(y, grp) {
  out <- list(p.value = NA_real_)
  if (length(y) < 2) return(out)
  if (length(unique(grp)) != 2) return(out)
  # require at least 5 cells per group to avoid nonsense
  tab <- table(grp)
  if (any(tab < 5)) return(out)
  res <- tryCatch(wilcox.test(y ~ grp), error = function(e) NULL)
  if (!is.null(res) && inherits(res, "htest")) out$p.value <- res$p.value
  out
}

results <- list()

for (gene in genes) {
  if (!(gene %in% rownames(expr))) {
    message(sprintf("[Skip] %s not found in assay", gene))
    next
  }

  # sgRNAs for this gene
  sg_list <- sort(unique(na.omit(meta$sgRNA[meta$targeted_gene == gene & !meta$is_nt])))
  if (length(sg_list) == 0) next

  # non-targeting cells
  cells_nt <- rownames(meta)[meta$is_nt]
  if (length(cells_nt) < 5) {
    message("[Skip] too few non-targeting cells overall: ", length(cells_nt))
    next
  }

  for (sg in sg_list) {
    cells_sg <- rownames(meta)[meta$sgRNA == sg & meta$targeted_gene == gene]
    if (length(cells_sg) < 5) next

    cells_use <- c(cells_nt, cells_sg)
    cells_use <- intersect(cells_use, colnames(expr))
    if (length(cells_use) < 10) next

    df <- data.frame(
      expr = as.numeric(expr[gene, cells_use]),
      group = factor(ifelse(cells_use %in% cells_nt, "Non-targeting", sg),
                     levels = c("Non-targeting", sg))
    )

    # drop NA/Inf
    df <- df[is.finite(df$expr), , drop = FALSE]
    if (nrow(df) < 10) next

    wt <- safe_wilcox(df$expr, df$group)

    mean_nt <- mean(df$expr[df$group == "Non-targeting"], na.rm = TRUE)
    mean_sg <- mean(df$expr[df$group == sg],             na.rm = TRUE)

    # On log-normalized scale: difference in means (more stable than ratio)
    delta_log <- mean_sg - mean_nt

    # Optional: approximate fold-change back on linear scale
    # (Seurat "data" is log1p-normalized counts; expm1 reverses log1p)
    mean_nt_lin <- mean(expm1(df$expr[df$group == "Non-targeting"]), na.rm = TRUE)
    mean_sg_lin <- mean(expm1(df$expr[df$group == sg]),             na.rm = TRUE)
    fold_change <- ifelse(mean_nt_lin > 0, mean_sg_lin / mean_nt_lin, NA_real_)
    percent_reduction <- ifelse(!is.na(fold_change), (1 - fold_change) * 100, NA_real_)

    # Plot (only add geom_signif when we actually have a p-value)
    p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.25, size = 0.6, color = "black") +
      labs(
        title = paste0(gene, " expression (log-normalized)"),
        subtitle = paste0("sgRNA: ", sg, " | delta(log)=", signif(delta_log, 3),
                          if (!is.na(wt$p.value)) paste0(" | p=", signif(wt$p.value, 3)) else ""),
        x = "",
        y = paste0(gene, " (RNA@data)")
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")

    if (is.finite(wt$p.value)) {
      p <- p + ggsignif::geom_signif(comparisons = list(c("Non-targeting", sg)),
                                    tip_length = 0.01, map_signif_level = TRUE)
    }

    plot_file <- file.path(plot_dir, paste0(gene, "__", sg, "__knockdown.pdf"))
    ggsave(plot_file, p, height = 5, width = 4)

    results[[paste(gene, sg, sep="__")]] <- data.frame(
      gene               = gene,
      sgRNA              = sg,
      n_non_targeting    = sum(df$group == "Non-targeting"),
      n_sgRNA            = sum(df$group == sg),
      mean_log_non_targeting = mean_nt,
      mean_log_sgRNA         = mean_sg,
      delta_log              = delta_log,
      mean_lin_non_targeting = mean_nt_lin,
      mean_lin_sgRNA         = mean_sg_lin,
      fold_change_lin        = fold_change,
      percent_reduction_lin  = percent_reduction,
      p_value_wilcox         = wt$p.value,
      plot_file              = plot_file,
      stringsAsFactors       = FALSE
    )
  }
}

results_df <- dplyr::bind_rows(results)

# FDR within gene (across sgRNAs) — or set group_by() to nothing for global FDR
results_df <- results_df %>%
  group_by(gene) %>%
  mutate(fdr = p.adjust(p_value_wilcox, method = "BH")) %>%
  ungroup()

out_csv <- "analyses/2026_PerturbSeq_WT_Kampmann/knockdown_results_summary.csv"
write.csv(results_df, out_csv, row.names = FALSE)

# ============================================================
# Knockdown QC summary (big-picture) from sgRNA-level results_df
# ============================================================

out_dir <- "analyses/2026_PerturbSeq_WT_Kampmann/knockdown_QC_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Clean minimal ----
df <- results_df %>%
  mutate(
    gene = as.character(gene),
    sgRNA = as.character(sgRNA),
    delta_log = as.numeric(delta_log),
    fdr = as.numeric(fdr),
    is_sig = is.finite(fdr) & fdr < 0.05 & is.finite(delta_log) & (delta_log < 0)
  ) %>%
  filter(is.finite(delta_log))

# Save cleaned input
fwrite(df, file.path(out_dir, "knockdown_results_cleaned.tsv"), sep = "\t")

# ============================================================
# 1) Primary summary: delta_log per gene (box + sgRNA dots)
# ============================================================
p_delta <- ggplot(df, aes(x = reorder(gene, delta_log, median, na.rm = TRUE), y = delta_log)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85") +
  geom_jitter(aes(color = is_sig), width = 0.2, size = 1, alpha = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
  theme_classic() +
  labs(
    title = "On-target knockdown across sgRNAs",
    subtitle = "Each dot = sgRNA; red = FDR < 0.05 and Δlog < 0",
    x = "Target gene",
    y = "Δ log-expression (sgRNA – non-targeting)",
    color = "Significant\nknockdown"
  )

ggsave(file.path(out_dir, "01_delta_log_per_gene_boxdot.png"), p_delta, width = 10, height = 20)

# ============================================================
# 2) Global distribution of knockdown (all sgRNAs)
# ============================================================
p_global <- ggplot(df, aes(x = delta_log)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    title = "Global distribution of on-target knockdown (all sgRNAs)",
    x = "Δ log-expression (sgRNA – non-targeting)",
    y = "Number of sgRNAs"
  )

ggsave(file.path(out_dir, "02_global_delta_log_hist.png"), p_global, width = 7, height = 5, dpi = 200)

# ============================================================
# 3) Compact volcano-like summary: effect size vs -log10(FDR)
# ============================================================
p_volc_summary <- ggplot(df,
                         aes(x = delta_log,
                             y = -log10(fdr + 1e-300))) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  theme_classic() +
  labs(
    title = "On-target knockdown summary (all sgRNAs)",
    subtitle = "Each dot = sgRNA; horizontal line = FDR 0.05",
    x = "Δ log-expression (sgRNA – non-targeting)",
    y = "-log10(FDR)"
  )

ggsave(file.path(out_dir, "03_knockdown_summary_volcano_like.png"), p_volc_summary, width = 7, height = 6, dpi = 200)

# ============================================================
# 4) A compact “top/bottom” table for quick inspection
# ============================================================
# Per-gene “best guide” and “median guide” summaries
df <- df %>%
  mutate(
    gene = as.character(gene),
    sgRNA = as.character(sgRNA)
  )

df2 <- df
df2$gene <- as.character(as.vector(df2$gene))
df2$sgRNA <- as.character(as.vector(df2$sgRNA))

best_guides <- df2[order(df2$gene, df2$delta_log), ] %>%
  group_by(gene) %>%
  summarise(
    best_sgRNA = sgRNA[1],
    best_delta_log = delta_log[1],
    best_fdr = fdr[1],
    .groups = "drop"
  )

median_guides <- df %>%
  group_by(gene) %>%
  summarise(
    median_delta_log = median(delta_log, na.rm = TRUE),
    frac_sig = mean(is_sig, na.rm = TRUE),
    n_sgRNA = n(),
    .groups = "drop"
  )

summary_table <- hit_summary %>%
  left_join(best_guides, by = "gene") %>%
  left_join(var_summary %>% select(gene, sd_delta, iqr_delta), by = "gene") %>%
  arrange(desc(frac_sig), median_delta_log)

fwrite(summary_table, file.path(out_dir, "knockdown_summary_by_gene.tsv"), sep = "\t")

# Top 20 strongest median knockdowns
top20 <- summary_table %>% arrange(median_delta_log) %>% head(20)
fwrite(top20, file.path(out_dir, "top20_genes_by_median_knockdown.tsv"), sep = "\t")

# Bottom 20 (worst / no knockdown)
bottom20 <- summary_table %>% arrange(desc(median_delta_log)) %>% head(20)
fwrite(bottom20, file.path(out_dir, "bottom20_genes_by_median_knockdown.tsv"), sep = "\t")

```

# Simplified version of the knockdown assessment.
```{r}



# -----------------------------
# Settings
# -----------------------------
set.seed(1)

out_dir <- "knockdown_optional_upgrade"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir_gene  <- file.path(out_dir, "plots_gene_level")
plot_dir_sgRNA <- file.path(out_dir, "plots_sgRNA_level")
dir.create(plot_dir_gene,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir_sgRNA, showWarnings = FALSE, recursive = TRUE)

nt_label <- "Non_Targeting_Human_CRi"
target_col <- "targeted_gene"
sg_col <- "sgRNA"

# Covariates (only used if present)
batch_cols_try <- c("orig.ident", "sample", "bc_wells")  # will choose first found
use_batch <- TRUE

min_cells_per_group_gene <- 50     # targeting vs NT
min_cells_per_group_sg   <- 20     # sgRNA vs NT
expr_detect_threshold    <- 0      # on RNA@data scale; >0 means detected

# -----------------------------
# Pull metadata & expression
# -----------------------------
DefaultAssay(scobj) <- "RNA"

# 1) Make sure the RNA counts layer has dimnames (Seurat v5 Assay5 sometimes stores NULL dimnames)
cnt <- GetAssayData(scobj, layer = "counts")

if (is.null(rownames(cnt))) rownames(cnt) <- rownames(scobj)  # features
if (is.null(colnames(cnt))) colnames(cnt) <- colnames(scobj)  # cells

# 2) Put it back (keeps everything aligned)
scobj[["RNA"]] <- SetAssayData(scobj[["RNA"]], layer = "counts", new.data = cnt)

# 3) Create normalized "data" layer if it doesn't exist yet
# (Your assay only has counts, so GetAssayData(layer="data") would be empty/not present.)
if (!("data" %in% Layers(scobj[["RNA"]]))) {
  scobj <- NormalizeData(scobj, assay = "RNA", normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)
}

# Use log-normalized expression
expr <- GetAssayData(scobj, layer = "data")

# Always rebuild meta fresh (keeps correct rownames)
meta <- scobj@meta.data

target_col <- "targeted_gene"
nt_label <- "Non_Targeting_Human_CRi"

# Robust non-targeting definition (VECTOR-SAFE)
meta$is_nt <- FALSE
if ("is_non_targeting" %in% colnames(meta)) {
  meta$is_nt <- meta$is_nt | (meta$is_non_targeting %in% TRUE)
}
if (target_col %in% colnames(meta)) {
  meta$is_nt <- meta$is_nt | (!is.na(meta[[target_col]]) & meta[[target_col]] == nt_label)
}

# Keep: NT OR (has targeted_gene AND has sgRNA)
keep <- meta$is_nt |
  (!is.na(meta[[target_col]]) & meta[[target_col]] != "" &
     !is.na(meta$sgRNA) & meta$sgRNA != "")

meta <- meta[keep, , drop = FALSE]

# Align safely (now colnames(expr) exists)
cells_common <- intersect(rownames(meta), colnames(expr))
meta <- meta[cells_common, , drop = FALSE]
expr <- expr[, cells_common, drop = FALSE]

message("Kept cells: ", nrow(meta))
message("Non-targeting cells: ", sum(meta$is_nt, na.rm = TRUE))
message("Targeting cells: ", sum(!meta$is_nt, na.rm = TRUE))

genes <- sort(unique(na.omit(as.character(meta[[target_col]])[!meta$is_nt])))
message("Unique targeted genes (non-NT): ", length(genes))


# Helper: safe model fit (glm/lm)
safe_fit <- function(expr_vec, group_vec, ncount_vec = NULL, batch_vec = NULL, family = NULL) {
  # group_vec should be factor with levels c("NT","TARGET")
  df <- data.frame(
    y = expr_vec,
    group = group_vec
  )
  if (!is.null(ncount_vec)) df$log10_nCount <- log10(pmax(ncount_vec, 1))
  if (!is.null(batch_vec))  df$batch <- as.factor(batch_vec)

  # Construct formula
  fml <- if (!is.null(batch_vec) && !is.null(ncount_vec)) {
    y ~ group + log10_nCount + batch
  } else if (!is.null(batch_vec) && is.null(ncount_vec)) {
    y ~ group + batch
  } else if (is.null(batch_vec) && !is.null(ncount_vec)) {
    y ~ group + log10_nCount
  } else {
    y ~ group
  }

  if (is.null(family)) {
    fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
  } else {
    fit <- tryCatch(glm(fml, data = df, family = family), error = function(e) NULL)
  }
  fit
}

# Helper: extract group coefficient p-value and estimate
coef_info <- function(fit) {
  out <- list(beta = NA_real_, p = NA_real_)
  if (is.null(fit)) return(out)
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm)) return(out)
  coefs <- sm$coefficients
  # groupTARGET term name depends on factor coding: "groupTARGET" if levels c("NT","TARGET")
  rn <- rownames(coefs)
  term <- rn[grepl("^group", rn)][1]
  if (is.na(term) || length(term) == 0) return(out)
  out$beta <- coefs[term, 1]
  out$p    <- coefs[term, 4]
  out
}

# ============================================================
# GENE-LEVEL analysis: cells targeting gene (any sgRNA) vs NT
# ============================================================
gene_results <- list()

for (g in genes) {
  if (!(g %in% rownames(expr))) {
    message("[Gene Skip] ", g, " not found in assay rownames")
    next
  }

  cells_t <- rownames(meta)[as.character(meta[[target_col]]) == g]
  cells_nt <- rownames(meta)[meta$is_nt]

  if (length(cells_t) < min_cells_per_group_gene || length(cells_nt) < min_cells_per_group_gene) {
    message("[Gene Skip] ", g, " too few cells. Targeting=", length(cells_t), " NT=", length(cells_nt))
    next
  }

  cells_use <- intersect(c(cells_nt, cells_t), colnames(expr))
  y <- as.numeric(expr[g, cells_use])

  grp <- ifelse(cells_use %in% cells_nt, "NT", "TARGET")
  grp <- factor(grp, levels = c("NT", "TARGET"))

  # Detection (on RNA@data): > 0 indicates detected after log1p norm
  det <- y > expr_detect_threshold

  # Detection rate summary + Fisher exact test
  tab <- table(grp, det)
  fisher_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)

  det_nt <- mean(det[grp == "NT"], na.rm = TRUE)
  det_t  <- mean(det[grp == "TARGET"], na.rm = TRUE)
  det_diff <- det_t - det_nt

  # Wilcoxon on expression
  wt_p <- tryCatch(wilcox.test(y ~ grp)$p.value, error = function(e) NA_real_)

  mean_nt_log <- mean(y[grp == "NT"], na.rm = TRUE)
  mean_t_log  <- mean(y[grp == "TARGET"], na.rm = TRUE)
  delta_log   <- mean_t_log - mean_nt_log

  # Back-transform approximate linear means (expm1)
  mean_nt_lin <- mean(expm1(y[grp == "NT"]), na.rm = TRUE)
  mean_t_lin  <- mean(expm1(y[grp == "TARGET"]), na.rm = TRUE)
  fold_change <- ifelse(mean_nt_lin > 0, mean_t_lin / mean_nt_lin, NA_real_)
  pct_reduction <- ifelse(!is.na(fold_change), (1 - fold_change) * 100, NA_real_)

  # Logistic regression: detection ~ group + log10(nCount_RNA) (+ batch)
  ncount_vec <- if ("nCount_RNA" %in% colnames(meta)) meta[cells_use, "nCount_RNA"] else NULL
  batch_vec  <- if (!is.null(batch_col)) meta[cells_use, batch_col] else NULL

  fit_logit <- safe_fit(expr_vec = as.numeric(det), group_vec = grp,
                        ncount_vec = ncount_vec, batch_vec = batch_vec,
                        family = binomial())

  logit_info <- coef_info(fit_logit)

  # Linear regression: expression ~ group + log10(nCount_RNA) (+ batch)
  fit_lm <- safe_fit(expr_vec = y, group_vec = grp,
                     ncount_vec = ncount_vec, batch_vec = batch_vec,
                     family = NULL)

  lm_info <- coef_info(fit_lm)

  # Plot: expression
  dfp <- data.frame(expr = y, group = grp)
  p_expr <- ggplot(dfp, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.25, size = 0.4, color = "black") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(g, " (gene-level KD)"),
      subtitle = paste0("Δlog=", signif(delta_log,3),
                        " | %red≈", signif(pct_reduction,3),
                        " | wilcox p=", signif(wt_p,3),
                        " | detNT=", signif(det_nt,3),
                        " detT=", signif(det_t,3),
                        " | fisher p=", signif(fisher_p,3)),
      x = "", y = paste0(g, " (RNA@data)")
    ) +
    theme(legend.position = "none") +
    ggsignif::geom_signif(comparisons = list(c("NT", "TARGET")), test = "wilcox.test", tip_length = 0.01)

  ggsave(file.path(plot_dir_gene, paste0(g, "__expr_genelevel.pdf")),
         p_expr, width = 4.3, height = 5, dpi = 200)

  # Plot: detection rate bar
  det_df <- data.frame(
    group = c("NT", "TARGET"),
    detect_rate = c(det_nt, det_t)
  )
  p_det <- ggplot(det_df, aes(x = group, y = detect_rate, fill = group)) +
    geom_col(alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(g, " detection rate"),
      subtitle = paste0("Δdetect=", signif(det_diff,3),
                        " | fisher p=", signif(fisher_p,3),
                        " | logit beta=", signif(logit_info$beta,3),
                        " p=", signif(logit_info$p,3)),
      x = "", y = paste0("Pr(", g, " detected)")
    ) +
    theme(legend.position = "none")

  ggsave(file.path(plot_dir_gene, paste0(g, "__detect_genelevel.pdf")),
         p_det, width = 4.3, height = 4.3, dpi = 200)

  gene_results[[g]] <- data.frame(
    gene = g,
    n_targeting = length(cells_t),
    n_nt = length(cells_nt),

    mean_log_nt = mean_nt_log,
    mean_log_targeting = mean_t_log,
    delta_log = delta_log,

    mean_lin_nt = mean_nt_lin,
    mean_lin_targeting = mean_t_lin,
    fold_change_lin = fold_change,
    percent_reduction_lin = pct_reduction,

    detect_rate_nt = det_nt,
    detect_rate_targeting = det_t,
    detect_rate_diff = det_diff,

    p_wilcox_expr = wt_p,
    p_fisher_detect = fisher_p,

    logit_beta_group = logit_info$beta,
    logit_p_group = logit_info$p,

    lm_beta_group = lm_info$beta,
    lm_p_group = lm_info$p,

    stringsAsFactors = FALSE
  )
}

gene_df <- bind_rows(gene_results)

# Adjust p-values across genes (optional)
if (nrow(gene_df) > 0) {
  gene_df$fdr_wilcox_expr   <- p.adjust(gene_df$p_wilcox_expr, method = "BH")
  gene_df$fdr_fisher_detect <- p.adjust(gene_df$p_fisher_detect, method = "BH")
  gene_df$fdr_logit_group   <- p.adjust(gene_df$logit_p_group, method = "BH")
  gene_df$fdr_lm_group      <- p.adjust(gene_df$lm_p_group, method = "BH")

  fwrite(as.data.table(gene_df), file.path(out_dir, "gene_level_knockdown_stats.tsv"), sep = "\t")
}

# ============================================================
# sgRNA-LEVEL analysis: each sgRNA vs NT (on-target gene expression)
# ============================================================
sg_results <- list()

# Precompute NT cells once
cells_nt_all <- rownames(meta)[meta$is_nt]
if (length(cells_nt_all) < min_cells_per_group_sg) {
  warning("Too few NT cells for sgRNA-level analysis: ", length(cells_nt_all))
}

for (g in genes) {
  if (!(g %in% rownames(expr))) next

  # sgRNAs targeting g
  sg_list <- sort(unique(na.omit(as.character(meta[[sg_col]][as.character(meta[[target_col]]) == g & !meta$is_nt]))))
  if (length(sg_list) == 0) next

  for (sg in sg_list) {
    cells_sg <- rownames(meta)[as.character(meta[[sg_col]]) == sg & as.character(meta[[target_col]]) == g]
    if (length(cells_sg) < min_cells_per_group_sg || length(cells_nt_all) < min_cells_per_group_sg) next

    cells_use <- intersect(c(cells_nt_all, cells_sg), colnames(expr))
    y <- as.numeric(expr[g, cells_use])

    grp <- ifelse(cells_use %in% cells_nt_all, "NT", "SG")
    grp <- factor(grp, levels = c("NT", "SG"))

    det <- y > expr_detect_threshold
    tab <- table(grp, det)
    fisher_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)

    det_nt <- mean(det[grp == "NT"], na.rm = TRUE)
    det_sg <- mean(det[grp == "SG"], na.rm = TRUE)

    wt_p <- tryCatch(wilcox.test(y ~ grp)$p.value, error = function(e) NA_real_)

    mean_nt_log <- mean(y[grp == "NT"], na.rm = TRUE)
    mean_sg_log <- mean(y[grp == "SG"], na.rm = TRUE)
    delta_log   <- mean_sg_log - mean_nt_log

    mean_nt_lin <- mean(expm1(y[grp == "NT"]), na.rm = TRUE)
    mean_sg_lin <- mean(expm1(y[grp == "SG"]), na.rm = TRUE)
    fold_change <- ifelse(mean_nt_lin > 0, mean_sg_lin / mean_nt_lin, NA_real_)
    pct_red <- ifelse(!is.na(fold_change), (1 - fold_change) * 100, NA_real_)

    ncount_vec <- if ("nCount_RNA" %in% colnames(meta)) meta[cells_use, "nCount_RNA"] else NULL
    batch_vec  <- if (!is.null(batch_col)) meta[cells_use, batch_col] else NULL

    fit_logit <- safe_fit(expr_vec = as.numeric(det), group_vec = grp,
                          ncount_vec = ncount_vec, batch_vec = batch_vec,
                          family = binomial())
    logit_info <- coef_info(fit_logit)

    fit_lm <- safe_fit(expr_vec = y, group_vec = grp,
                       ncount_vec = ncount_vec, batch_vec = batch_vec,
                       family = NULL)
    lm_info <- coef_info(fit_lm)

    # Plot expression (sgRNA vs NT)
    dfp <- data.frame(expr = y, group = grp)
    p_expr <- ggplot(dfp, aes(x = group, y = expr, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.25, size = 0.4, color = "black") +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(g, " KD: ", sg),
        subtitle = paste0("Δlog=", signif(delta_log,3),
                          " | %red≈", signif(pct_red,3),
                          " | wilcox p=", signif(wt_p,3)),
        x = "", y = paste0(g, " (RNA@data)")
      ) +
      theme(legend.position = "none")

    if (is.finite(wt_p)) {
      p_expr <- p_expr +
        ggsignif::geom_signif(comparisons = list(c("NT","SG")), test = "wilcox.test", tip_length = 0.01)
    }

    ggsave(file.path(plot_dir_sgRNA, paste0(g, "__", sg, "__expr_sgRNAlevel.pdf")),
           p_expr, width = 4.3, height = 5, dpi = 200)

    sg_results[[paste(g, sg, sep="__")]] <- data.frame(
      gene = g,
      sgRNA = sg,
      n_sg = length(cells_sg),
      n_nt = length(cells_nt_all),

      delta_log = delta_log,
      fold_change_lin = fold_change,
      percent_reduction_lin = pct_red,

      detect_rate_nt = det_nt,
      detect_rate_sg = det_sg,

      p_wilcox_expr = wt_p,
      p_fisher_detect = fisher_p,

      logit_beta_group = logit_info$beta,
      logit_p_group = logit_info$p,

      lm_beta_group = lm_info$beta,
      lm_p_group = lm_info$p,

      stringsAsFactors = FALSE
    )
  }
}

sg_df <- bind_rows(sg_results)

# Adjust p-values across sgRNAs globally (optional)
if (nrow(sg_df) > 0) {
  sg_df$fdr_wilcox_expr   <- p.adjust(sg_df$p_wilcox_expr, method = "BH")
  sg_df$fdr_fisher_detect <- p.adjust(sg_df$p_fisher_detect, method = "BH")
  sg_df$fdr_logit_group   <- p.adjust(sg_df$logit_p_group, method = "BH")
  sg_df$fdr_lm_group      <- p.adjust(sg_df$lm_p_group, method = "BH")

  fwrite(as.data.table(sg_df), file.path(out_dir, "sgRNA_level_knockdown_stats.tsv"), sep = "\t")
}

# ============================================================
# GLOBAL KNOCKDOWN SUMMARY FIGURE
# Uses sg_df (sgRNA-level results)
# ============================================================

# Clean
df_global <- sg_df %>%
  filter(is.finite(delta_log)) %>%
  mutate(
    sig = ifelse(!is.na(fdr_wilcox_expr) & fdr_wilcox_expr < 0.05 & delta_log < 0,
                 "Significant KD", "Not significant")
  )

# ------------------------------------------------------------
# Panel 1: Histogram + density
# ------------------------------------------------------------
p1 <- ggplot(df_global, aes(x = delta_log)) +
  geom_histogram(aes(y = ..density..),
                 bins = 60,
                 fill = "grey80",
                 color = "grey40") +
  geom_density(color = "black", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    title = "Global on-target knockdown (sgRNA-level)",
    x = expression(Delta~"log-expression (sgRNA - NT)"),
    y = "Density"
  )

# ------------------------------------------------------------
# Panel 2: Distribution split by significance
# ------------------------------------------------------------
p2 <- ggplot(df_global, aes(x = delta_log, fill = sig)) +
  geom_histogram(position = "identity",
                 alpha = 0.6,
                 bins = 60) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Significant KD" = "#D55E00",
                               "Not significant" = "grey70")) +
  theme_classic() +
  labs(
    title = "Knockdown by statistical support",
    x = expression(Delta~"log-expression"),
    y = "Count",
    fill = ""
  )

# ------------------------------------------------------------
# Panel 3: Cumulative distribution (very interpretable)
# ------------------------------------------------------------
p3 <- ggplot(df_global, aes(x = delta_log, color = sig)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Significant KD" = "#D55E00",
                                "Not significant" = "grey40")) +
  theme_classic() +
  labs(
    title = "Cumulative distribution of knockdown",
    x = expression(Delta~"log-expression"),
    y = "Cumulative fraction",
    color = ""
  )

# ------------------------------------------------------------
# Panel 4: Percent reduction (linear scale)
# ------------------------------------------------------------
if ("percent_reduction_lin" %in% colnames(df_global)) {
  p4 <- ggplot(df_global, aes(x = percent_reduction_lin)) +
    geom_histogram(bins = 60, fill = "#56B4E9", color = "black") +
    theme_classic() +
    labs(
      title = "Approximate percent reduction (linear scale)",
      x = "% reduction",
      y = "Count"
    )
} else {
  p4 <- ggplot() + theme_void() + labs(title = "Percent reduction not available")
}

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------
p_combined <- (p1 | p2) / (p3 | p4)

ggsave(
  file.path(out_dir, "GLOBAL_KNOCKDOWN_SUMMARY.png"),
  p_combined,
  width = 12,
  height = 10,
  dpi = 300
)


```

# Analysis groupping all sgRNAs targetting the same gene.
```{r}
meta <- scobj@meta.data
expr <- GetAssayData(scobj, slot = "data")   # normalized/log data

# Keep only cells that have either a targeting sgRNA or are non-targeting
keep <- (!is.na(meta$sgRNA) & !is.na(meta$targeted_gene)) | (meta$is_non_targeting %in% TRUE)
meta <- meta[keep, , drop = FALSE]
expr <- expr[, keep, drop = FALSE]

# Directory for plots
plot_dir <- "knockdown_plots_grouped"
dir.create(plot_dir, showWarnings = FALSE)

# Enumerate genes
genes <- sort(unique(na.omit(meta$targeted_gene[!meta$is_non_targeting])))
results_grouped <- vector("list", length = 0)

# Helper: safe Wilcoxon
safe_wilcox <- function(y, grp) {
  out <- list(p.value = NA_real_, method = "wilcox")
  if (length(y) < 2) return(out)
  if (length(unique(grp)) != 2) return(out)
  res <- tryCatch(wilcox.test(y ~ grp), error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(res) && inherits(res, "htest") && !is.null(res$p.value)) {
    out$p.value <- res$p.value
  }
  out
}

# Iterate: for each targeted gene, group all sgRNAs targeting it
for (gene in genes) {
  if (!(gene %in% rownames(expr))) {
    message(sprintf("[Skip] %s not found in expression matrix", gene))
    next
  }
  gene_expr <- as.numeric(expr[gene, ])

  # Cells with any sgRNA targeting this gene
  cells_targeting <- rownames(meta)[meta$targeted_gene == gene & !meta$is_non_targeting]
  # Cells with non-targeting sgRNA
  cells_nt <- rownames(meta)[meta$is_non_targeting %in% TRUE]

  # Build data frame for plotting/testing
  idx <- c(match(cells_targeting, colnames(expr)), match(cells_nt, colnames(expr)))
  df <- data.frame(
    expr  = gene_expr[idx],
    group = factor(c(rep("Targeting", length(cells_targeting)), rep("Non-targeting", length(cells_nt))),
                   levels = c("Non-targeting", "Targeting"))
  )

  # Plot
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.25, size = 0.6, color = "black") +
    labs(
      title = paste(gene, "expression by targeting sgRNA presence"),
      x = "",
      y = paste(gene, "normalized expression")
    ) +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
    geom_signif(comparisons = list(c("Non-targeting", "Targeting")), tip_length = 0.01) +
    theme(legend.position = "none")

  plot_file <- file.path(plot_dir, paste0(gene, "_grouped_knockdown.pdf"))
  ggsave(plot_file, p, height = 5, width = 4)

  # Statistics
  wt <- safe_wilcox(df$expr, df$group)

  mean_nt <- mean(df$expr[df$group == "Non-targeting"], na.rm = TRUE)
  mean_targeting <- mean(df$expr[df$group == "Targeting"], na.rm = TRUE)
  fold_change <- mean_targeting / mean_nt
  percent_reduction <- (1 - fold_change) * 100

  results_grouped[[gene]] <- data.frame(
    gene               = gene,
    n_non_targeting    = sum(df$group == "Non-targeting"),
    n_targeting        = sum(df$group == "Targeting"),
    mean_non_targeting = mean_nt,
    mean_targeting     = mean_targeting,
    fold_change        = fold_change,
    percent_reduction  = percent_reduction,
    p_value_wilcox     = wt$p.value,
    test_method        = wt$method,
    plot_file          = plot_file,
    stringsAsFactors   = FALSE
  )
}

# Collate results
results_df_grouped <- dplyr::bind_rows(results_grouped)
results_df_grouped <- results_df_grouped %>%
    mutate(fdr = p.adjust(p_value_wilcox, method = "BH"))

write.csv(results_df_grouped, "data/Seqmatic_112025_Parse_100K_Kampmann/knockdown_results_grouped_summary.csv", row.names = FALSE)

```
