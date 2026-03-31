# Obtaining transcriptomic signals for each sgRNA and target

## ------------------------------------------------------------
## STEP 0: Load libraries
## ------------------------------------------------------------

```{r}
library(Seurat) # make sure v5 is loaded
library(BPCells)
library(dplyr)
library(RANN)
library(readr)
library(tibble)
library(stringr)
library(edgeR)
library(Matrix)
library(tidyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(forcats)
library(patchwork)
library(scales)
library(ggrepel)
```

## ------------------------------------------------------------
## STEP 1: Preprocessing
## ------------------------------------------------------------

```{r}
md <- scobj_filt@meta.data
md$orig.ident <- as.factor(md$orig.ident)
md$gRNA   <- as.character(md$gRNA)
md$target <- as.character(md$target)
DefaultAssay(scobj_filt) <- "RNA"
```

## ------------------------------------------------------------
## STEP 2: Run the Mixscape pipeline
## See: https://satijalab.org/seurat/articles/mixscape_vignette.html
## ------------------------------------------------------------

```{r}
## ------------------------------------------------------------
## Normalize and scale the counts
## ------------------------------------------------------------

scobj_filt$mix_label <- ifelse(scobj_filt$target == "NTC", "NTC", scobj_filt$target)
scobj_filt <- NormalizeData(scobj_filt, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj_filt <- FindVariableFeatures(scobj_filt, selection.method = "vst", nfeatures = 3000)
scobj_filt <- ScaleData(scobj_filt)

## ------------------------------------------------------------
## Perform PC reduction
## ------------------------------------------------------------

scobj_filt <- RunPCA(scobj_filt, npcs = 50)

# Note: sometimes the regular RunPCA() fails due to memory-restrictions so can perform it using BPCell functions. Leaving code below just in case.
# See: https://bnprks.github.io/BPCells/articles/pbmc3k.html#rna-normalization-pca-and-umap
#counts_norm <- GetAssayData(scobj_filt, layer = "scale.data")
#svd <- BPCells::svds(mat_norm, k=50)
#pc_embed   <- multiply_cols(svd$v, svd$d) # cells × dimensions
#pc_loadings <- svd$u # genes × dimensions
#cells_bpc <- colnames(counts_norm)
#genes_bpc <- rownames(counts_norm)        
#pc_embed <- pc_embed[match(colnames(scobj_filt), cells_bpc), , drop = FALSE]
#rownames(pc_embed) <- colnames(scobj_filt)
#hvg <- VariableFeatures(scobj_filt)
#hvg_in <- intersect(hvg, genes_bpc)
#pc_loadings <- pc_loadings[match(hvg_in, genes_bpc), , drop = FALSE]
#rownames(pc_loadings) <- hvg_in                               
#pc_names <- paste0("PC_", seq_len(ncol(pc_embed)))
#colnames(pc_embed)    <- pc_names
#colnames(pc_loadings) <- pc_names
#k <- ncol(pc_embed)
#scobj_filt[["pca"]] <- CreateDimReducObject(embeddings = pc_embed, loadings= pc_loadings, stdev= svd$d[seq_len(k)], key= "PC_", assay= "RNA")

## ------------------------------------------------------------
## Calculate 'Perturbation Score'
## ------------------------------------------------------------

# The CalcPerturbSig() function failed due to memory-restrictions so will run it manually.
# scobj_filt <- CalcPerturbSig(object = scobj_filt, assay = "RNA", gd.class = "mix_label", nt.cell.class = "NTC", new.assay.name = "PRTB", reduction = "pca", num.neighbors = 20, ndims= 50)

DefaultAssay(scobj_filt) <- "RNA"

# Use HVGs to keep memory bounded (you can change nfeatures if needed)
hvg <- VariableFeatures(scobj_filt)

# Pull sparse log-normalized data (genes x cells), then subset to HVGs
E_all  <- GetAssayData(scobj_filt, layer = "data")       # dgCMatrix
E_all  <- E_all[hvg, , drop = FALSE]

# Nearest NTC neighbors from your PCA (cells x PCs)
PC_all  <- Embeddings(scobj_filt, "pca")
cells    <- colnames(scobj_filt)
ntc_mask <- scobj_filt$mix_label == "NTC"
PC_ntc <- PC_all[ntc_mask, , drop = FALSE]
PC_qry <- PC_all
k      <- 20
nn <- RANN::nn2(data = PC_ntc, query = PC_qry, k = k, searchtype = "standard")
nn_idx_rel  <- nn$nn.idx                 # indices into NTC subset
ntc_cells   <- which(ntc_mask)           # positions of NTC in all cells
nn_idx_glob <- matrix(ntc_cells[nn_idx_rel], nrow = nrow(nn_idx_rel), ncol = ncol(nn_idx_rel))

# Expression of NTC cells (HVGs x NTC) once
E_ntc <- E_all[, ntc_cells, drop = FALSE]
n_cells <- ncol(E_all)
n_genes <- nrow(E_all)

# Choose a reasonably small chunk size
chunk_n <- 2000L
n_blocks <- ceiling(n_cells / chunk_n)
res_blocks <- vector("list", n_blocks)

# Process cells in blocks
b <- 0L
for (start in seq(1L, n_cells, by = chunk_n)) {
  idx <- start:min(start + chunk_n - 1L, n_cells)
  b <- b + 1L

  # Neighbor indices for this block (relative to NTC set)
  nn_chunk_glob <- nn_idx_glob[idx, , drop = FALSE]
  nn_chunk_rel  <- match(nn_chunk_glob, ntc_cells)  # map to 1..#NTC
  if (anyNA(nn_chunk_rel)) stop("NA in neighbor mapping; check mix_label/PC indices.")

  # Initialize a sparse sum matrix (HVGs x block_size)
  sum_ntc <- Matrix(0, n_genes, length(idx), sparse = TRUE)
  # Accumulate over k neighbors: each step pulls an HVGs x block_size matrix
  for (m in seq_len(k)) {
    # columns-of-interest in E_ntc for the m-th neighbor of each cell in the block
    cols_m <- nn_chunk_rel[, m]
    # Subset returns HVGs x block_size; accumulate
    sum_ntc <- sum_ntc + E_ntc[, cols_m, drop = FALSE]
  }
  mean_ntc <- sum_ntc / k

  # Current block of the full matrix
  curr <- E_all[, idx, drop = FALSE]

  # Residuals block (still sparse); drop structural zeros
  res_block <- Matrix::drop0(curr - mean_ntc)

  # Collect the block; DO NOT assign into a giant preallocated object
  res_blocks[[b]] <- res_block

  # free interm to keep peak RSS down
  rm(sum_ntc, mean_ntc, curr, res_block); gc()
}

# Single concatenation to create the full residual matrix (HVGs x cells)
residuals <- do.call(Matrix::cBind, res_blocks)
dimnames(residuals) <- list(rownames(E_all), colnames(E_all))

# Store residuals into a new assay 'PRTB' (layer 'data')
prt <- CreateAssayObject(counts = Matrix(0, nrow(residuals), ncol(residuals), sparse = TRUE, dimnames = dimnames(residuals)))
prt$data <- residuals
scobj_filt[["PRTB"]] <- prt

# Save object
saveRDS(scobj_filt, file="output_combined_filtered/scobj_filt.rds")

## ------------------------------------------------------------
## Run Mixscape to classify cells into NP (non-perturbed), KO (Perturbed), NT (controls)
## ------------------------------------------------------------

scobj_filt <- RunMixscape(object = scobj_filt, assay = "PRTB", labels = "mix_label", nt.class.name = "NTC", 
    new.class.name = "mixscape_class", fine.mode = TRUE, fine.mode.labels = "gRNA", 
    verbose = TRUE)

## ------------------------------------------------------------
## Plot summaries of perturbation classifications
## ------------------------------------------------------------

obj <- scobj_filt
col_global   <- "mixscape_class.global"
col_fine     <- "mixscape_class" 
col_pko      <- "mixscape_class_p_ko"
col_target   <- "target"
col_label    <- "mix_label"
col_gRNA     <- "gRNA"
ntc_name <- "NTC"

# Plot sizing and selection
top_n_targets_for_facets <- 30
min_cells_per_target     <- 200 
outdir <- "output_combined_filtered/mixscape_plots"
dir.create(outdir, showWarnings = FALSE)

# Grab metadata and make clean target label

md <- obj@meta.data %>%
  mutate(
    target_id = .data[[col_target]],
    global = .data[[col_global]],
    p_ko = suppressWarnings(as.numeric(.data[[col_pko]])),
    gRNA = .data[[col_gRNA]]
  )

# Plot #1: Overall composition

p_overall <- md %>%
  count(global) %>%
  mutate(frac = n / sum(n)) %>%
  ggplot(aes(x = fct_relevel(global, ntc_name, "NP", "KO"), y = frac, fill = global)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Mixscape global class composition",
    x = NULL, y = "Fraction of cells"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "01_overall_global_composition.png"),
       p_overall, width = 6, height = 4, dpi = 200)

# Plot #2: Per-target KO rate

per_target <- md %>%
  filter(!is.na(target_id)) %>%
  group_by(target_id) %>%
  summarise(
    n_cells = n(),
    frac_KO = mean(global == "KO"),
    frac_NP = mean(global == "NP"),
    frac_NTC = mean(global == ntc_name),
    med_pko = median(p_ko, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_target) %>%
  arrange(desc(frac_KO))

p_rank <- per_target %>%
  mutate(target_id = fct_inorder(target_id)) %>%
  ggplot(aes(x = target_id, y = frac_KO, size = n_cells, color = n_cells)) +
  geom_point(alpha = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = sprintf("KO fraction per target (n ≥ %d cells/target)", min_cells_per_target),
    x = NULL, y = "KO fraction", size = "# cells", color = "# cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outdir, "02_ranked_KO_fraction_per_target.png"),
       p_rank, width = 10, height = 4.5, dpi = 200)

# Plot #3: Show the top 50 targets with labels

top_label_n <- min(50, nrow(per_target))
p_rank_labeled <- per_target %>%
  slice_head(n = top_label_n) %>%
  mutate(target_id = fct_reorder(target_id, frac_KO)) %>%
  ggplot(aes(x = target_id, y = frac_KO)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = sprintf("Top %d targets by KO fraction", top_label_n),
    x = NULL, y = "KO fraction"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "03_top_targets_by_KO_fraction.png"),
       p_rank_labeled, width = 7, height = 10, dpi = 200)

# Plot #4: Top 50 targets, faceted

top_targets <- per_target %>%
  filter(target_id != ntc_name) %>%
  slice_head(n = top_n_targets_for_facets) %>%
  pull(target_id)

df_long <- md %>%
  filter(target_id %in% top_targets) %>%
  group_by(target_id, gRNA, global) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    global = factor(global, levels = c(ntc_name, "NP", "KO"))
  )

p_facets <- ggplot(df_long, aes(x = gRNA, y = frac * 100, fill = global)) +
  geom_col(width = 0.9) +
  facet_wrap(~ target_id, scales = "free_x", ncol = 5) +
  theme_classic(base_size = 11) +
  labs(
    title = sprintf("Mixscape class per gRNA (top %d targets by KO rate)", length(top_targets)),
    x = "gRNA", y = "% of cells", fill = "Mixscape class"
  ) +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "04_faceted_class_per_gRNA_top_targets.png"),
       p_facets, width = 14, height = 9, dpi = 200)

# Plot #5: KO confidence

# Compare p_ko distributions by global class, and per-target medians

p_pko_by_class <- md %>%
  filter(!is.na(p_ko)) %>%
  ggplot(aes(x = global, y = p_ko)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  theme_classic(base_size = 12) +
  labs(
    title = "KO confidence distribution by global class",
    x = NULL, y = "mixscape_class_p_ko"
  )

ggsave(file.path(outdir, "05_pko_distribution_by_global_class.png"),
       p_pko_by_class, width = 7, height = 4.5, dpi = 200)

# Plot #6: KO confidence per target


p_target_pko <- per_target %>%
  filter(target_id != ntc_name) %>%
  slice_head(n = min(60, nrow(per_target))) %>%
  mutate(target_id = fct_reorder(target_id, med_pko)) %>%
  ggplot(aes(x = target_id, y = med_pko)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Median KO confidence per target (top targets by KO rate shown)",
    x = NULL, y = "Median mixscape_class_p_ko"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "06_median_pko_per_target.png"),
       p_target_pko, width = 7, height = 9, dpi = 200)

```

## ------------------------------------------------------------
## STEP 3: Pseudobulk aggregation
##  A) by gRNA × orig.ident (per-guide vs NTC)
##  B) by target × orig.ident (per-target vs NTC)
## ------------------------------------------------------------

### Define helper functions and set up

```{r}
obj <- scobj_filt
assay_use <- "RNA"
layer_counts <- "counts"           # counts layer to use
outdir <- "pseudoRep_edgeR_results"
dir.create(outdir, showWarnings = FALSE)

# Pseudo-replicate settings
K <- 4                             # number of pseudo-replicates
min_cells_per_rep <- 50            # minimum cells per pseudo-rep (tunable)
max_groups_to_run <- Inf           # set to limit number of guides/targets processed (for testing)

# QC / plotting
plots_dir <- file.path(outdir, "plots")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(file.path(outdir, "concordance_plots"), showWarnings = FALSE)

# -----------------------------
# Helpers: groups and chunked pseudobulk
# -----------------------------
make_groups <- function(f) {
  f <- droplevels(factor(f))
  split(seq_along(f), f)
}

# chunked column-sum for a set of cell indices
sum_cols_chunked <- function(counts_obj, idx, chunk_size = 5000) {
  # returns numeric vector length = nrow(counts_obj)
  gnames <- rownames(counts_obj)
  if (length(idx) == 0) return(numeric(length(gnames)))
  if (length(idx) == 1) {
    tmp <- as(counts_obj[, idx, drop = FALSE], "dgCMatrix")
    return(as.numeric(Matrix::rowSums(tmp)))
  }
  # else chunk
  chunks <- split(idx, ceiling(seq_along(idx) / chunk_size))
  svec <- numeric(length(gnames))
  for (ch in chunks) {
    tmp <- as(counts_obj[, ch, drop = FALSE], "dgCMatrix")
    svec <- svec + as.numeric(Matrix::rowSums(tmp))
  }
  svec
}

# Build genes x K pseudobulk for a given set of cell indices and fold labels
build_pseudorep_matrix <- function(counts_obj, cell_idx, folds, K, chunk_size = 5000) {
  gnames <- rownames(counts_obj)
  out_cols <- vector("list", K)
  for (k in seq_len(K)) {
    ids <- cell_idx[which(folds == k)]
    if (length(ids) == 0) return(NULL)   # missing fold
    svec <- sum_cols_chunked(counts_obj, ids, chunk_size = chunk_size)
    out_cols[[k]] <- Matrix(svec, ncol = 1, sparse = TRUE)
  }
  pb <- do.call(cbind, out_cols)
  rownames(pb) <- gnames
  colnames(pb) <- paste0("rep", seq_len(K))
  pb
}

make_folds <- function(n, K = 4, seed = 1) {
  set.seed(seed)
  if (n < K) return(rep(1, n))
  sample(rep(seq_len(K), length.out = n))
}

# -----------------------------
# Prep: counts (BPCells RenameDims) and metadata
# -----------------------------
DefaultAssay(obj) <- assay_use
counts <- GetAssayData(obj, assay = assay_use, layer = layer_counts)  # RenameDims / BPCells
md <- obj@meta.data
stopifnot(ncol(counts) == nrow(md))

# Collapse NTC guides into "NTC"
md$gRNA_group <- ifelse(md$target == "NTC" | grepl("^NTC", md$gRNA), "NTC", md$gRNA)
obj@meta.data$gRNA_group <- md$gRNA_group

# -----------------------------
# Build cell-index lists
# -----------------------------
cells_all <- seq_len(nrow(md))
cells_by_gRNA <- split(cells_all, md$gRNA_group)
cells_by_target <- split(cells_all, md$target)

# -----------------------------
# Core: function to test one group vs NTC using pseudo-reps
# -----------------------------
run_pseudorep_edger_one <- function(group_name, group_idx, ntc_idx,
                                   counts_obj, K = 4, min_cells_per_rep = 50,
                                   seed = 1, chunk_size = 5000) {
  # require enough cells overall
  if (length(group_idx) < K * min_cells_per_rep) return(NULL)
  if (length(ntc_idx) < K * min_cells_per_rep) return(NULL)

  folds_grp <- make_folds(length(group_idx), K = K, seed = seed)
  folds_ntc <- make_folds(length(ntc_idx),   K = K, seed = seed + 999)

  # ensure each fold has enough cells
  if (any(tabulate(folds_grp, nbins = K) < min_cells_per_rep)) return(NULL)
  if (any(tabulate(folds_ntc, nbins = K) < min_cells_per_rep)) return(NULL)

  pb_grp <- build_pseudorep_matrix(counts_obj, group_idx, folds_grp, K, chunk_size = chunk_size)
  pb_ntc <- build_pseudorep_matrix(counts_obj, ntc_idx, folds_ntc, K, chunk_size = chunk_size)
  if (is.null(pb_grp) || is.null(pb_ntc)) return(NULL)

  pb <- cbind(pb_ntc, pb_grp)
  colnames(pb) <- c(paste0("NTC.rep", seq_len(K)), paste0(group_name, ".rep", seq_len(K)))

  group_fac <- factor(c(rep("NTC", K), rep(group_name, K)), levels = c("NTC", group_name))

  y <- DGEList(counts = pb)
  keep <- filterByExpr(y, group = group_fac)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) < 200) return(NULL)

  y <- calcNormFactors(y)
  design <- model.matrix(~ group_fac)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  coef_name <- paste0("group_fac", make.names(group_name))
  if (!(coef_name %in% colnames(design))) return(NULL)

  qlf <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)
  tt$tested_group <- group_name
  tt$K <- K
  tt$n_cells_group <- length(group_idx)
  tt$n_cells_ntc <- length(ntc_idx)
  tt
}

# -----------------------------
# Wrapper: run many groups (guides or targets)
# -----------------------------
run_many_groups <- function(group_to_idx, ntc_idx,
                            counts_obj, K = 4, min_cells_per_rep = 50,
                            seed = 1, progress_every = 50, limit = Inf) {
  groups <- names(group_to_idx)
  res_list <- list()
  tot <- min(length(groups), limit)
  for (i in seq_len(tot)) {
    g <- groups[i]
    if (g == "NTC") next
    res <- tryCatch({
      run_pseudorep_edger_one(
        group_name = g,
        group_idx = group_to_idx[[g]],
        ntc_idx = ntc_idx,
        counts_obj = counts_obj,
        K = K,
        min_cells_per_rep = min_cells_per_rep,
        seed = seed + i
      )
    }, error = function(e) {
      message("Error for group ", g, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(res)) res_list[[g]] <- res
    if (i %% progress_every == 0) message("Processed ", i, " / ", tot, " groups")
  }
  if (!length(res_list)) return(NULL)
  bind_rows(res_list)
}
```

### Run DE test per sgRNA and per targeted region

```{r}
# -----------------------------
# Run per-guide (gRNA_group) vs NTC
# -----------------------------
ntc_idx_guides <- cells_by_gRNA[["NTC"]]
guide_res <- run_many_groups(
  group_to_idx = cells_by_gRNA,
  ntc_idx = ntc_idx_guides,
  counts_obj = counts,
  K = K,
  min_cells_per_rep = min_cells_per_rep,
  seed = 1,
  progress_every = 100,
  limit = max_groups_to_run
)
if (!is.null(guide_res)) {
  write.csv(guide_res, file.path(outdir, "DE_perGuide_vsNTC_edgeR_pseudoreps.csv"), row.names = FALSE)
  message("Wrote DE_perGuide_vsNTC_edgeR_pseudoreps.csv")
} else {
  message("No guide-level contrasts produced (too few cells / folds).")
}

# -----------------------------
# Run per-target vs NTC
# -----------------------------
ntc_idx_targets <- cells_by_gRNA[["NTC"]]  # same NTC pool
target_res <- run_many_groups(
  group_to_idx = cells_by_target,
  ntc_idx = ntc_idx_targets,
  counts_obj = counts,
  K = K,
  min_cells_per_rep = min_cells_per_rep,
  seed = 20000,
  progress_every = 50,
  limit = max_groups_to_run
)
if (!is.null(target_res)) {
  write.csv(target_res, file.path(outdir, "DE_perTarget_vsNTC_edgeR_pseudoreps.csv"), row.names = FALSE)
  message("Wrote DE_perTarget_vsNTC_edgeR_pseudoreps.csv")
} else {
  message("No target-level contrasts produced (too few cells / folds).")
}
```

### Perform guide concordance tests

```{r}
outdir <- "sgRNA_concordance"
dir.create(outdir, showWarnings = FALSE)

# significance definition
fdr_thresh <- 0.05
lfc_thresh <- 0.0 

# correlation settings
cor_method <- "pearson"
min_genes_all <- 2000     # min overlapping genes for "all genes" correlations
min_genes_sig <- 10      # min genes for "sig genes only" correlations (usually smaller)

# significant gene set mode: "union" or "intersection"
sig_mode <- "union"

# -----------------------------
# 1) Sanitize guide_res
# -----------------------------
gr <- as.data.frame(guide_res, stringsAsFactors = FALSE)
rownames(gr) <- NULL
needed <- c("gene", "logFC", "FDR", "tested_group")
missing <- setdiff(needed, names(gr))
if (length(missing)) stop("guide_res missing columns: ", paste(missing, collapse = ", "))

gr$gene <- as.character(gr$gene)
gr$tested_group <- as.character(gr$tested_group)
gr$logFC <- as.numeric(gr$logFC)
gr$FDR <- as.numeric(gr$FDR)

# -----------------------------
# 2) Map gRNA -> target (majority rule)
# -----------------------------
g2t <- obj@meta.data %>%
  dplyr::count(gRNA, target, name = "n") %>%
  dplyr::group_by(gRNA) %>%
  dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gRNA, target)

guide2target <- setNames(as.character(g2t$target), as.character(g2t$gRNA))

# -----------------------------
# 3) Build per-guide tables (fast lookup)
# -----------------------------
gr_split <- split(gr, gr$tested_group)

guide_logFC <- lapply(gr_split, function(df) {
  setNames(df$logFC, df$gene)
})

guide_is_sig <- lapply(gr_split, function(df) {
  sig <- (df$FDR <= fdr_thresh) & (abs(df$logFC) >= lfc_thresh)
  setNames(sig, df$gene)
})

all_guides <- names(guide_logFC)
all_targets <- guide2target[all_guides]
keep_guides <- !is.na(all_targets) & all_targets != "NTC" & all_guides != "NTC"
all_guides <- all_guides[keep_guides]
all_targets <- all_targets[keep_guides]

# group guides by target
guides_by_target <- split(all_guides, all_targets)

# -----------------------------
# 4) Correlation functions
# -----------------------------
pair_cor_one <- function(g1, g2,
                         use_sig = FALSE,
                         sig_mode = c("union","intersection"),
                         min_genes = 200) {
  sig_mode <- match.arg(sig_mode)

  v1 <- guide_logFC[[g1]]
  v2 <- guide_logFC[[g2]]
  common <- intersect(names(v1), names(v2))
  if (length(common) < min_genes) return(NULL)

  if (!use_sig) {
    genes_use <- common
  } else {
    s1 <- guide_is_sig[[g1]][common]
    s2 <- guide_is_sig[[g2]][common]
    s1[is.na(s1)] <- FALSE
    s2[is.na(s2)] <- FALSE
    if (sig_mode == "union") {
      genes_use <- common[s1 | s2]
    } else {
      genes_use <- common[s1 & s2]
    }
    if (length(genes_use) < min_genes) return(NULL)
  }

  r <- suppressWarnings(cor(v1[genes_use], v2[genes_use], method = cor_method))
  data.frame(
    guide1 = g1,
    guide2 = g2,
    n_genes = length(genes_use),
    cor = r,
    stringsAsFactors = FALSE
  )
}

compute_correlations <- function(guides_by_target,
                                 use_sig = FALSE,
                                 sig_mode = "union",
                                 min_genes = 200) {
  out <- list()
  k <- 1
  for (tgt in names(guides_by_target)) {
    gs <- guides_by_target[[tgt]]
    if (length(gs) < 2) next

    # pairwise across guides for this target
    for (i in 1:(length(gs) - 1)) {
      for (j in (i + 1):length(gs)) {
        rec <- pair_cor_one(gs[i], gs[j],
                            use_sig = use_sig,
                            sig_mode = sig_mode,
                            min_genes = min_genes)
        if (is.null(rec)) next
        rec$target <- tgt
        out[[k]] <- rec
        k <- k + 1
      }
    }
  }
  if (!length(out)) return(NULL)
  bind_rows(out)
}

# -----------------------------
# 5) Compute correlations: all genes
# -----------------------------
message("Computing correlations using ALL genes (overlap threshold = ", min_genes_all, ") ...")
cor_all <- compute_correlations(
  guides_by_target = guides_by_target,
  use_sig = FALSE,
  min_genes = min_genes_all
)
if (is.null(cor_all) || nrow(cor_all) == 0) {
  stop("No ALL-gene correlations computed. Lower min_genes_all or check guide_res.")
}
write.csv(cor_all, file.path(outdir, "sgRNA_pairwise_cor_allgenes.csv"), row.names = FALSE)

# -----------------------------
# 6) Compute correlations: significant genes only
# -----------------------------
message("Computing correlations using SIGNIFICANT genes only (", sig_mode,
        "; overlap threshold = ", min_genes_sig, ") ...")
cor_sig <- compute_correlations(
  guides_by_target = guides_by_target,
  use_sig = TRUE,
  sig_mode = sig_mode,
  min_genes = min_genes_sig
)
if (is.null(cor_sig) || nrow(cor_sig) == 0) {
  warning("No SIG-gene correlations computed. Try lowering min_genes_sig or using sig_mode='union'.")
} else {
  write.csv(cor_sig,
            file.path(outdir, paste0("sgRNA_pairwise_cor_siggenes_", sig_mode, ".csv")),
            row.names = FALSE)
}

# -----------------------------
# 7) Histograms
# -----------------------------
p1 <- ggplot(cor_all, aes(x = cor)) +
  geom_histogram(bins = 60) +
  theme_classic(base_size = 12) +
  labs(
    title = "Within-target sgRNA concordance (all genes)",
    subtitle = paste0("Pairwise ", cor_method, " correlation of logFC; min_genes = ", min_genes_all),
    x = "Correlation", y = "Number of sgRNA pairs"
  )
ggsave(file.path(outdir, "hist_cor_allgenes.png"), p1, width = 7, height = 5, dpi = 200)

if (!is.null(cor_sig) && nrow(cor_sig) > 0) {
  p2 <- ggplot(cor_sig, aes(x = cor)) +
    geom_histogram(bins = 60) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Within-target sgRNA concordance (significant genes only; ", sig_mode, ")"),
      subtitle = paste0("Sig = FDR≤", fdr_thresh,
                        if (lfc_thresh > 0) paste0(" & |logFC|≥", lfc_thresh) else "",
                        "; min_genes = ", min_genes_sig),
      x = "Correlation", y = "Number of sgRNA pairs"
    )
  ggsave(file.path(outdir, paste0("hist_cor_siggenes_", sig_mode, ".png")),
         p2, width = 7, height = 5, dpi = 200)
}

# -----------------------------
# 8) Target-level summaries (median/mean)
# -----------------------------
summarize_targets <- function(df) {
  df %>%
    group_by(target) %>%
    summarize(
      n_pairs = n(),
      median_cor = median(cor, na.rm = TRUE),
      mean_cor = mean(cor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_cor))
}

sum_all <- summarize_targets(cor_all)
write.csv(sum_all, file.path(outdir, "target_concordance_summary_allgenes.csv"), row.names = FALSE)

if (!is.null(cor_sig) && nrow(cor_sig) > 0) {
  sum_sig <- summarize_targets(cor_sig)
  write.csv(sum_sig, file.path(outdir, paste0("target_concordance_summary_siggenes_", sig_mode, ".csv")),
            row.names = FALSE)
}
```

### Volcano plots and DE genes summaries

```{r}
# -----------------------------
# Volcano plots
# -----------------------------
plot_volcano_fast <- function(df, which_label, prefix, outdir,
                              fdr_cutoff = 0.05, lfc_cutoff = 0.25,
                              label_top_n = 0) {
  x <- df[df$tested_group == which_label, , drop = FALSE]
  if (nrow(x) == 0) return(invisible(NULL))

  # -log10(FDR); avoid Inf
  x$FDR <- as.numeric(x$FDR)
  x$logFC <- as.numeric(x$logFC)
  x$neglog10FDR <- -log10(pmax(x$FDR, .Machine$double.xmin))

  x$sig <- (x$FDR <= fdr_cutoff) & (abs(x$logFC) >= lfc_cutoff)

  p <- ggplot(x, aes(x = logFC, y = neglog10FDR)) +
    geom_point(size = 0.35, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2) +
    labs(
      title = paste0(prefix, " vs NTC: ", which_label),
      x = "logFC", y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12)

  # Optionally label top hits by FDR
  if (label_top_n > 0) {
    top <- x[order(x$FDR), , drop = FALSE]
    top <- head(top, label_top_n)
    if ("gene" %in% names(top)) {
      p <- p + ggrepel::geom_text_repel(
        data = top,
        aes(label = gene),
        size = 2.5,
        max.overlaps = Inf
      )
    }
  }

  fn <- file.path(outdir, paste0("volcano_", prefix, "_", gsub("[^A-Za-z0-9]+","_", which_label), ".png"))
  ggsave(fn, p, width = 6, height = 5, dpi = 200)
  invisible(NULL)
}

  ggsave(
    file.path(outdir,
              paste0("volcano_", prefix, "_", gsub("[^A-Za-z0-9]+", "_", which_label), ".png")),
    width = 6, height = 5, dpi = 200
  )
  invisible(NULL)
}

# Make volcano plot for all sgRNAs:

dir.create(file.path(outdir, "volcano_guides"), showWarnings = FALSE)

gr <- as.data.frame(guide_res, stringsAsFactors = FALSE)
rownames(gr) <- NULL
gr$tested_group <- as.character(gr$tested_group)

all_guides <- sort(unique(gr$tested_group))
all_guides <- setdiff(all_guides, "NTC")   # don’t plot control-vs-control

for (g in all_guides) {
  plot_volcano_fast(
    df = gr,
    which_label = g,
    prefix = "guide",
    outdir = file.path(outdir, "volcano_guides"),
    label_top_n = 0   # set to e.g. 10 if you want labels
  )
}

# Make volcano plot for all targets
dir.create(file.path(outdir, "volcano_targets"), showWarnings = FALSE)

tr <- as.data.frame(target_res, stringsAsFactors = FALSE)
rownames(tr) <- NULL
tr$tested_group <- as.character(tr$tested_group)

all_targets <- sort(unique(tr$tested_group))
all_targets <- setdiff(all_targets, "NTC")

for (tgt in all_targets) {
  plot_volcano_fast(
    df = tr,
    which_label = tgt,
    prefix = "target",
    outdir = file.path(outdir, "volcano_targets"),
    label_top_n = 0
  )
}

# Summary tables
# per sgRNA summary
guide_summary <- gr %>%
  group_by(tested_group) %>%
  summarize(
    n_genes_FDR05 = sum(FDR <= 0.05, na.rm = TRUE),
    minFDR = min(FDR, na.rm = TRUE),
    top_gene = gene[which.min(FDR)],
    .groups = "drop"
  ) %>%
  arrange(minFDR)

write.csv(guide_summary, file.path(outdir, "guide_summary.csv"), row.names = FALSE)

# per target summary
target_summary <- tr %>%
  group_by(tested_group) %>%
  summarize(
    n_genes_FDR05 = sum(FDR <= 0.05, na.rm = TRUE),
    minFDR = min(FDR, na.rm = TRUE),
    top_gene = gene[which.min(FDR)],
    .groups = "drop"
  ) %>%
  arrange(minFDR)

write.csv(target_summary, file.path(outdir, "target_summary.csv"), row.names = FALSE)
```

## ------------------------------------------------------------
## STEP 4: Power / detectability diagnostics for sgRNA effects
## ------------------------------------------------------------

```{r}
# -----------------------------
# SETTINGS
# -----------------------------
outdir_power <- file.path(outdir, "power_diagnostics")
dir.create(outdir_power, showWarnings = FALSE)

fdr_cut <- 0.05
lfc_cut <- 0.25         # used only for "big effect" flag; keep 0.25 or set 0
bin_width <- 0.05       # |logFC| bins for detection curve

# For downsampling check:
do_downsample_check <- TRUE
downsample_guides_n <- 20     # number of guides to test (top by minFDR); set Inf for all (not recommended)
ntc_downsample_sizes <- c(100, 200, 500, 1000, 5000)  # compare to full NTC
downsample_seed <- 1

# NOTE: downsample check requires you have:
#   - obj (Seurat)
#   - counts (RenameDims/BPCells counts used earlier)
#   - run_pseudorep_edger_one() from the pseudo-rep script
# If you don't, set do_downsample_check <- FALSE

# -----------------------------
# INPUT: guide_res
# -----------------------------
gr <- as.data.frame(guide_res, stringsAsFactors = FALSE)
rownames(gr) <- NULL

# enforce types
gr$tested_group <- as.character(gr$tested_group)
gr$gene <- as.character(gr$gene)
gr$logFC <- as.numeric(gr$logFC)
gr$FDR <- as.numeric(gr$FDR)

# remove NTC if present
gr <- gr[gr$tested_group != "NTC", , drop = FALSE]

# -----------------------------
# 1) Effect-size summaries per sgRNA
# -----------------------------
effect_summary <- gr %>%
  group_by(tested_group) %>%
  summarize(
    n_genes = dplyr::n(),
    n_cells = dplyr::first(n_cells_group),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    p90_abs_logFC = as.numeric(quantile(abs(logFC), 0.90, na.rm = TRUE)),
    p95_abs_logFC = as.numeric(quantile(abs(logFC), 0.95, na.rm = TRUE)),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    minFDR = min(FDR, na.rm = TRUE),
    n_sig_FDR05 = sum(FDR <= fdr_cut, na.rm = TRUE),
    n_big = sum((FDR <= fdr_cut) & (abs(logFC) >= lfc_cut), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(minFDR)

write.csv(effect_summary, file.path(outdir_power, "guide_effect_summary.csv"), row.names = FALSE)

# Histograms
p_med <- ggplot(effect_summary, aes(x = median_abs_logFC)) +
  geom_histogram(bins = 60) +
  theme_classic(base_size = 12) +
  labs(title = "Per-sgRNA effect size (median |logFC| across genes)",
       x = "Median |logFC|", y = "Number of sgRNAs")
ggsave(file.path(outdir_power, "hist_median_abs_logFC.png"), p_med, width = 7, height = 5, dpi = 200)

p_p90 <- ggplot(effect_summary, aes(x = p90_abs_logFC)) +
  geom_histogram(bins = 60) +
  theme_classic(base_size = 12) +
  labs(title = "Per-sgRNA effect size (90th percentile |logFC| across genes)",
       x = "P90 |logFC|", y = "Number of sgRNAs")
ggsave(file.path(outdir_power, "hist_p90_abs_logFC.png"), p_p90, width = 7, height = 5, dpi = 200)


# -----------------------------
# 2) Optional: NTC downsampling sensitivity check
# -----------------------------
# This answers: "Are we only significant because NTC is huge?"
# For a set of guides, re-run the pseudo-rep test using different NTC sizes
# and track minFDR and #sig genes.

if (isTRUE(do_downsample_check)) {
  needed_objs <- c("obj", "counts", "run_pseudorep_edger_one", "cells_by_gRNA")
  missing2 <- needed_objs[!vapply(needed_objs, exists, logical(1))]
  if (length(missing2)) {
    warning("Skipping downsample check; missing objects/functions: ",
            paste(missing2, collapse = ", "),
            "\nSet do_downsample_check <- FALSE or ensure these exist.")
  } else {
    message("Running NTC downsampling check (guide subset = ", downsample_guides_n, ") ...")

    # guides to test: take top by minFDR (strongest effects)
    guides_to_test <- effect_summary$tested_group
    if (is.finite(downsample_guides_n)) {
      guides_to_test <- head(guides_to_test, downsample_guides_n)
    }

    ntc_full <- cells_by_gRNA[["NTC"]]
    if (is.null(ntc_full)) stop("cells_by_gRNA[['NTC']] not found")

    set.seed(downsample_seed)

    down_rows <- list()
    kk <- 1
    for (g in guides_to_test) {
      grp_idx <- cells_by_gRNA[[g]]
      if (is.null(grp_idx)) next

      for (n_ntc in ntc_downsample_sizes) {
        if (n_ntc > length(ntc_full)) next

        ntc_idx <- sample(ntc_full, n_ntc)

        tt <- run_pseudorep_edger_one(
          group_name = g,
          group_idx = grp_idx,
          ntc_idx = ntc_idx,
          counts_obj = counts,
          K = K,
          min_cells_per_rep = min_cells_per_rep,
          seed = downsample_seed + 10000
        )

        if (is.null(tt) || nrow(tt) == 0) next

        tt <- as.data.frame(tt, stringsAsFactors = FALSE)
        tt$FDR <- as.numeric(tt$FDR)
        tt$logFC <- as.numeric(tt$logFC)

        down_rows[[kk]] <- data.frame(
          tested_group = g,
          n_cells_group = length(grp_idx),
          n_cells_ntc = n_ntc,
          minFDR = min(tt$FDR, na.rm = TRUE),
          n_sig = sum(tt$FDR <= fdr_cut, na.rm = TRUE),
          n_big = sum((tt$FDR <= fdr_cut) & (abs(tt$logFC) >= lfc_cut), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        kk <- kk + 1
        message("  guide=", g, " NTC=", n_ntc, " done")
      }
    }

    down_df <- if (length(down_rows)) bind_rows(down_rows) else NULL
    if (!is.null(down_df) && nrow(down_df) > 0) {
      write.csv(down_df, file.path(outdir_power, "ntc_downsampling_check.csv"), row.names = FALSE)

      p_ds1 <- ggplot(down_df, aes(x = n_cells_ntc, y = -log10(pmax(minFDR, .Machine$double.xmin)),
                                   group = tested_group)) +
        geom_line(alpha = 0.5) +
        geom_point(alpha = 0.7) +
        theme_classic(base_size = 12) +
        labs(title = "NTC downsampling: best significance vs NTC size",
             x = "NTC cells used", y = "-log10(min FDR)")
      ggsave(file.path(outdir_power, "ntc_downsampling_minFDR.png"),
             p_ds1, width = 7, height = 5, dpi = 200)

      p_ds2 <- ggplot(down_df, aes(x = n_cells_ntc, y = n_sig, group = tested_group)) +
        geom_line(alpha = 0.5) +
        geom_point(alpha = 0.7) +
        theme_classic(base_size = 12) +
        labs(title = "NTC downsampling: # significant genes vs NTC size",
             x = "NTC cells used", y = paste0("# genes with FDR ≤ ", fdr_cut))
      ggsave(file.path(outdir_power, "ntc_downsampling_nSig.png"),
             p_ds2, width = 7, height = 5, dpi = 200)
    } else {
      warning("Downsampling check produced no results (try fewer constraints or verify run_pseudorep_edger_one).")
    }
  }
}

cat("\nDONE. Wrote diagnostics to:\n", normalizePath(outdir_power), "\n",
    "Key files:\n",
    " - guide_effect_summary.csv\n",
    " - hist_median_abs_logFC.png, hist_p90_abs_logFC.png\n",
    " - detection_curve_fracSig_vs_absLogFC.csv, curve_fracSig_vs_absLogFC.png\n",
    " - scatter_cells_vs_medianAbsLogFC.png\n",
    " - scatter_cells_vs_nSig.png\n",
    " - scatter_cells_vs_neglog10minFDR.png\n",
    if (do_downsample_check) " - ntc_downsampling_check.csv + plots (if computed)\n" else "",
    sep = "")

  ```

### Additional evaluations

  ```{r}
# ============================================================
# Follow-up analyses:
#   A) Cell-count–adjust effect sizes and re-rank sgRNAs
#   B) Overlay guide concordance on the effect-vs-cells plot
#
# Assumes you have:
#   - guide_res (edgeR pseudo-rep results) with: tested_group, gene, logFC, FDR, n_cells_group
#   - obj (Seurat object) for mapping gRNA -> target
#
# Outputs written to: outdir_followup
# ============================================================

outdir_followup <- file.path(outdir, "followup_sgRNA_assessment")
dir.create(outdir_followup, showWarnings = FALSE)

# -----------------------------
# 0) Sanitize guide_res and compute per-sgRNA summaries
# -----------------------------
stopifnot(exists("guide_res"))
stopifnot(exists("obj"))

gr <- as.data.frame(guide_res, stringsAsFactors = FALSE)
rownames(gr) <- NULL
gr$tested_group <- as.character(gr$tested_group)
gr$gene <- as.character(gr$gene)
gr$logFC <- as.numeric(gr$logFC)
gr$FDR <- as.numeric(gr$FDR)
gr$n_cells_group <- as.integer(gr$n_cells_group)

gr <- gr[gr$tested_group != "NTC", , drop = FALSE]

# per-sgRNA summary metrics
guide_summary <- gr %>%
  group_by(tested_group) %>%
  summarize(
    n_cells = first(n_cells_group),
    n_genes = n(),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    p90_abs_logFC = as.numeric(quantile(abs(logFC), 0.90, na.rm = TRUE)),
    p95_abs_logFC = as.numeric(quantile(abs(logFC), 0.95, na.rm = TRUE)),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    minFDR = min(FDR, na.rm = TRUE),
    n_sig_FDR05 = sum(FDR <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(guide_summary, file.path(outdir_followup, "guide_summary_raw.csv"), row.names = FALSE)

# -----------------------------
# A) Correct effect sizes for cell count and re-rank sgRNAs
# -----------------------------
# We adjust median_abs_logFC for cell number by regressing on log(n_cells).
# The residual is: "stronger/weaker than expected given its cell count".

fit_adj <- lm(median_abs_logFC ~ log(n_cells), data = guide_summary)
guide_summary$adj_median_abs_logFC <- residuals(fit_adj)

# Also adjust p90 if you prefer a more tail-based metric
fit_adj90 <- lm(p90_abs_logFC ~ log(n_cells), data = guide_summary)
guide_summary$adj_p90_abs_logFC <- residuals(fit_adj90)

# Re-rank (descending = strongest beyond expectation)
guide_ranked <- guide_summary %>%
  arrange(desc(adj_median_abs_logFC))

write.csv(guide_ranked, file.path(outdir_followup, "guide_ranked_adjusted.csv"), row.names = FALSE)

# Plot: raw vs adjusted
p_adj <- ggplot(guide_summary, aes(x = log(n_cells), y = median_abs_logFC)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "Median |logFC| vs log(cells) with linear fit",
       x = "log(cells per sgRNA)", y = "Median |logFC|")
ggsave(file.path(outdir_followup, "lm_fit_medianAbsLogFC_vs_logCells.png"),
       p_adj, width = 7, height = 5, dpi = 200)

p_resid <- ggplot(guide_ranked, aes(x = adj_median_abs_logFC)) +
  geom_histogram(bins = 60) +
  theme_classic(base_size = 12) +
  labs(title = "Cell-count–adjusted sgRNA effect strength (residuals)",
       x = "Residual of median |logFC| after regressing on log(cells)",
       y = "Number of sgRNAs")
ggsave(file.path(outdir_followup, "hist_adjusted_effect_residuals.png"),
       p_resid, width = 7, height = 5, dpi = 200)

# -----------------------------
# Helper: map gRNA -> target (majority rule)
# -----------------------------
g2t <- obj@meta.data %>%
  dplyr::count(gRNA, target, name = "n") %>%
  dplyr::group_by(gRNA) %>%
  dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gRNA, target)

guide2target <- setNames(as.character(g2t$target), as.character(g2t$gRNA))
guide_summary$target <- guide2target[guide_summary$tested_group]
guide_summary$target <- as.character(guide_summary$target)

# -----------------------------
# B) Overlay guide concordance on the effect-vs-cells plot
# -----------------------------
# Concordance definition:
#   For each target with >=2 guides, compute all pairwise correlations
#   between guides' logFC vectors across genes, then assign each guide
#   its median correlation to sibling guides.
#
# This is intended as a QC metric: higher = more consistent perturbation signal.

# Build per-guide logFC vectors (named by gene)
gr_split <- split(gr, gr$tested_group)
guide_logFC <- lapply(gr_split, function(df) {
  setNames(df$logFC, df$gene)
})

pair_cor <- function(g1, g2, min_genes = 200, method = "pearson") {
  v1 <- guide_logFC[[g1]]
  v2 <- guide_logFC[[g2]]
  common <- intersect(names(v1), names(v2))
  if (length(common) < min_genes) return(NA_real_)
  suppressWarnings(cor(v1[common], v2[common], method = method))
}

min_genes_cor <- 200

# compute per-guide median correlation to sibling guides
guide_median_cor <- rep(NA_real_, length(guide_logFC))
names(guide_median_cor) <- names(guide_logFC)

targets <- guide2target[names(guide_logFC)]
ok <- !is.na(targets) & targets != "NTC"
targets <- targets[ok]
guides_ok <- names(targets)

guides_by_target <- split(guides_ok, targets)

for (tgt in names(guides_by_target)) {
  gs <- guides_by_target[[tgt]]
  if (length(gs) < 2) next
  # pairwise correlations
  cors <- list()
  for (i in seq_along(gs)) {
    gi <- gs[i]
    other <- gs[-i]
    vals <- vapply(other, function(gj) pair_cor(gi, gj, min_genes = min_genes_cor), numeric(1))
    guide_median_cor[gi] <- median(vals, na.rm = TRUE)
  }
}

guide_summary$concordance_medCor <- guide_median_cor[guide_summary$tested_group]

write.csv(guide_summary, file.path(outdir_followup, "guide_summary_with_concordance.csv"),
          row.names = FALSE)

# Plot: effect vs cells colored by concordance
p_conc <- ggplot(guide_summary, aes(x = n_cells, y = median_abs_logFC, color = concordance_medCor)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  theme_classic(base_size = 12) +
  labs(title = "Effect size vs cells per sgRNA (colored by within-target concordance)",
       subtitle = paste0("Concordance = median Pearson r to sibling guides; min_genes=", min_genes_cor),
       x = "Cells per sgRNA", y = "Median |logFC| across genes", color = "Median r") +
  scale_color_viridis_c(option = "C", na.value = "grey70")
ggsave(file.path(outdir_followup, "scatter_cells_vs_medianAbsLogFC_coloredByConcordance.png"),
       p_conc, width = 7, height = 5, dpi = 200)

```