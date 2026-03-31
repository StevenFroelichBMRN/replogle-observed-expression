# Evaluate knockdowns and potential impact in genes of interest

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
# Load unfiltered gene x cell matrix object
scobj_filt <- readRDS("analyses/2026_PerturbSeq_Mega_druggable/scobj_filt.rds")

# Normalize and scale:
scobj_filt <- NormalizeData(scobj_filt, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj_filt <- FindVariableFeatures(scobj_filt, selection.method = "vst", nfeatures = 3000)
scobj_filt <- ScaleData(scobj_filt)
```

```{r}
############################################################
## Perturb-seq promoter targeting analysis in Seurat v5
## Object: scobj_filt
##
## Goals:
## 1) Did the gRNAs work?
##    - gRNA-focused
##    - agreement across gRNAs targeting same promoter
##    - promoter-focused
##    - gene-focused
##
## 2) Did any perturbation change expression of genes of interest?
##    - using a CSV with no header
##
## Notes:
## - Uses Wilcoxon tests on per-cell expression
## - Uses Seurat v5 syntax: layer = "data"
## - Keeps code explicit/simple, without many helper functions
############################################################

############################################################
## 0) INPUTS
############################################################

# Seurat object
seu <- scobj_filt

# Metadata columns
gRNA_col <- "gRNA"
target_col <- "target"

# File with genes of interest (NO HEADER)
goi_file <- "analyses/2026_PerturbSeq_Mega_druggable/target_touch_genes.csv"

# Analysis settings
min_cells <- 15
ntc_pattern <- "^Non_Targeting_Human_CRi_"

# Assay / layer
DefaultAssay(seu) <- "RNA"
expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")

############################################################
## 1) PREP METADATA
############################################################

meta <- seu@meta.data
meta$cell_id <- rownames(meta)

stopifnot(gRNA_col %in% colnames(meta))
stopifnot(target_col %in% colnames(meta))

meta$gRNA_name <- meta[[gRNA_col]]
meta$target_name <- meta[[target_col]]

# Flag NTCs
meta$is_ntc <- grepl(ntc_pattern, meta$target_name)
meta$target_collapsed <- ifelse(meta$is_ntc, "NTC", meta$target_name)

# Parse target gene and promoter from target column
# Example: CACNA1B_P1 -> gene = CACNA1B ; promoter = CACNA1B_P1
meta$target_gene <- ifelse(
  meta$is_ntc,
  NA,
  sub("_P[0-9]+$", "", meta$target_name)
)

meta$target_promoter <- ifelse(
  meta$is_ntc,
  NA,
  meta$target_name
)

# Parse gene and promoter from gRNA column
# Example: CACNA1B_P1_4 -> promoter = CACNA1B_P1 ; gene = CACNA1B
meta$gRNA_promoter <- ifelse(
  meta$is_ntc,
  NA,
  sub("_[0-9]+$", "", meta$gRNA_name)
)

meta$gRNA_gene <- ifelse(
  meta$is_ntc,
  NA,
  sub("_P[0-9]+_[0-9]+$", "", meta$gRNA_name)
)

# Add back to object
seu <- AddMetaData(
  seu,
  metadata = meta[, c(
    "gRNA_name", "target_name", "is_ntc", "target_collapsed",
    "target_gene", "target_promoter", "gRNA_promoter", "gRNA_gene"
  )]
)

# NTC cells
ntc_cells <- rownames(meta)[meta$is_ntc]

############################################################
## 2) LOAD GENES OF INTEREST (CSV HAS NO HEADER)
############################################################

goi_df <- read.csv(goi_file, header = FALSE, stringsAsFactors = FALSE)
genes_of_interest <- unique(goi_df[[1]])
genes_of_interest <- genes_of_interest[!is.na(genes_of_interest)]
genes_of_interest <- genes_of_interest[genes_of_interest != ""]

genes_of_interest_present <- intersect(genes_of_interest, rownames(expr_mat))

############################################################
## 3) QUICK OVERVIEW OF CELL COUNTS
############################################################

cells_per_gRNA <- meta %>%
  count(gRNA_name, sort = TRUE)

cells_per_promoter <- meta %>%
  filter(!is_ntc) %>%
  count(target_promoter, sort = TRUE)

cells_per_gene <- meta %>%
  filter(!is_ntc) %>%
  count(target_gene, sort = TRUE)

p_cells_gRNA <- ggplot(cells_per_gRNA, aes(x = reorder(gRNA_name, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Cells per gRNA", x = "gRNA", y = "Number of cells")

p_cells_promoter <- ggplot(cells_per_promoter, aes(x = reorder(target_promoter, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Cells per promoter", x = "Promoter", y = "Number of cells")

p_cells_gene <- ggplot(cells_per_gene, aes(x = reorder(target_gene, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Cells per targeted gene", x = "Gene", y = "Number of cells")

print(p_cells_gRNA)
print(p_cells_promoter)
print(p_cells_gene)

############################################################
## 4) gRNA-FOCUSED KNOCKDOWN ASSESSMENT
## Compare each guide against NTC for the corresponding target gene
############################################################

gRNA_results <- data.frame()

all_guides <- unique(meta$gRNA_name[!meta$is_ntc])

for (gr in all_guides) {

  guide_cells <- rownames(meta)[meta$gRNA_name == gr]

  if (length(guide_cells) < min_cells) next

  promoter_here <- unique(meta$gRNA_promoter[meta$gRNA_name == gr])
  gene_here <- unique(meta$gRNA_gene[meta$gRNA_name == gr])

  if (length(promoter_here) != 1) next
  if (length(gene_here) != 1) next
  if (is.na(gene_here)) next
  if (!(gene_here %in% rownames(expr_mat))) next

  x <- as.numeric(Matrix::as.matrix(expr_mat[gene_here, guide_cells, drop = FALSE]))
  y <- as.numeric(Matrix::as.matrix(expr_mat[gene_here, ntc_cells, drop = FALSE]))

  wt <- wilcox.test(x, y)

 gRNA_results <- rbind(
    gRNA_results,
    data.frame(
      gRNA = gr,
      promoter = promoter_here,
      target_gene = gene_here,
      n_gRNA = length(guide_cells),
      n_ntc = length(ntc_cells),
      mean_expr_gRNA = mean(x),
      mean_expr_ntc = mean(y),
      median_expr_gRNA = median(x),
      median_expr_ntc = median(y),
      logFC_mean = mean(x) - mean(y),
      logFC_median = median(x) - median(y),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  )
}

gRNA_results$p_adj <- p.adjust(gRNA_results$p_value, method = "fdr")
gRNA_results$significant_kd <- gRNA_results$p_adj < 0.1 & gRNA_results$logFC_mean < 0
gRNA_results <- gRNA_results %>% arrange(p_adj, logFC_mean)

write.csv(gRNA_results, "analyses/2026_PerturbSeq_Mega_druggable/gRNA_target_gene_knockdown_results.csv", row.names = FALSE)

############################################################
## 5) VISUALIZE gRNA-LEVEL KNOCKDOWN
############################################################

p_gRNA_volcano <- ggplot(gRNA_results, aes(x = logFC_mean, y = -log10(p_adj))) +
  geom_point(aes(color = significant_kd), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  theme_bw() +
  labs(
    title = "gRNA-focused knockdown assessment",
    x = "Mean expression difference vs NTC",
    y = "-log10(FDR)"
  )
ggsave("analyses/2026_PerturbSeq_Mega_druggable/gRNA_knockdown_volcano.pdf", p_gRNA_volcano, width = 6, height = 5)

############################################################
## 6) AGREEMENT BETWEEN gRNAs TARGETING THE SAME PROMOTER
############################################################

guide_agreement <- data.frame()

if (nrow(gRNA_results) > 0) {
  guide_agreement <- gRNA_results %>%
    group_by(promoter, target_gene) %>%
    summarise(
      n_guides_tested = n(),
      n_sig_kd = sum(significant_kd),
      frac_sig_kd = mean(significant_kd),
      mean_logFC = mean(logFC_mean),
      sd_logFC = sd(logFC_mean),
      all_negative = all(logFC_mean < 0),
      .groups = "drop"
    ) %>%
    arrange(desc(frac_sig_kd), mean_logFC)
}

write.csv(guide_agreement, "analyses/2026_PerturbSeq_Mega_druggable/gRNA_agreement_by_promoter.csv", row.names = FALSE)

############################################################
## 7) PROMOTER-FOCUSED KNOCKDOWN ASSESSMENT
## Combine all cells for guides targeting the same promoter
############################################################

promoter_results <- data.frame()

all_promoters <- unique(meta$target_promoter[!meta$is_ntc])

for (pr in all_promoters) {

  idx <- !is.na(meta$target_promoter) & meta$target_promoter == pr
  pert_cells <- rownames(meta)[idx]

  if (length(pert_cells) < min_cells) next

  gene_here <- unique(meta$target_gene[idx])

  if (length(gene_here) != 1) next
  if (is.na(gene_here)) next
  if (!(gene_here %in% rownames(expr_mat))) next

  x <- as.numeric(Matrix::as.matrix(expr_mat[gene_here, pert_cells, drop = FALSE]))
  y <- as.numeric(Matrix::as.matrix(expr_mat[gene_here, ntc_cells, drop = FALSE]))

  wt <- wilcox.test(x, y)

  promoter_results <- rbind(
    promoter_results,
    data.frame(
      promoter = pr,
      target_gene = gene_here,
      n_promoter = length(pert_cells),
      n_ntc = length(ntc_cells),
      mean_expr_promoter = mean(x),
      mean_expr_ntc = mean(y),
      median_expr_promoter = median(x),
      median_expr_ntc = median(y),
      logFC_mean = mean(x) - mean(y),
      logFC_median = median(x) - median(y),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  )
}
promoter_results$p_adj <- p.adjust(promoter_results$p_value, method = "fdr")
promoter_results$significant_kd <- promoter_results$p_adj < 0.1 & promoter_results$logFC_mean < 0
promoter_results <- promoter_results %>% arrange(p_adj, logFC_mean)

write.csv(promoter_results, "analyses/2026_PerturbSeq_Mega_druggable/promoter_target_gene_knockdown_results.csv", row.names = FALSE)

p_promoter_volcano <- ggplot(promoter_results, aes(x = logFC_mean, y = -log10(p_adj))) +
    geom_point(aes(color = significant_kd), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
    theme_bw() +
    labs(
      title = "Promoter-focused knockdown assessment",
      x = "Mean expression difference vs NTC",
      y = "-log10(FDR)"
    )

ggsave("analyses/2026_PerturbSeq_Mega_druggable/promoter_knockdown_volcano.pdf", p_promoter_volcano, width = 6, height = 5)

############################################################
## 8) GENE-FOCUSED KNOCKDOWN ASSESSMENT
## Combine all cells targeting promoters of the same gene
############################################################

gene_results <- data.frame()

all_target_genes <- unique(meta$target_gene[!meta$is_ntc])

for (gn in all_target_genes) {

  idx <- !is.na(meta$target_gene) & meta$target_gene == gn
  pert_cells <- unique(rownames(meta)[idx])

  if (length(pert_cells) < min_cells) next
  if (!(gn %in% rownames(expr_mat))) next

  x <- as.numeric(Matrix::as.matrix(expr_mat[gn, pert_cells, drop = FALSE]))
  y <- as.numeric(Matrix::as.matrix(expr_mat[gn, ntc_cells, drop = FALSE]))

  wt <- wilcox.test(x, y)

  gene_results <- rbind(
    gene_results,
    data.frame(
      target_gene = gn,
      n_gene = length(pert_cells),
      n_ntc = length(ntc_cells),
      mean_expr_gene = mean(x),
      mean_expr_ntc = mean(y),
      median_expr_gene = median(x),
      median_expr_ntc = median(y),
      logFC_mean = mean(x) - mean(y),
      logFC_median = median(x) - median(y),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  )
}

gene_results$p_adj <- p.adjust(gene_results$p_value, method = "fdr")
gene_results$significant_kd <- gene_results$p_adj < 0.05 & gene_results$logFC_mean < 0
gene_results <- gene_results %>% arrange(p_adj, logFC_mean)

write.csv(gene_results, "analyses/2026_PerturbSeq_Mega_druggable/gene_target_gene_knockdown_results.csv", row.names = FALSE)

p_gene_volcano <- ggplot(gene_results, aes(x = logFC_mean, y = -log10(p_adj))) +
    geom_point(aes(color = significant_kd), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
    theme_bw() +
    labs(
      title = "Gene-focused knockdown assessment",
      x = "Mean expression difference vs NTC",
      y = "-log10(FDR)"
  )

ggsave("analyses/2026_PerturbSeq_Mega_druggable/gene_knockdown_volcano.pdf", p_gene_volcano, width = 6, height = 5)

############################################################
## 9) TEST EFFECTS ON GENES OF INTEREST
## Main analysis: promoter-level perturbation vs NTC
############################################################

goi_results <- data.frame()

perturbations_to_test <- unique(
  meta$target_promoter[!is.na(meta$target_promoter) & !meta$is_ntc]
)

ntc_use <- intersect(
  unique(rownames(meta)[meta$is_ntc]),
  colnames(expr_mat)
)

for (pr in perturbations_to_test) {

  idx <- !is.na(meta$target_promoter) & meta$target_promoter == pr

  pert_cells <- intersect(
    unique(rownames(meta)[idx]),
    colnames(expr_mat)
  )

  if (length(pert_cells) < min_cells) next

  for (gene_test in genes_of_interest_present) {

    if (!(gene_test %in% rownames(expr_mat))) next

    x <- as.numeric(Matrix::as.matrix(
      expr_mat[gene_test, pert_cells, drop = FALSE]
    ))

    y <- as.numeric(Matrix::as.matrix(
      expr_mat[gene_test, ntc_use, drop = FALSE]
    ))

    wt <- wilcox.test(x, y)

    goi_results <- rbind(
      goi_results,
      data.frame(
        perturbation = pr,
        tested_gene = gene_test,
        n_perturb = length(pert_cells),
        n_ntc = length(ntc_use),
        mean_expr_perturb = mean(x),
        mean_expr_ntc = mean(y),
        median_expr_perturb = median(x),
        median_expr_ntc = median(y),
        logFC_mean = mean(x) - mean(y),
        logFC_median = median(x) - median(y),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}


goi_results$p_adj <- p.adjust(goi_results$p_value, method = "fdr")
goi_results$significant <- goi_results$p_adj < 0.05
goi_results$direction <- ifelse(goi_results$logFC_mean > 0, "up", "down")
goi_results <- goi_results %>% arrange(p_adj)

write.csv(goi_results, "promoter_vs_genes_of_interest_results.csv", row.names = FALSE)

############################################################
## 11) VISUALIZE GENES-OF-INTEREST RESULTS
############################################################

if (nrow(goi_results) > 0) {

  # Full effect-size heatmap
  goi_heat_df <- goi_results %>%
    select(perturbation, tested_gene, logFC_mean) %>%
    tidyr::pivot_wider(names_from = tested_gene, values_from = logFC_mean)

  if (nrow(goi_heat_df) > 0) {
    goi_mat <- as.matrix(goi_heat_df[, -1])
    rownames(goi_mat) <- goi_heat_df$perturbation

    pdf("promoter_vs_genes_of_interest_heatmap.pdf", width = 12, height = 10)
    pheatmap(
      goi_mat,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      border_color = NA,
      main = "Promoter perturbation effects on genes of interest"
    )
    dev.off()
  }

  # Significant-only heatmap: non-significant set to 0
  goi_sig_df <- goi_results %>%
    mutate(signed_sig = ifelse(significant, logFC_mean, 0)) %>%
    select(perturbation, tested_gene, signed_sig) %>%
    tidyr::pivot_wider(names_from = tested_gene, values_from = signed_sig)

  if (nrow(goi_sig_df) > 0) {
    goi_sig_mat <- as.matrix(goi_sig_df[, -1])
    rownames(goi_sig_mat) <- goi_sig_df$perturbation

    pdf("promoter_vs_genes_of_interest_significant_heatmap.pdf", width = 12, height = 10)
    pheatmap(
      goi_sig_mat,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      border_color = NA,
      main = "Significant promoter effects on genes of interest (NS = 0)"
    )
    dev.off()
  }

  # Number of significant GOI hits per perturbation
  goi_hits_per_pert <- goi_results %>%
    group_by(perturbation) %>%
    summarise(
      n_sig_hits = sum(significant),
      mean_abs_effect = mean(abs(logFC_mean)),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_hits))

  write.csv(goi_hits_per_pert, "promoter_genes_of_interest_hit_counts.csv", row.names = FALSE)

  p_goi_hits <- ggplot(goi_hits_per_pert, aes(x = reorder(perturbation, n_sig_hits), y = n_sig_hits)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Number of significant genes-of-interest changes per perturbation",
      x = "Perturbation",
      y = "Number of significant genes"
    )

  print(p_goi_hits)
  ggsave("promoter_genes_of_interest_hit_counts.pdf", p_goi_hits, width = 7, height = 10)

  # Global volcano-like plot
  top_goi_changes <- goi_results %>%
    filter(significant) %>%
    arrange(p_adj) %>%
    slice_head(n = min(30, n()))

  p_top_goi <- ggplot(goi_results, aes(x = logFC_mean, y = -log10(p_adj))) +
    geom_point(alpha = 0.35) +
    geom_point(data = subset(goi_results, significant), color = "red", alpha = 0.7) +
    geom_text_repel(
      data = top_goi_changes,
      aes(label = paste(perturbation, tested_gene, sep = " -> ")),
      size = 3,
      max.overlaps = 30
    ) +
    theme_bw() +
    labs(
      title = "Promoter perturbation effects on genes of interest",
      x = "Mean expression difference vs NTC",
      y = "-log10(FDR)"
    )

  print(p_top_goi)
  ggsave("promoter_genes_of_interest_volcano.pdf", p_top_goi, width = 8, height = 6)
}

############################################################
## 12) OPTIONAL: GENE-LEVEL PERTURBATION ANALYSIS FOR GOI
## Combine all promoters for the same gene
############################################################

goi_gene_level_results <- data.frame()

for (gn in unique(meta$target_gene[!meta$is_ntc])) {

  pert_cells <- rownames(meta)[meta$target_gene == gn]

  if (length(pert_cells) < min_cells) next

  for (gene_test in genes_of_interest_present) {

    x <- as.numeric(expr_mat[gene_test, pert_cells])
    y <- as.numeric(expr_mat[gene_test, ntc_cells])

    wt <- wilcox.test(x, y)

    goi_gene_level_results <- rbind(
      goi_gene_level_results,
      data.frame(
        perturbation_gene = gn,
        tested_gene = gene_test,
        n_perturb = length(pert_cells),
        n_ntc = length(ntc_cells),
        mean_expr_perturb = mean(x),
        mean_expr_ntc = mean(y),
        logFC_mean = mean(x) - mean(y),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}

if (nrow(goi_gene_level_results) > 0) {
  goi_gene_level_results$p_adj <- p.adjust(goi_gene_level_results$p_value, method = "fdr")
  goi_gene_level_results$significant <- goi_gene_level_results$p_adj < 0.05
  goi_gene_level_results <- goi_gene_level_results %>% arrange(p_adj)
}

write.csv(goi_gene_level_results, "gene_level_vs_genes_of_interest_results.csv", row.names = FALSE)

############################################################
## 13) SUMMARY TABLES
############################################################

summary_table <- data.frame(
  metric = c(
    "n_cells_total",
    "n_ntc_cells",
    "n_guides_tested",
    "n_promoters_tested",
    "n_target_genes_tested",
    "n_significant_guide_knockdowns",
    "n_significant_promoter_knockdowns",
    "n_significant_gene_knockdowns",
    "n_genes_of_interest_in_file",
    "n_genes_of_interest_present",
    "n_significant_promoter_GOI_pairs",
    "n_significant_gene_GOI_pairs"
  ),
  value = c(
    nrow(meta),
    length(ntc_cells),
    nrow(gRNA_results),
    nrow(promoter_results),
    nrow(gene_results),
    ifelse(nrow(gRNA_results) > 0, sum(gRNA_results$significant_kd), 0),
    ifelse(nrow(promoter_results) > 0, sum(promoter_results$significant_kd), 0),
    ifelse(nrow(gene_results) > 0, sum(gene_results$significant_kd), 0),
    length(genes_of_interest),
    length(genes_of_interest_present),
    ifelse(nrow(goi_results) > 0, sum(goi_results$significant), 0),
    ifelse(nrow(goi_gene_level_results) > 0, sum(goi_gene_level_results$significant), 0)
  )
)

print(summary_table)
write.csv(summary_table, "analysis_summary_table.csv", row.names = FALSE)

############################################################
## 14) OPTIONAL: SAVE TOP TABLES
############################################################

if (nrow(gRNA_results) > 0) {
  write.csv(
    gRNA_results %>% arrange(p_adj, logFC_mean) %>% slice_head(n = 50),
    "top_50_gRNA_knockdown_hits.csv",
    row.names = FALSE
  )
}

if (nrow(promoter_results) > 0) {
  write.csv(
    promoter_results %>% arrange(p_adj, logFC_mean) %>% slice_head(n = 50),
    "top_50_promoter_knockdown_hits.csv",
    row.names = FALSE
  )
}

if (nrow(gene_results) > 0) {
  write.csv(
    gene_results %>% arrange(p_adj, logFC_mean) %>% slice_head(n = 50),
    "top_50_gene_knockdown_hits.csv",
    row.names = FALSE
  )
}

if (nrow(goi_results) > 0) {
  write.csv(
    goi_results %>% arrange(p_adj) %>% slice_head(n = 100),
    "top_100_promoter_GOI_hits.csv",
    row.names = FALSE
  )
}




```