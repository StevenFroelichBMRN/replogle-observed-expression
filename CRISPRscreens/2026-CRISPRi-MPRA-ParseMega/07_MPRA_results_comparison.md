# Compare results to the MPRA results from iGluts and iGABAs
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
scobj <- readRDS("analyses/2026_PerturbSeq_Mega_MPRA/scobj.rds")

# Filter to keep high-quality cells
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
lower_threshold= 1000

# Filtering:
scobj_filt <- subset(scobj, subset = nFeature_RNA > lower_threshold & 
percent.mt < 15 & 
gRNA != 'negative' & 
gRNA != 'ambiguous')
# Total cells retained= 618,617

# Normalize and scale:
scobj_filt <- NormalizeData(scobj_filt, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj_filt <- FindVariableFeatures(scobj_filt, selection.method = "vst", nfeatures = 3000)
scobj_filt <- ScaleData(scobj_filt)
```

```{r}
# =========================================================
# 1. Load metadata and expression matrix
# =========================================================

meta <- read.csv("analyses/2026_PerturbSeq_Mega_MPRA/MPRA_metadata.csv", head=T, stringsAsFactors = T)

meta <- meta %>%
  select(
    target, gene1, dist_gene1, gene2, dist_gene2,
    delta.gaba_glut, fdr, direction
  ) %>%
  distinct()

meta$dist_gene1 <- as.numeric(meta$dist_gene1)
meta$dist_gene2 <- as.numeric(meta$dist_gene2)

DefaultAssay(scobj_filt) <- "RNA"
expr_mat <- GetAssayData(scobj_filt, layer = "data")

cell_meta <- scobj_filt@meta.data
cell_meta$cell_id <- colnames(scobj_filt)

target_col <- "target"
ntc_label <- "NTC"

# =========================================================
# 2. Prepare sgRNA and element information
# =========================================================

cell_meta <- cell_meta %>%
  mutate(
    sgrna = gRNA,
    element = .data[[target_col]],
    sgrna_num = str_extract(gRNA, "[0-9]+$"),
    sgrna_num = as.numeric(sgrna_num)
  )

ntc_cells <- cell_meta$cell_id[cell_meta$element == ntc_label] # 21709 NTC cells

# number of cells per sgRNA
sgrna_counts <- cell_meta %>%
  filter(element != ntc_label) %>%
  count(element, sgrna, name = "n_cells_sgrna")

# number of cells per element
element_counts <- cell_meta %>%
  filter(element != ntc_label) %>%
  count(element, name = "n_cells_element")

# =========================================================
# 3. Build element-gene pairs
# =========================================================

pairs1 <- meta %>%
  transmute(
    element = target,
    nearby_gene = gene1,
    distance = dist_gene1,
    gene_rank = "gene1",
    delta.gaba_glut,
    mpra_fdr = fdr,
    mpra_direction = direction
  )

pairs2 <- meta %>%
  transmute(
    element = target,
    nearby_gene = gene2,
    distance = dist_gene2,
    gene_rank = "gene2",
    delta.gaba_glut,
    mpra_fdr = fdr,
    mpra_direction = direction
  )

pairs_long <- bind_rows(pairs1, pairs2) %>%
  filter(!is.na(nearby_gene), nearby_gene != "") %>%
  distinct(element, nearby_gene, .keep_all = TRUE)

pairs_long <- pairs_long %>%
  group_by(element) %>%
  mutate(
    min_dist = min(abs(distance), na.rm = TRUE),
    is_closest = abs(distance) == min_dist
  ) %>%
  ungroup() %>%
  left_join(element_counts, by = "element")

# =========================================================
# 4. sgRNA-level testing: each sgRNA vs NTC
# =========================================================

all_genes <- rownames(expr_mat)

results_sgrna <- list()
k <- 1

for (i in seq_len(nrow(pairs_long))) {

  this_element <- pairs_long$element[i]
  this_gene <- pairs_long$nearby_gene[i]
  this_dist <- pairs_long$distance[i]

  if (!(this_gene %in% all_genes)) next

  sgrnas_this_element <- unique(
    cell_meta$sgrna[cell_meta$element == this_element]
  )

  for (sg in sgrnas_this_element) {

    sg_cells <- cell_meta$cell_id[cell_meta$sgrna == sg]

    if (length(sg_cells) < 15 || length(ntc_cells) < 20) next

    x_sg <- as.numeric(Matrix::as.matrix(
      expr_mat[this_gene, sg_cells, drop = FALSE]
    ))

    x_ntc <- as.numeric(Matrix::as.matrix(
      expr_mat[this_gene, ntc_cells, drop = FALSE]
    ))

    wt <- suppressWarnings(wilcox.test(x_sg, x_ntc))

    results_sgrna[[k]] <- data.frame(
      element = this_element,
      sgrna = sg,
      nearby_gene = this_gene,
      distance = this_dist,
      abs_distance = abs(this_dist),
      gene_rank = pairs_long$gene_rank[i],
      is_closest = pairs_long$is_closest[i],
      n_sg = length(sg_cells),
      n_ntc = length(ntc_cells),
      mean_sg = mean(x_sg),
      mean_ntc = mean(x_ntc),
      logFC_mean = mean(x_sg) - mean(x_ntc),
      median_sg = median(x_sg),
      median_ntc = median(x_ntc),
      delta_median = median(x_sg) - median(x_ntc),
      pct_sg = mean(x_sg > 0),
      pct_ntc = mean(x_ntc > 0),
      p_val = wt$p.value,
      delta.gaba_glut = pairs_long$delta.gaba_glut[i],
      mpra_fdr = pairs_long$mpra_fdr[i],
      mpra_direction = pairs_long$mpra_direction[i],
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}

results_sgrna <- bind_rows(results_sgrna)

results_sgrna <- results_sgrna %>%
  left_join(sgrna_counts, by = c("element", "sgrna"))

results_sgrna$p_adj <- p.adjust(results_sgrna$p_val, method = "BH")
results_sgrna$sig <- results_sgrna$p_adj < 0.05

write.csv(results_sgrna, "analyses/2026_PerturbSeq_Mega_MPRA/CRISPRi_sgRNA_vs_NTC_results.csv", row.names = FALSE)

# =========================================================
# 5. Consistency across sgRNAs targeting the same element
# =========================================================

sgrna_consistency <- results_sgrna %>%
  group_by(element, nearby_gene) %>%
  summarise(
    n_sgrna = n(),
    mean_effect = mean(logFC_mean, na.rm = TRUE),
    median_effect = median(logFC_mean, na.rm = TRUE),
    sd_effect = sd(logFC_mean, na.rm = TRUE),
    mean_abs_effect = mean(abs(logFC_mean), na.rm = TRUE),
    frac_same_sign = max(
      mean(logFC_mean > 0, na.rm = TRUE),
      mean(logFC_mean < 0, na.rm = TRUE)
    ),
    frac_sig = mean(sig, na.rm = TRUE),
    n_sig = sum(sig, na.rm = TRUE),
    closest_gene = first(is_closest),
    min_distance = min(abs_distance, na.rm = TRUE),
    mpra_fdr = first(mpra_fdr),
    delta.gaba_glut = first(delta.gaba_glut),
    .groups = "drop"
  )

write.csv(sgrna_consistency, "analyses/2026_PerturbSeq_Mega_MPRA/CRISPRi_sgRNA_consistency_summary.csv", row.names = FALSE)

# =========================================================
# 6. Pairwise agreement among sgRNAs within each element-gene pair
# =========================================================

pairwise_list <- list()
k <- 1

for (nm in unique(paste(results_sgrna$element, results_sgrna$nearby_gene, sep = "___"))) {

  df_sub <- results_sgrna %>%
    filter(paste(element, nearby_gene, sep = "___") == nm) %>%
    distinct(sgrna, .keep_all = TRUE)

  if (nrow(df_sub) < 2) next

  combs <- combn(seq_len(nrow(df_sub)), 2, simplify = FALSE)

  for (cc in combs) {
    a <- df_sub[cc[1], ]
    b <- df_sub[cc[2], ]

    pairwise_list[[k]] <- data.frame(
      element = a$element,
      nearby_gene = a$nearby_gene,
      sgrna_1 = a$sgrna,
      sgrna_2 = b$sgrna,
      effect_1 = a$logFC_mean,
      effect_2 = b$logFC_mean,
      abs_diff = abs(a$logFC_mean - b$logFC_mean),
      same_sign = sign(a$logFC_mean) == sign(b$logFC_mean),
      both_sig = a$sig & b$sig,
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}

pairwise_agreement <- bind_rows(pairwise_list)

if (nrow(pairwise_agreement) > 0) {
  write.csv(pairwise_agreement, "analyses/2026_PerturbSeq_Mega_MPRA/CRISPRi_sgRNA_pairwise_agreement.csv", row.names = FALSE)
}

# =========================================================
# 7. Element-level test after pooling concordant sgRNAs
#    Pool only sgRNAs with the same direction of effect
# =========================================================

min_sgrna_per_pool <- 2
min_cells_pool <- 20

element_directional_results <- list()
k <- 1

for (nm in unique(paste(results_sgrna$element, results_sgrna$nearby_gene, sep = "___"))) {

  df_sub <- results_sgrna %>%
    filter(paste(element, nearby_gene, sep = "___") == nm) %>%
    distinct(sgrna, .keep_all = TRUE)

  if (nrow(df_sub) == 0) next

  this_element <- df_sub$element[1]
  this_gene <- df_sub$nearby_gene[1]

  if (!(this_gene %in% rownames(expr_mat))) next

  neg_sgrnas <- df_sub$sgrna[df_sub$logFC_mean < 0]
  pos_sgrnas <- df_sub$sgrna[df_sub$logFC_mean > 0]

  run_pooled_test <- function(sgrna_set, direction_label) {

    if (length(sgrna_set) < min_sgrna_per_pool) return(NULL)

    pooled_cells <- cell_meta %>%
      filter(sgrna %in% sgrna_set) %>%
      pull(cell_id) %>%
      unique()

    if (length(pooled_cells) < min_cells_pool || length(ntc_cells) < 20) return(NULL)

    x_pool <- as.numeric(Matrix::as.matrix(
      expr_mat[this_gene, pooled_cells, drop = FALSE]
    ))

    x_ntc <- as.numeric(Matrix::as.matrix(
      expr_mat[this_gene, ntc_cells, drop = FALSE]
    ))

    wt <- suppressWarnings(wilcox.test(x_pool, x_ntc))

    data.frame(
      element = this_element,
      nearby_gene = this_gene,
      direction_pool = direction_label,
      n_sgrna_pooled = length(sgrna_set),
      sgrna_pooled = paste(sort(sgrna_set), collapse = ";"),
      n_cells_pooled = length(pooled_cells),
      n_ntc = length(ntc_cells),
      mean_pool = mean(x_pool),
      mean_ntc = mean(x_ntc),
      logFC_mean = mean(x_pool) - mean(x_ntc),
      median_pool = median(x_pool),
      median_ntc = median(x_ntc),
      delta_median = median(x_pool) - median(x_ntc),
      pct_pool = mean(x_pool > 0),
      pct_ntc = mean(x_ntc > 0),
      p_val = wt$p.value,
      gene_rank = df_sub$gene_rank[1],
      is_closest = df_sub$is_closest[1],
      abs_distance = df_sub$abs_distance[1],
      distance = df_sub$distance[1],
      delta.gaba_glut = df_sub$delta.gaba_glut[1],
      mpra_fdr = df_sub$mpra_fdr[1],
      mpra_direction = df_sub$mpra_direction[1],
      stringsAsFactors = FALSE
    )
  }

  res_neg <- run_pooled_test(neg_sgrnas, "negative")
  res_pos <- run_pooled_test(pos_sgrnas, "positive")

  if (!is.null(res_neg)) {
    element_directional_results[[k]] <- res_neg
    k <- k + 1
  }

  if (!is.null(res_pos)) {
    element_directional_results[[k]] <- res_pos
    k <- k + 1
  }
}

element_directional_results <- bind_rows(element_directional_results)

if (nrow(element_directional_results) > 0) {
  element_directional_results$p_adj <- p.adjust(element_directional_results$p_val, method = "BH")
  element_directional_results$sig <- element_directional_results$p_adj < 0.1
} else {
  element_directional_results <- data.frame()
}

write.csv(
  element_directional_results,
  "analyses/2026_PerturbSeq_Mega_MPRA/CRISPRi_element_directional_pool_vs_NTC.csv",
  row.names = FALSE
)

# Compact summary per element-gene pair
if (nrow(element_directional_results) > 0) {
  element_directional_summary <- element_directional_results %>%
    group_by(element, nearby_gene) %>%
    summarise(
      n_direction_pools_tested = n(),
      any_sig = any(sig, na.rm = TRUE),
      best_p_adj = min(p_adj, na.rm = TRUE),
      best_abs_effect = max(abs(logFC_mean), na.rm = TRUE),
      best_direction = direction_pool[which.min(p_adj)],
      max_sgrna_pooled = max(n_sgrna_pooled, na.rm = TRUE),
      max_cells_pooled = max(n_cells_pooled, na.rm = TRUE),
      gene_rank = first(gene_rank),
      is_closest = first(is_closest),
      distance = first(distance),
      abs_distance = first(abs_distance),
      .groups = "drop"
    )
} else {
  element_directional_summary <- data.frame()
}

write.csv(
  element_directional_summary,
  "analyses/2026_PerturbSeq_Mega_MPRA/CRISPRi_element_directional_pool_summary.csv",
  row.names = FALSE
)

# =========================================================
# 8. Interesting statistical summaries
# =========================================================

# A. Does distance relate to sgRNA-level effect size?
cor_all <- cor.test(
  results_sgrna$abs_distance,
  abs(results_sgrna$logFC_mean),
  method = "spearman"
)

print(cor_all) # rho= -0.1133002, p-value < 2.2e-16

# B. Does closest gene show bigger sgRNA-level effects?
closest_test_df <- results_sgrna %>%
  group_by(element, nearby_gene) %>%
  summarise(
    mean_abs_effect = mean(abs(logFC_mean), na.rm = TRUE),
    is_closest = first(is_closest),
    .groups = "drop"
  ) %>%
  group_by(element) %>%
  filter(n() == 2) %>%
  ungroup()

if (nrow(closest_test_df) > 0) {
  closest_wide <- closest_test_df %>%
    arrange(element, desc(is_closest)) %>%
    group_by(element) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(element, rank, mean_abs_effect) %>%
    pivot_wider(names_from = rank, values_from = mean_abs_effect)

  if (all(c("1", "2") %in% colnames(closest_wide))) {
    closest_test <- wilcox.test(closest_wide$`1`, closest_wide$`2`, paired = TRUE)
    print(closest_test)
  }
}

# V = 77427, p-value < 2.2e-16

# C. Do elements with stronger effects have better sgRNA agreement?
agreement_cor <- cor.test(
  sgrna_consistency$mean_abs_effect,
  sgrna_consistency$frac_same_sign,
  method = "spearman"
)

print(agreement_cor) # rho= 0.03424306 , p= 0.2793

# D. Does distance relate to pooled directional element-level effect?
if (nrow(element_directional_results) > 2) {
  cor_directional_pool <- cor.test(
    element_directional_results$abs_distance,
    abs(element_directional_results$logFC_mean),
    method = "spearman"
  )
  print(cor_directional_pool)
}

# rho= -0.006793336, p= 0.7923

# =========================================================
# 9. Plots
# =========================================================

# Plot 1: stronger effects and sgRNA agreement
p_consistency <- ggplot(sgrna_consistency, aes(x = mean_abs_effect, y = frac_same_sign)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Mean absolute sgRNA effect",
    y = "Fraction of sgRNAs with same sign",
    title = "Do stronger element-gene effects show better sgRNA agreement?"
  )

ggsave("analyses/2026_PerturbSeq_Mega_MPRA/sgRNA_consistency_vs_effect.png", p_consistency, width = 7, height = 5, dpi = 300)

# Plot 2: distance vs sgRNA-level effect size
p_distance <- ggplot(results_sgrna, aes(x = abs_distance, y = abs(logFC_mean), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(aes(color = sig), method = "lm", se = TRUE) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Absolute distance from element to gene (bp, log10 scale)",
    y = "Absolute sgRNA effect size",
    color = "FDR < 0.05",
    title = "Distance vs sgRNA-level effect size"
  )

ggsave("analyses/2026_PerturbSeq_Mega_MPRA/sgRNA_distance_vs_effect.png", p_distance, width = 7, height = 5, dpi = 300)

# Plot 3: closest vs non-closest gene summary

sgrna_consistency <- sgrna_consistency %>%
  mutate(closest_gene = ifelse(closest_gene, "Closest", "Not closest"))

p_closest <- ggplot(sgrna_consistency, aes(x = closest_gene, y = mean_abs_effect)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  ggsignif::geom_signif(
    comparisons = list(c("Closest", "Not closest")),
    test = "wilcox.test"
  ) +
  theme_bw() +
  labs(
    x = "Closest gene",
    y = "Mean absolute sgRNA effect",
    title = "Are closest genes affected more strongly?"
  )


ggsave("analyses/2026_PerturbSeq_Mega_MPRA/closest_gene_vs_effect_sgRNA_summary.png", p_closest, width = 6, height = 5, dpi = 300)

# Heatmap:

heat_df_cat <- results_df %>%
  mutate(
    change_cat = case_when(
      p_adj < 0.05 & logFC_mean > 0 ~ "Increased",
      p_adj < 0.05 & logFC_mean < 0 ~ "Decreased",
      TRUE ~ "No change"
    ),
    gene_rank = factor(gene_rank, levels = c("gene1", "gene2")),
    change_cat = factor(change_cat, levels = c("Decreased", "No change", "Increased"))
  )

element_order <- heat_df_cat %>%
  group_by(target) %>%
  summarise(best_effect = max(abs(logFC_mean), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(best_effect)) %>%
  pull(target)

heat_df_cat$target <- factor(heat_df_cat$target, levels = rev(element_order))

p_heat_cat <- ggplot(heat_df_cat, aes(x = gene_rank, y = target, fill = change_cat)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Decreased" = "blue",
      "No change" = "gray80",
      "Increased" = "red"
    )
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "Element",
    fill = "",
    title = "Nearby-gene response across all 540 tested elements"
  ) +
  theme(
    axis.text.y = element_text(size = 4),
    panel.grid = element_blank()
  )

p_heat_cat


# Compare to MPRA results

df_plot <- results_df %>%
  filter(!is.na(delta.gaba_glut), !is.na(logFC_mean)) %>%
  mutate(
    crispri_sig = p_adj < 0.05,
    mpra_sig = mpra_fdr < 0.05,
    category = case_when(
      crispri_sig & mpra_sig ~ "Both",
      crispri_sig & !mpra_sig ~ "CRISPRi only",
      !crispri_sig & mpra_sig ~ "MPRA only",
      TRUE ~ "Neither"
    )
  )

p_scatter <- ggplot(df_plot, aes(x = delta.gaba_glut, y = logFC_mean)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(color = category), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(
    x = "MPRA effect (delta.gaba_glut)",
    y = "CRISPRi effect (Target vs NTC)",
    color = "",
    title = "MPRA vs CRISPRi effects"
  )

p_scatter

#
# =========================================================
# 10. Plot for specific elements
# =========================================================

# -----------------------------
# settings
# -----------------------------
this_element <- "chr2_165203510_165203780"
genes_to_plot <- c("SCN2A", "SCN3A")
target_col <- "target"
ntc_label <- "NTC"

cell_meta <- scobj_filt@meta.data
cell_meta$cell_id <- colnames(scobj_filt)

# -----------------------------
# choose cells
# -----------------------------
plot_cells <- cell_meta %>%
  filter(.data[[target_col]] %in% c(this_element, ntc_label)) %>%
  mutate(group = ifelse(.data[[target_col]] == this_element, "Targeted", "NTC"))

# -----------------------------
# build plotting data frame
# -----------------------------
plot_list <- list()

for (g in genes_to_plot) {
  if (!(g %in% rownames(expr_mat))) next
  
  x <- as.numeric(Matrix::as.matrix(
    expr_mat[g, plot_cells$cell_id, drop = FALSE]
  ))
  
  plot_list[[g]] <- data.frame(
    cell_id = plot_cells$cell_id,
    group = plot_cells$group,
    gene = g,
    expr = x,
    stringsAsFactors = FALSE
  )
}

plot_df <- bind_rows(plot_list)

plot_df$group <- factor(plot_df$group, levels = c("NTC", "Targeted"))

# -----------------------------
# optional statistics
# -----------------------------
stats_df <- plot_df %>%
  group_by(gene) %>%
  summarise(
    p_val = wilcox.test(expr ~ group)$p.value,
    mean_targeted = mean(expr[group == "Targeted"]),
    mean_ntc = mean(expr[group == "NTC"]),
    logFC_mean = mean_targeted - mean_ntc,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Wilcoxon p = ", signif(p_val, 3),
                   "\nmean diff = ", round(logFC_mean, 3))
  )

label_y <- plot_df %>%
  group_by(gene) %>%
  summarise(y = max(expr, na.rm = TRUE) * 1.05, .groups = "drop")

stats_df <- left_join(stats_df, label_y, by = "gene")

# -----------------------------
# plot
# -----------------------------
plot_list[[g]] <- data.frame(
    gene = g,
    group = plot_cells$group,
    expr = x,
    stringsAsFactors = FALSE
  )

plot_df <- bind_rows(plot_list) %>%
  filter(!is.na(expr)) %>%
  mutate(group = factor(group, levels = c("NTC", "Targeted")))

p_gene_effect <- ggplot(plot_df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("NTC", "Targeted")),
    test = "wilcox.test"
  ) +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  labs(
    x = "",
    y = "Log-normalized expression",
    title = paste0("Effect of targeting ", this_element)
  ) +
  theme(legend.position = "none")

p_gene_effect
```

```{r}

# =========================================================
# 1. Load files
# =========================================================

vitro <- read_csv("analyses/2026_PerturbSeq_Mega_MPRA/Regions_tested_in_vitro.csv", show_col_types = FALSE)

# optional: join MPRA metadata if available
# keep only useful columns from MPRA metadata
mpra_meta <- meta %>%
  select(
    target,
    gene1, dist_gene1,
    gene2, dist_gene2,
    delta.gaba_glut,
    fdr,
    direction
  ) %>%
  distinct()

# =========================================================
# 2. Prepare in vitro subset table
# =========================================================

# bring in distances / MPRA annotations when region matches target
vitro_meta <- vitro %>%
  rename(target = region) %>%
  left_join(mpra_meta, by = "target", suffix = c(".vitro", ".mpra"))

# use gene names from in vitro file as primary source
pairs1 <- vitro_meta %>%
  transmute(
    target,
    nearby_gene = gene1.vitro,
    gene_rank = "gene1",
    distance = dist_gene1,
    delta.gaba_glut,
    mpra_fdr = fdr,
    mpra_direction = direction
  )

pairs2 <- vitro_meta %>%
  transmute(
    target,
    nearby_gene = gene2.vitro,
    gene_rank = "gene2",
    distance = dist_gene2,
    delta.gaba_glut,
    mpra_fdr = fdr,
    mpra_direction = direction
  )

pairs_vitro <- bind_rows(pairs1, pairs2) %>%
  filter(!is.na(nearby_gene), nearby_gene != "") %>%
  distinct(target, nearby_gene, .keep_all = TRUE)

# closest-gene flag if distance is available
pairs_vitro <- pairs_vitro %>%
  group_by(target) %>%
  mutate(
    abs_distance = abs(distance),
    is_closest = if (all(is.na(abs_distance))) NA else abs_distance == min(abs_distance, na.rm = TRUE)
  ) %>%
  ungroup()

# =========================================================
# 3. Expression matrix and metadata
# =========================================================

DefaultAssay(scobj_filt) <- "RNA"
expr_mat <- GetAssayData(scobj_filt, layer = "data")

cell_meta <- scobj_filt@meta.data
cell_meta$cell_id <- colnames(scobj_filt)

target_col <- "target"
ntc_label <- "NTC"

ntc_cells <- rownames(cell_meta)[cell_meta[[target_col]] == ntc_label]
all_genes <- rownames(expr_mat)

# =========================================================
# 4. Test each in vitro region-gene pair vs NTC
# =========================================================

results_vitro <- list()
k <- 1

for (i in seq_len(nrow(pairs_vitro))) {

  this_target <- pairs_vitro$target[i]
  this_gene <- pairs_vitro$nearby_gene[i]

  if (!(this_gene %in% all_genes)) next

  # exact target match first
  target_cells <- rownames(cell_meta)[cell_meta[[target_col]] == this_target]

  # if no cells match exactly, try matching via gRNA prefix
  if (length(target_cells) == 0 && "gRNA" %in% colnames(cell_meta)) {
    target_cells <- rownames(cell_meta)[grepl(paste0("^", this_target, "_"), cell_meta$gRNA)]
  }

  if (length(target_cells) < 20 || length(ntc_cells) < 20) next

  x_target <- as.numeric(Matrix::as.matrix(
    expr_mat[this_gene, target_cells, drop = FALSE]
  ))
  x_ntc <- as.numeric(Matrix::as.matrix(
    expr_mat[this_gene, ntc_cells, drop = FALSE]
  ))

  wt <- suppressWarnings(wilcox.test(x_target, x_ntc))

  results_vitro[[k]] <- data.frame(
    target = this_target,
    nearby_gene = this_gene,
    gene_rank = pairs_vitro$gene_rank[i],
    distance = pairs_vitro$distance[i],
    abs_distance = pairs_vitro$abs_distance[i],
    is_closest = pairs_vitro$is_closest[i],
    n_target = length(target_cells),
    n_ntc = length(ntc_cells),
    mean_target = mean(x_target),
    mean_ntc = mean(x_ntc),
    logFC_mean = mean(x_target) - mean(x_ntc),
    median_target = median(x_target),
    median_ntc = median(x_ntc),
    delta_median = median(x_target) - median(x_ntc),
    pct_target = mean(x_target > 0),
    pct_ntc = mean(x_ntc > 0),
    p_val = wt$p.value,
    delta.gaba_glut = pairs_vitro$delta.gaba_glut[i],
    mpra_fdr = pairs_vitro$mpra_fdr[i],
    mpra_direction = pairs_vitro$mpra_direction[i],
    stringsAsFactors = FALSE
  )

  k <- k + 1
}

results_vitro <- bind_rows(results_vitro)
results_vitro$p_adj <- p.adjust(results_vitro$p_val, method = "BH")
results_vitro$sig <- results_vitro$p_adj < 0.1

write.csv(results_vitro, "analyses/2026_PerturbSeq_Mega_MPRA/In_vitro_regions_CRISPRi_results.csv", row.names = FALSE)

# =========================================================
# 5. Summary statistics
# =========================================================

summary_vitro <- results_vitro %>%
  group_by(target) %>%
  summarise(
    n_genes_tested = n(),
    n_sig = sum(sig, na.rm = TRUE),
    any_sig = any(sig, na.rm = TRUE),
    best_p_adj = min(p_adj, na.rm = TRUE),
    best_abs_effect = max(abs(logFC_mean), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_vitro, "analyses/2026_PerturbSeq_Mega_MPRA/In_vitro_regions_summary.csv", row.names = FALSE)

cat("\nNumber of tested region-gene pairs:", nrow(results_vitro), "\n")
cat("Number significant at BH FDR < 0.1:", sum(results_vitro$sig, na.rm = TRUE), "\n")
cat("Number of regions with at least one significant nearby gene:",
    sum(summary_vitro$any_sig, na.rm = TRUE), "\n")

# =========================================================
# 6. Statistical tests for this subset
# =========================================================

# A. distance vs effect size
if (sum(!is.na(results_vitro$abs_distance)) > 3) {
  cor_dist <- cor.test(
    results_vitro$abs_distance,
    abs(results_vitro$logFC_mean),
    method = "spearman"
  )
  print(cor_dist)
}

# B. closest vs not-closest
plot_df_closest <- results_vitro %>%
  filter(!is.na(is_closest), !is.na(logFC_mean)) %>%
  mutate(
    closest_gene = ifelse(is_closest, "Closest", "Not closest"),
    closest_gene = factor(closest_gene, levels = c("Closest", "Not closest"))
  )

if (nrow(plot_df_closest) > 1 && n_distinct(plot_df_closest$closest_gene) == 2) {
  wt_closest <- wilcox.test(abs(logFC_mean) ~ closest_gene, data = plot_df_closest)
  print(wt_closest)
}

# C. MPRA vs CRISPRi
plot_df_mpra <- results_vitro %>%
  filter(!is.na(delta.gaba_glut), !is.na(logFC_mean))

if (nrow(plot_df_mpra) > 3) {
  cor_mpra <- cor.test(
    plot_df_mpra$delta.gaba_glut,
    plot_df_mpra$logFC_mean,
    method = "spearman"
  )
  print(cor_mpra)
}

# =========================================================
# 7. Visualizations
# =========================================================

# -------------------------
# Plot 1. Discrete heatmap
# -------------------------
heat_df <- results_vitro %>%
  mutate(
    change_cat = case_when(
      p_adj < 0.05 & logFC_mean > 0 ~ "Increased",
      p_adj < 0.05 & logFC_mean < 0 ~ "Decreased",
      TRUE ~ "No change"
    ),
    gene_rank = factor(gene_rank, levels = c("gene1", "gene2")),
    change_cat = factor(change_cat, levels = c("Decreased", "No change", "Increased"))
  )

region_order <- heat_df %>%
  group_by(target) %>%
  summarise(best_abs_effect = max(abs(logFC_mean), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(best_abs_effect)) %>%
  pull(target)

heat_df$target <- factor(heat_df$target, levels = rev(region_order))

p_heat <- ggplot(heat_df, aes(x = gene_rank, y = target, fill = change_cat)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Decreased" = "blue",
      "No change" = "gray80",
      "Increased" = "red"
    )
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "In vitro-tested region",
    fill = "",
    title = "CRISPRi effect of in vitro-tested regions on nearby genes"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 5)
  )

# -------------------------
# Plot 2. Top hits
# -------------------------
top_hits <- results_vitro %>%
  arrange(p_adj, desc(abs(logFC_mean))) %>%
  slice_head(n = 20)

p_top <- ggplot(
  top_hits,
  aes(
    x = reorder(paste(target, nearby_gene, sep = " | "), logFC_mean),
    y = logFC_mean,
    fill = sig
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "",
    y = "Target - NTC expression difference",
    fill = "FDR < 0.05",
    title = "Top CRISPRi effects in the in vitro-tested subset"
  )

ggsave("In_vitro_regions_top_hits.png", p_top, width = 8, height = 6, dpi = 300)

# -------------------------
# Plot 3. Closest vs not-closest
# -------------------------
if (exists("wt_closest")) {
  p_closest <- ggplot(plot_df_closest, aes(x = closest_gene, y = abs(logFC_mean))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    ggsignif::geom_signif(
      comparisons = list(c("Closest", "Not closest")),
      test = "wilcox.test"
    ) +
    theme_bw() +
    labs(
      x = "",
      y = "Absolute CRISPRi effect",
      title = "Are closest genes more affected in the in vitro-tested subset?"
    )

  ggsave("In_vitro_regions_closest_vs_not.png", p_closest, width = 5, height = 4, dpi = 300)
}

# -------------------------
# Plot 4. MPRA vs CRISPRi
# -------------------------
if (nrow(plot_df_mpra) > 3) {

  plot_df_mpra <- plot_df_mpra %>%
    mutate(
      crispri_sig = p_adj < 0.05,
      mpra_sig = mpra_fdr < 0.05,
      category = case_when(
        crispri_sig & mpra_sig ~ "Both",
        crispri_sig & !mpra_sig ~ "CRISPRi only",
        !crispri_sig & mpra_sig ~ "MPRA only",
        TRUE ~ "Neither"
      )
    )

  p_mpra <- ggplot(plot_df_mpra, aes(x = delta.gaba_glut, y = logFC_mean)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(color = category), alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    theme_bw() +
    labs(
      x = "MPRA effect (delta.gaba_glut)",
      y = "CRISPRi effect (Target - NTC)",
      color = "",
      title = "MPRA vs CRISPRi for in vitro-tested regions"
    )

  ggsave("In_vitro_regions_MPRA_vs_CRISPRi.png", p_mpra, width = 6, height = 5, dpi = 300)
}

# -------------------------
# Plot 5. Distance vs effect
# -------------------------
if (exists("cor_dist")) {
  p_dist <- ggplot(results_vitro, aes(x = abs_distance, y = abs(logFC_mean), color = sig)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_x_log10() +
    theme_bw() +
    labs(
      x = "Absolute distance to nearby gene (bp, log10 scale)",
      y = "Absolute CRISPRi effect",
      color = "FDR < 0.05",
      title = "Does distance matter in the in vitro-tested subset?"
    )

  ggsave("In_vitro_regions_distance_vs_effect.png", p_dist, width = 6, height = 5, dpi = 300)
}

# =========================================================
# 8. Optional: per-region boxplots for selected examples
# =========================================================



```