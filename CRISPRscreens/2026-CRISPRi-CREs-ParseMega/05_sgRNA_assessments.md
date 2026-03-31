## ============================================================
## Obtaining transcriptomic signals for each sgRNA and target
## Streamlined + optimized version
## ============================================================

## ------------------------------------------------------------
## STEP 0: Load libraries
## ------------------------------------------------------------
```{r}
suppressPackageStartupMessages({
  library(Seurat)
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
  library(limma)
  library(splines)
  library(cowplot)
})
```


## ------------------------------------------------------------
## STEP 1: Global settings
## ------------------------------------------------------------
```{r}
outdir <- "analyses/2026_PerturbSeq_Mega_CRE"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fdr_cut <- 0.10
promoter_fdr_cut <- 0.1
effect_thr <- 0.10
min_cells_per_guide <- 30
min_cells_wilcox_sgRNA <- 10
min_cells_wilcox_CRE <- 20
min_total_cells_per_group <- 30
min_cells_per_pb <- 20
mixed_strategy <- "split"   # "split" or "merge" for CREs with both up/down guide trends
```

## ------------------------------------------------------------
## STEP 2: Normalize and scale
## ------------------------------------------------------------
```{r}
scobj_filt <- NormalizeData(scobj_filt, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj_filt <- FindVariableFeatures(scobj_filt, selection.method = "vst", nfeatures = 3000)
scobj_filt <- ScaleData(scobj_filt)
```

## ------------------------------------------------------------
## STEP 3: Helper functions
## ------------------------------------------------------------
```{r}
fill_promoter_target_gene <- function(meta) {
  is_prom <- meta$region_type == "promoter"
  is_prom[is.na(is_prom)] <- FALSE

  missing_tg <- is.na(meta$target_gene) | meta$target_gene == ""
  missing_tg[is.na(missing_tg)] <- TRUE

  prom_src <- meta$target_name
  bad_src <- is.na(prom_src) | prom_src == ""
  prom_src[bad_src] <- meta$target[bad_src]

  prom_tg <- sub("_promoter.*$", "", prom_src)
  prom_tg <- sub("_Promoter.*$", "", prom_tg)

  idx_fill <- which(is_prom & missing_tg & !is.na(prom_tg) & prom_tg != "")
  meta$target_gene[idx_fill] <- prom_tg[idx_fill]
  meta
}

classify_cre <- function(target_name) {
  case_when(
    grepl("control",  target_name, ignore.case = TRUE) ~ "control",
    grepl("silencer", target_name, ignore.case = TRUE) ~ "silencer",
    grepl("enhancer", target_name, ignore.case = TRUE) ~ "enhancer",
    grepl("promoter", target_name, ignore.case = TRUE) ~ "promoter",
    TRUE ~ "other"
  )
}

build_target_expr_df <- function(meta, cells, expr_all_targets, keep_region = NULL) {
  cells <- unique(cells)
  cells <- cells[!is.na(cells)]
  cells <- intersect(cells, colnames(expr_all_targets))
  if (!is.null(keep_region)) {
    cells <- cells[meta[cells, "region_type"] %in% keep_region]
  }
  cells <- cells[
    !is.na(meta[cells, "target_gene"]) &
      meta[cells, "target_gene"] != "" &
      meta[cells, "target_gene"] %in% rownames(expr_all_targets)
  ]

  target_gene_expr <- rep(NA_real_, length(cells))
  names(target_gene_expr) <- cells

  tg_by_cell <- meta[cells, "target_gene"]
  for (g in unique(tg_by_cell)) {
    idx_cells <- cells[tg_by_cell == g]
    if (length(idx_cells) == 0) next
    target_gene_expr[idx_cells] <- as.numeric(as.matrix(expr_all_targets[g, idx_cells, drop = FALSE]))
  }

  out <- meta[cells, c("gRNA", "target", "target_name", "region_type", "target_gene"), drop = FALSE]
  out$cell <- rownames(out)
  out$target_gene_expr <- target_gene_expr[out$cell]

  out %>% filter(!is.na(target_gene_expr))
}

make_ntc_baseline <- function(target_genes, expr_all_targets, cells_ntc) {
  bind_rows(lapply(target_genes, function(g) {
    y <- as.numeric(as.matrix(expr_all_targets[g, cells_ntc, drop = FALSE]))
    data.frame(
      target_gene = g,
      n_ntc = length(y),
      mean_ntc = mean(y),
      median_ntc = median(y),
      stringsAsFactors = FALSE
    )
  }))
}

make_gmap <- function(meta) {
  meta %>%
    filter(!is.na(gRNA), gRNA != "") %>%
    count(gRNA, target_name, target, name = "n") %>%
    group_by(gRNA) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gRNA, target_name, target)
}

make_tmap <- function(meta) {
  meta %>%
    filter(!is.na(target), target != "") %>%
    count(target, target_name, name = "n") %>%
    group_by(target) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(target, target_name)
}

plot_effect_panel <- function(df, title, xlab) {
  ggplot(df %>% filter(!is.na(p_adj), !is.na(logFC_mean)),
         aes(logFC_mean, -log10(p_adj))) +
    geom_point(aes(color = fdr_lab, size = abs_logFC), alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
    geom_text_repel(
      data = df %>% filter(fdr_hit),
      aes(label = target_gene),
      size = 3,
      max.overlaps = 50,
      box.padding = 0.3,
      segment.size = 0.3
    ) +
    theme_bw() +
    scale_color_manual(values = c("NS" = "gray70", "FDR<0.1" = "red")) +
    labs(title = title,
         x = xlab,
         y = "-log10(BH-adjusted p)",
         color = "FDR",
         size = "|logFC|")
}

make_directional_pool_def <- function(guide_tbl, effect_thr = 0.10, mixed_strategy = "split") {
  mixed_strategy <- match.arg(mixed_strategy, c("split", "merge"))

  gdir <- guide_tbl %>%
    mutate(
      direction = case_when(
        logFC_mean <= -effect_thr ~ "down",
        logFC_mean >=  effect_thr ~ "up",
        TRUE ~ "neutral"
      )
    )

  cre_direction_summary <- gdir %>%
    group_by(target_gene, target) %>%
    summarise(
      n_guides = n(),
      n_down = sum(direction == "down"),
      n_up = sum(direction == "up"),
      n_neutral = sum(direction == "neutral"),
      .groups = "drop"
    ) %>%
    mutate(
      mode = case_when(
        n_down > 0 & n_up == 0 ~ "down_only",
        n_up > 0 & n_down == 0 ~ "up_only",
        n_up > 0 & n_down > 0 ~ "mixed",
        TRUE ~ "all_neutral"
      )
    )

  pool_def <- gdir %>%
    left_join(cre_direction_summary, by = c("target_gene", "target")) %>%
    filter(mode != "all_neutral") %>%
    mutate(
      pooled_group = case_when(
        mode == "down_only" & direction == "down" ~ "down",
        mode == "up_only"   & direction == "up"   ~ "up",
        mode == "mixed" & mixed_strategy == "split" & direction == "down" ~ "down",
        mode == "mixed" & mixed_strategy == "split" & direction == "up"   ~ "up",
        mode == "mixed" & mixed_strategy == "merge" & direction %in% c("down", "up") ~ "mixed",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(pooled_group))

  list(
    pool_def = pool_def,
    summary = cre_direction_summary
  )
}

log2fc_to_pct <- function(log2fc) {
  100 * (2^log2fc - 1)
}

pct_to_log2fc <- function(pct) {
  log2(1 + pct / 100)
}

pretty_pct <- function(x, digits = 1) {
  paste0(round(x, digits), "%")
}

```
## ------------------------------------------------------------
## STEP 4: Prepare metadata and expression matrices
## ------------------------------------------------------------
```{r}
obj <- scobj_filt
assay_use <- DefaultAssay(obj)
layer_use <- "data"

meta <- obj@meta.data
meta$cell <- rownames(meta)
meta <- fill_promoter_target_gene(meta)

all_cells <- colnames(obj)

cells_ntc <- unique(meta$cell[meta$region_type == "NTC"])
cells_prom <- unique(meta$cell[meta$region_type == "promoter"])
cells_cre <- unique(meta$cell[meta$region_type == "CRE"])

cells_ntc <- intersect(cells_ntc, all_cells)
cells_prom <- intersect(cells_prom, all_cells)
cells_cre <- intersect(cells_cre, all_cells)

target_genes <- unique(meta$target_gene[meta$region_type %in% c("CRE", "promoter")])
target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]
target_genes <- intersect(target_genes, rownames(obj[[assay_use]]))

if (length(target_genes) == 0) stop("No target genes found in assay.")

expr_all_targets <- GetAssayData(obj, assay = assay_use, layer = layer_use)[target_genes, , drop = FALSE]

cells_ntc <- intersect(cells_ntc, colnames(expr_all_targets))
cells_prom <- intersect(cells_prom, colnames(expr_all_targets))
cells_cre <- intersect(cells_cre, colnames(expr_all_targets))

message("Cells total: ", nrow(meta))
message("NTC cells: ", length(cells_ntc))
message("Promoter cells: ", length(cells_prom))
message("CRE cells: ", length(cells_cre))
message("Target genes in assay: ", length(target_genes))

ntc_baseline <- make_ntc_baseline(target_genes, expr_all_targets, cells_ntc)
gmap <- make_gmap(meta)
tmap <- make_tmap(meta)
```

## ------------------------------------------------------------
## STEP 5: Promoter sanity check
## ------------------------------------------------------------
```{r}
prom_guides <- unique(meta[cells_prom, "gRNA"])
prom_guides <- prom_guides[!is.na(prom_guides) & prom_guides != ""]

promoter_by_gRNA <- bind_rows(lapply(prom_guides, function(gr) {
  prom_cells_gr <- unique(cells_prom[meta[cells_prom, "gRNA"] == gr])
  if (length(prom_cells_gr) == 0) return(NULL)

  g0 <- unique(meta[prom_cells_gr, "target_gene"])
  g0 <- g0[!is.na(g0) & g0 != ""]
  if (length(g0) == 0) return(NULL)
  g0 <- g0[1]
  if (!(g0 %in% target_genes)) return(NULL)

  x <- as.numeric(as.matrix(expr_all_targets[g0, prom_cells_gr, drop = FALSE]))
  y <- as.numeric(as.matrix(expr_all_targets[g0, cells_ntc, drop = FALSE]))

  data.frame(
    target_gene = g0,
    gRNA = gr,
    n_prom = length(x),
    n_ntc = length(y),
    mean_prom = mean(x),
    mean_ntc = mean(y),
    logFC_mean = log2((mean(x)+1e-3)/(mean(y)+1e-3)),
percent_change = log2fc_to_pct(log2((mean(x)+1e-3)/(mean(y)+1e-3))),
    p = ifelse(length(x) >= min_cells_wilcox_sgRNA && length(y) >= min_cells_wilcox_sgRNA,
               wilcox.test(x, y)$p.value, NA_real_),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

write.csv(promoter_by_gRNA, file.path(outdir, "promoter_by_gRNA.csv"), row.names = FALSE)

p_prom_gRNA <- promoter_by_gRNA %>%
  filter(!is.na(p)) %>%
  mutate(significant = p_adj < promoter_fdr_cut) %>%
  ggplot(aes(x = reorder(gRNA, logFC_mean), y = logFC_mean, fill = significant)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  labs(
    title = "Promoter sanity per sgRNA (vs NTC)",
    x = "promoter sgRNA",
    y = "log2(mean promoter sgRNA / mean NTC)",
    fill = "FDR < 0.1"
  )

  p_prom_gRNA2 <- promoter_by_gRNA %>%
  filter(!is.na(p)) %>%
  mutate(significant = p_adj < promoter_fdr_cut) %>%
  ggplot(aes(x = reorder(gRNA, percent_change), y = percent_change, fill = significant)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  labs(
    title = "Promoter sanity per sgRNA (vs NTC)",
    x = "promoter sgRNA",
    y = "Mean % change vs NTC",
    fill = "FDR < 0.1"
  )

ggsave(file.path(outdir, "Check_promoters_sgRNA.pdf"), p_prom_gRNA, height = 8, width = 8)

ggsave(file.path(outdir, "Check_promoters_sgRNA_per_change.pdf"), p_prom_gRNA2, height = 8, width = 8)


ntc_means <- sapply(target_genes, function(g) {
  y <- as.numeric(as.matrix(expr_all_targets[g, cells_ntc, drop = FALSE]))
  mean(y)
}, USE.NAMES = TRUE)

promoter_pooled <- bind_rows(lapply(target_genes, function(g) {
  prom_cells_g <- unique(cells_prom[meta[cells_prom, "target_gene"] == g])
  x <- if (length(prom_cells_g) > 0) as.numeric(as.matrix(expr_all_targets[g, prom_cells_g, drop = FALSE])) else numeric(0)
  y_mean <- ntc_means[g]

  # compute log2FC in a local variable first
  lf <- if (length(x) > 0) log2((mean(x) + 1e-3) / (y_mean + 1e-3)) else NA_real_
  pct <- if (!is.na(lf)) log2fc_to_pct(lf) else NA_real_

  pval <- if (length(x) >= min_cells_wilcox_sgRNA && length(cells_ntc) >= min_cells_wilcox_sgRNA) {
    wilcox.test(x, as.numeric(as.matrix(expr_all_targets[g, cells_ntc, drop = FALSE])))$p.value
  } else NA_real_

  data.frame(
    target_gene = g,
    n_prom = length(x),
    n_ntc = length(cells_ntc),
    mean_prom = ifelse(length(x) > 0, mean(x), NA_real_),
    mean_ntc = y_mean,
    logFC_mean = lf,
    percent_change = pct,
    p = pval,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

write.csv(promoter_pooled, file.path(outdir, "promoter_pooled.csv"), row.names = FALSE)

p_prom_pooled <- promoter_pooled %>%
  filter(!is.na(p)) %>%
  mutate(significant = factor(ifelse(p_adj < promoter_fdr_cut, "Significant", "NS"),
                              levels = c("NS", "Significant"))) %>%
  ggplot(aes(x = reorder(target_gene, logFC_mean), y = logFC_mean, color = significant)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "Significant" = "red")) +
  labs(
    title = "Promoter sanity pooled (all promoter sgRNAs per gene vs NTC)",
    x = "target_gene",
    y = "log2(mean pooled promoter / mean NTC)",
    color = "BH < 0.1"
  )

  p_prom_pooled2 <- promoter_pooled %>%
  filter(!is.na(p)) %>%
  mutate(significant = factor(ifelse(p_adj < promoter_fdr_cut, "Significant", "NS"),
                              levels = c("NS", "Significant"))) %>%
  ggplot(aes(x = reorder(target_gene, percent_change), y = percent_change, color = significant)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "Significant" = "red")) +
  labs(
    title = "Promoter sanity pooled (all promoter sgRNAs per gene vs NTC)",
    x = "target_gene",
    y = "Mean % change vs NTC",
    color = "BH < 0.1"
  )

ggsave(file.path(outdir, "Check_promoters_pooled.pdf"), p_prom_pooled, height = 8, width = 8)

ggsave(file.path(outdir, "Check_promoters_pooled_per_change.pdf"), p_prom_pooled2, height = 8, width = 8)

```

## ------------------------------------------------------------
## STEP 6: Build per-cell target_gene expression tables
## ------------------------------------------------------------
```{r}
df_cre <- build_target_expr_df(meta, cells_cre, expr_all_targets, keep_region = "CRE")
df_prom <- build_target_expr_df(meta, cells_prom, expr_all_targets, keep_region = "promoter")

# combined table, useful for guide-level diagnostics
df_cp <- bind_rows(df_cre, df_prom)
```

## ------------------------------------------------------------
## STEP 7: sgRNA focus (CRE only)
## ------------------------------------------------------------
```{r}
sgRNA_summary <- df_cre %>%
  group_by(target_gene, gRNA, target, target_name) %>%
  summarise(
    n_g = n(),
    mean_expr = mean(target_gene_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(ntc_baseline, by = "target_gene") %>%
  mutate(logFC_mean = log2((mean_expr + 1e-3) / (mean_ntc + 1e-3)))

sgRNA_p <- df_cre %>%
  group_by(target_gene, gRNA) %>%
  summarise(
    p = {
      g0 <- first(target_gene)
      x <- target_gene_expr
      y <- as.numeric(as.matrix(expr_all_targets[g0, cells_ntc, drop = FALSE]))
      if (length(x) >= min_cells_wilcox_sgRNA && length(y) >= min_cells_wilcox_sgRNA) {
        wilcox.test(x, y)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

sgRNA_results <- sgRNA_summary %>%
  left_join(sgRNA_p, by = c("target_gene", "gRNA")) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    fdr_hit = !is.na(p_adj) & p_adj < fdr_cut,
    fdr_lab = factor(ifelse(fdr_hit, "FDR<0.1", "NS"), levels = c("NS", "FDR<0.1")),
    abs_logFC = abs(logFC_mean),
    percent_change = log2fc_to_pct(logFC_mean),
    abs_percent_change = abs(percent_change),
    direction = case_when(
      logFC_mean <= -effect_thr ~ "down",
      logFC_mean >=  effect_thr ~ "up",
      TRUE ~ "neutral"
    ),
    cre_class = factor(classify_cre(target_name),
                       levels = c("enhancer", "silencer", "control", "other"))
  ) %>%
  arrange(p_adj)

write.csv(sgRNA_results, file.path(outdir, "sgRNA_results_CRE_only.csv"), row.names = FALSE)

df_sg_plot <- sgRNA_results %>%
  filter(!is.na(p_adj), !is.na(logFC_mean)) %>%
  mutate(y = -log10(p_adj))

p1a <- ggplot(df_sg_plot, aes(logFC_mean, y)) +
  geom_point(aes(color = fdr_lab, size = abs_logFC), alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  geom_text_repel(
    data = df_sg_plot %>% filter(fdr_hit),
    aes(label = target_gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    segment.size = 0.3
  ) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "FDR<0.1" = "red")) +
  labs(
    title = "sgRNA-focus: FDR vs effect size",
    x = "log2(mean CRE / mean NTC)",
    y = "-log10(BH-adjusted p)",
    color = "FDR",
    size = "|logFC|"
  )

  p1b <- ggplot(df_sg_plot, aes(percent_change, y)) +
  geom_point(aes(color = fdr_lab, size = abs_percent_change), alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  geom_text_repel(
    data = df_sg_plot %>% filter(fdr_hit),
    aes(label = target_gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    segment.size = 0.3
  ) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "FDR<0.1" = "red")) +
  labs(
    title = "sgRNA-focus: FDR vs effect size",
    x = "Mean % change vs NTC",
    y = "-log10(BH-adjusted p)",
    color = "FDR",
    size = "|logFC|"
  )

p2a <- ggplot(df_sg_plot, aes(logFC_mean, y)) +
  geom_point(color = "gray85", alpha = 0.6, size = 3) +
  geom_point(
    data = df_sg_plot %>% filter(cre_class == "control"),
    color = "red",
    alpha = 0.85,
    size = 3
  ) +
  geom_text_repel(
    data = df_sg_plot %>% filter(cre_class == "control", fdr_hit),
    aes(label = target_name),
    size = 3,
    max.overlaps = 60,
    box.padding = 0.3,
    segment.size = 0.25
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  theme_bw() +
  labs(
    title = "Panel 2: Controls highlighted",
    x = "log2(mean gRNA / mean NTC)",
    y = "-log10(BH-adjusted p)"
  )

  p2b <- ggplot(df_sg_plot, aes(percent_change, y)) +
  geom_point(color = "gray85", alpha = 0.6, size = 3) +
  geom_point(
    data = df_sg_plot %>% filter(cre_class == "control"),
    color = "red",
    alpha = 0.85,
    size = 3
  ) +
  geom_text_repel(
    data = df_sg_plot %>% filter(cre_class == "control", fdr_hit),
    aes(label = target_name),
    size = 3,
    max.overlaps = 60,
    box.padding = 0.3,
    segment.size = 0.25
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  theme_bw() +
  labs(
    title = "Panel 2: Controls highlighted",
    x = "Mean % change vs NTC",
    y = "-log10(BH-adjusted p)"
  )

p3a <- ggplot(df_sg_plot, aes(logFC_mean, y)) +
  geom_point(aes(color = cre_class), alpha = 0.75, size = 3) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  theme_bw() +
  labs(
    title = "Panel 3: Colored by CRE class",
    x = "log2(mean gRNA / mean NTC)",
    y = "-log10(BH-adjusted p)",
    color = "cre_class"
  )

  p3b <- ggplot(df_sg_plot, aes(percent_change, y)) +
  geom_point(aes(color = cre_class), alpha = 0.75, size = 3) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  theme_bw() +
  labs(
    title = "Panel 3: Colored by CRE class",
    x = "Mean % change vs NTC",
    y = "-log10(BH-adjusted p)",
    color = "cre_class"
  )

hits_sg <- sgRNA_results %>%
  filter(fdr_hit, !is.na(logFC_mean)) %>%
  mutate(
    expected_dir = case_when(
      cre_class %in% c("enhancer", "control") ~ "Decrease",
      cre_class == "silencer" ~ "Increase",
      TRUE ~ NA_character_
    ),
    observed_dir = case_when(
      logFC_mean < 0 ~ "Decrease",
      logFC_mean > 0 ~ "Increase",
      TRUE ~ "Zero"
    ),
    match = case_when(
      is.na(expected_dir) ~ NA_character_,
      observed_dir == expected_dir ~ "Match",
      TRUE ~ "Mismatch"
    )
  )

denom_sg <- nrow(hits_sg)

bar_df_sg <- hits_sg %>%
  count(cre_class, match, name = "n") %>%
  group_by(cre_class) %>%
  mutate(
    n_cat = sum(n),
    pct_of_hits = 100 * n_cat / denom_sg,
    frac_within_cat = n / n_cat,
    seg_height = pct_of_hits * frac_within_cat
  ) %>%
  ungroup() %>%
  filter(!is.na(match)) %>%
  mutate(match = factor(match, levels = c("Match", "Mismatch")))

p4 <- ggplot(bar_df_sg, aes(cre_class, seg_height, fill = match)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Match" = "red", "Mismatch" = "gray70")) +
  labs(
    title = "C) FDR<0.1 CRE sgRNAs by category (direction match)",
    subtitle = paste0("Bar height = % of all FDR<0.1 CRE sgRNAs (n=", denom_sg, ")."),
    x = "", y = "% of FDR<0.1 CRE sgRNAs", fill = "Match"
  )

sgRNA_panel1 <- (p1a / p2a / p3a / p4) + plot_layout(heights = c(1, 1, 1, 0.9))
sgRNA_panel2 <- (p1b / p2b / p3b / p4) + plot_layout(heights = c(1, 1, 1, 0.9))

ggsave(file.path(outdir, "Plots_sgRNA_focus_CRE_only.pdf"), sgRNA_panel1, height = 14, width = 10)
ggsave(file.path(outdir, "Plots_sgRNA_focus_CRE_only_perc_change.pdf"), sgRNA_panel2, height = 14, width = 10)

```

## ------------------------------------------------------------
## STEP 8: Direction-aware CRE focus (standard Wilcoxon)
## ------------------------------------------------------------
```{r}
directional <- make_directional_pool_def(
  guide_tbl = sgRNA_results %>% select(target_gene, target, gRNA, target_name, logFC_mean),
  effect_thr = effect_thr,
  mixed_strategy = mixed_strategy
)

cre_direction_summary <- directional$summary
cre_pool_def <- directional$pool_def

write.csv(cre_direction_summary, file.path(outdir, "CRE_direction_summary.csv"), row.names = FALSE)
write.csv(cre_pool_def, file.path(outdir, "CRE_directional_pool_definition.csv"), row.names = FALSE)

df_cre_pool <- df_cre %>%
  inner_join(
    cre_pool_def %>% select(target_gene, target, gRNA, pooled_group),
    by = c("target_gene", "target", "gRNA")
  )

cre_pooled_summary <- df_cre_pool %>%
  group_by(target_gene, target, pooled_group) %>%
  summarise(
    n_cells = n(),
    mean_expr = mean(target_gene_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(ntc_baseline, by = "target_gene") %>%
  mutate(logFC_mean = log2((mean_expr + 1e-3) / (mean_ntc + 1e-3)))

cre_pooled_p <- df_cre_pool %>%
  group_by(target_gene, target, pooled_group) %>%
  summarise(
    p = {
      g0 <- first(target_gene)
      x <- target_gene_expr
      y <- as.numeric(as.matrix(expr_all_targets[g0, cells_ntc, drop = FALSE]))
      if (length(x) >= min_cells_wilcox_CRE && length(y) >= min_cells_wilcox_CRE) {
        wilcox.test(x, y)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

cre_results_directional <- cre_pooled_summary %>%
  left_join(cre_pooled_p, by = c("target_gene", "target", "pooled_group")) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    fdr_hit = !is.na(p_adj) & p_adj < fdr_cut,
    fdr_lab = factor(ifelse(fdr_hit, "FDR<0.1", "NS"), levels = c("NS", "FDR<0.1")),
    abs_logFC = abs(logFC_mean),
    percent_change = log2fc_to_pct(logFC_mean),
    abs_percent_change = abs(percent_change)
  ) %>%
  left_join(tmap, by = "target") %>%
  mutate(
    cre_class = factor(classify_cre(target_name),
                       levels = c("enhancer", "silencer", "control", "other"))
  ) %>%
  arrange(p_adj)

write.csv(cre_results_directional, file.path(outdir, "cre_results_directional_wilcox.csv"), row.names = FALSE)

pA_cre <- cre_results_directional %>%
  filter(!is.na(p_adj), !is.na(logFC_mean)) %>%
  ggplot(aes(logFC_mean, -log10(p_adj))) +
  geom_point(aes(color = fdr_lab, size = abs_logFC), alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  geom_text_repel(
    data = cre_results_directional %>% filter(fdr_hit),
    aes(label = target_gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    segment.size = 0.3
  ) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "FDR<0.1" = "red")) +
  labs(
    title = "A) Directional CRE-focus: FDR vs effect size",
    x = "log2(mean pooled CRE / mean NTC)",
    y = "-log10(BH-adjusted p)",
    color = "FDR",
    size = "|logFC|"
  )

  pAb_cre <- cre_results_directional %>%
  filter(!is.na(p_adj), !is.na(logFC_mean)) %>%
  ggplot(aes(percent_change, -log10(p_adj))) +
  geom_point(aes(color = fdr_lab, size = abs_percent_change), alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  geom_text_repel(
    data = cre_results_directional %>% filter(fdr_hit),
    aes(label = target_gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    segment.size = 0.3
  ) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "gray70", "FDR<0.1" = "red")) +
  labs(
    title = "A) Directional CRE-focus: FDR vs effect size",
    x = "Mean % change vs NTC",
    y = "-log10(BH-adjusted p)",
    color = "FDR",
    size = "|logFC|"
  )

pB_cre <- cre_results_directional %>%
  filter(!is.na(p_adj), !is.na(logFC_mean)) %>%
  ggplot(aes(logFC_mean, -log10(p_adj))) +
  geom_point(aes(color = cre_class, size = abs_logFC), alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = 2) +
  geom_text_repel(
    data = cre_results_directional %>% filter(fdr_hit),
    aes(label = target_gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    segment.size = 0.3
  ) +
  theme_bw() +
  labs(
    title = "B) Directional CRE-focus: CRE class vs effect size",
    x = "log2(mean pooled CRE / mean NTC)",
    y = "-log10(BH-adjusted p)",
    color = "CRE class",,
    size = "|logFC|"
  )

hits_cre <- cre_results_directional %>%
  filter(fdr_hit, !is.na(logFC_mean)) %>%
  mutate(
    expected_dir = case_when(
      cre_class %in% c("enhancer", "control") ~ "Decrease",
      cre_class == "silencer" ~ "Increase",
      TRUE ~ NA_character_
    ),
    observed_dir = case_when(
      logFC_mean < 0 ~ "Decrease",
      logFC_mean > 0 ~ "Increase",
      TRUE ~ "Zero"
    ),
    match = case_when(
      is.na(expected_dir) ~ NA_character_,
      observed_dir == expected_dir ~ "Match",
      TRUE ~ "Mismatch"
    )
  )

denom_cre <- nrow(hits_cre)

bar_df_cre <- hits_cre %>%
  count(cre_class, match, name = "n") %>%
  group_by(cre_class) %>%
  mutate(
    pct_of_hits = 100 * sum(n) / denom_cre,
    frac_within_cat = n / sum(n),
    seg_height = pct_of_hits * frac_within_cat
  ) %>%
  ungroup() %>%
  filter(!is.na(match)) %>%
  mutate(match = factor(match, levels = c("Match", "Mismatch")))

pC_cre <- ggplot(bar_df_cre, aes(cre_class, seg_height, fill = match)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Match" = "red", "Mismatch" = "gray70")) +
  labs(
    title = "C) FDR<0.1 directional CRE targets by category",
    subtitle = paste0("Bar height = % of all FDR<0.1 directional CRE targets (n=", denom_cre, ")."),
    x = "", y = "% of FDR<0.1 directional CRE targets", fill = "Match"
  )

cre_panel <- (pA_cre / pB_cre / pC_cre) + plot_layout(heights = c(1, 1, 0.9))
ggsave(file.path(outdir, "Plots_CREs_focus_directional_wilcox.pdf"), cre_panel, height = 12, width = 10)
```

## ------------------------------------------------------------
## STEP 9: Concordance summary between guides targeting the same CRE
## ------------------------------------------------------------
```{r}
guide_level <- sgRNA_results %>%
  select(target_gene, target, gRNA, n_g, mean_expr, logFC_mean, percent_change, p, p_adj, fdr_hit) %>%
  mutate(
    direction_all = case_when(
      logFC_mean >=  effect_thr ~ "Up",
      logFC_mean <= -effect_thr ~ "Down",
      TRUE ~ "Neutral"
    )
  )

cre_dir_all <- guide_level %>%
  group_by(target_gene, target) %>%
  summarise(
    n_guides = n(),
    n_up = sum(direction_all == "Up"),
    n_down = sum(direction_all == "Down"),
    n_neutral = sum(direction_all == "Neutral"),
    pct_up = 100 * n_up / n_guides,
    pct_down = 100 * n_down / n_guides,
    pct_neutral = 100 * n_neutral / n_guides,
    median_logFC_all = median(logFC_mean, na.rm = TRUE),
    .groups = "drop"
  )

sort_by <- "pct_down"

cre_dir_all_sorted <- cre_dir_all %>%
  mutate(label = paste(target_gene, target, sep = " | ")) %>%
  arrange(desc(.data[[sort_by]])) %>%
  mutate(label = factor(label, levels = unique(label)))

cre_long_all <- cre_dir_all_sorted %>%
  pivot_longer(
    cols = c(pct_up, pct_down, pct_neutral),
    names_to = "direction",
    values_to = "pct"
  ) %>%
  mutate(
    direction = recode(direction,
                       pct_up = "Up",
                       pct_down = "Down",
                       pct_neutral = "Neutral"),
    direction = factor(direction, levels = c("Up", "Down", "Neutral"))
  )

p_stack_all <- ggplot(cre_long_all, aes(x = label, y = pct, fill = direction)) +
  geom_col(width = 0.9) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "Neutral" = "gray70")) +
  labs(
    title = paste0("Per-CRE sgRNA direction composition (sorted by ", sort_by, ")"),
    subtitle = paste0("Neutral if |logFC| < ", effect_thr),
    x = "",
    y = "% of sgRNAs",
    fill = "Direction"
  ) +
  theme(axis.text.y = element_text(size = 4))

ggsave(file.path(outdir, "CRE_direction_composition_ALL_sorted.pdf"), p_stack_all, width = 10, height = 40)

write.csv(guide_level, file.path(outdir, "concordance_guide_level_with_direction.csv"), row.names = FALSE)
write.csv(cre_dir_all, file.path(outdir, "concordance_CRE_direction_summary.csv"), row.names = FALSE)
```

## ------------------------------------------------------------
## STEP 10: Evaluate if closer sgRNAs have same response (sgRNA-focus)
## ------------------------------------------------------------
```{r}
## ============================================================
## STEP X: Are nearby sgRNAs within a CRE more concordant?
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# ----------------------------
# User inputs
# ----------------------------
# Replace with your guide annotation spreadsheet containing sgRNA sequences
guide_annot_file <- "path/to/your_sgRNA_annotation.csv"

# Use the same effect threshold as in the rest of your pipeline
effect_thr <- 0.10

# If TRUE, only use uniquely mapped perfect matches
keep_unique_only <- TRUE

# ----------------------------
# 1) Read guide annotation
# ----------------------------
guide_annot <- read_csv(guide_annot_file, show_col_types = FALSE)

stopifnot(all(c("gRNA", "sgRNA_sequence") %in% colnames(guide_annot)))

guide_annot <- guide_annot %>%
  mutate(
    sgRNA_sequence = toupper(gsub("[^ACGT]", "", sgRNA_sequence))
  ) %>%
  filter(!is.na(gRNA), gRNA != "", !is.na(sgRNA_sequence), nchar(sgRNA_sequence) > 0)

# ----------------------------
# 2) Map sgRNA sequences to hg38 by exact match
#    (search both forward and reverse-complement strands)
# ----------------------------
hg38 <- BSgenome.Hsapiens.UCSC.hg38

map_one_guide <- function(seq, guide_name, genome) {
  dna <- DNAString(seq)
  dna_rc <- reverseComplement(dna)

  out <- list()

  for (chr_name in seqnames(genome)) {
    chr_seq <- genome[[chr_name]]

    # forward-strand match
    m_fwd <- matchPattern(dna, chr_seq)
    if (length(m_fwd) > 0) {
      out[[length(out) + 1]] <- data.frame(
        gRNA = guide_name,
        sgRNA_sequence = seq,
        chr = as.character(chr_name),
        start = start(m_fwd),
        end = end(m_fwd),
        strand = "+",
        stringsAsFactors = FALSE
      )
    }

    # reverse-strand match
    m_rev <- matchPattern(dna_rc, chr_seq)
    if (length(m_rev) > 0) {
      out[[length(out) + 1]] <- data.frame(
        gRNA = guide_name,
        sgRNA_sequence = seq,
        chr = as.character(chr_name),
        start = start(m_rev),
        end = end(m_rev),
        strand = "-",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

guide_maps <- bind_rows(lapply(seq_len(nrow(guide_annot)), function(i) {
  map_one_guide(
    seq = guide_annot$sgRNA_sequence[i],
    guide_name = guide_annot$gRNA[i],
    genome = hg38
  )
}))

if (nrow(guide_maps) == 0) {
  stop("No guides mapped to hg38.")
}

guide_maps <- guide_maps %>%
  mutate(midpoint = (start + end) / 2)

# Count matches per guide
guide_map_counts <- guide_maps %>%
  count(gRNA, name = "n_matches")

guide_maps <- guide_maps %>%
  left_join(guide_map_counts, by = "gRNA")

if (keep_unique_only) {
  guide_maps <- guide_maps %>% filter(n_matches == 1)
}

write.csv(
  guide_maps,
  file.path(outdir, "sgRNA_hg38_mapping.csv"),
  row.names = FALSE
)

# ----------------------------
# 3) Merge mapping with sgRNA results
# ----------------------------
sgRNA_results_geo <- sgRNA_results %>%
  left_join(
    guide_annot %>% select(gRNA, sgRNA_sequence),
    by = "gRNA"
  ) %>%
  left_join(
    guide_maps %>% select(gRNA, chr, start, end, midpoint, strand, n_matches),
    by = "gRNA"
  ) %>%
  mutate(
    direction = case_when(
      logFC_mean <= -effect_thr ~ "down",
      logFC_mean >=  effect_thr ~ "up",
      TRUE ~ "neutral"
    )
  )

write.csv(
  sgRNA_results_geo,
  file.path(outdir, "sgRNA_results_with_coordinates.csv"),
  row.names = FALSE
)

# ----------------------------
# 4) Build pairwise sgRNA comparisons within the same CRE
# ----------------------------
pairwise_df <- sgRNA_results_geo %>%
  filter(
    !is.na(target), !is.na(gRNA), !is.na(midpoint), !is.na(chr),
    region_type == "CRE"
  ) %>%
  select(target_gene, target, gRNA, chr, midpoint, direction, logFC_mean, p_adj) %>%
  inner_join(
    .,
    .,
    by = c("target_gene", "target", "chr"),
    suffix = c("_1", "_2")
  ) %>%
  filter(gRNA_1 < gRNA_2) %>%
  mutate(
    distance_bp = abs(midpoint_1 - midpoint_2),
    same_direction = direction_1 == direction_2,
    both_non_neutral = direction_1 != "neutral" & direction_2 != "neutral",
    same_non_neutral_direction = both_non_neutral & (direction_1 == direction_2),
    opposite_non_neutral_direction = both_non_neutral & (direction_1 != direction_2)
  )

write.csv(
  pairwise_df,
  file.path(outdir, "sgRNA_pairwise_distance_direction.csv"),
  row.names = FALSE
)

# ----------------------------
# 5) Summary statistics
# ----------------------------
pairwise_summary <- pairwise_df %>%
  summarise(
    n_pairs = n(),
    median_distance_bp = median(distance_bp, na.rm = TRUE),
    pct_same_direction = mean(same_direction, na.rm = TRUE) * 100,
    pct_both_non_neutral = mean(both_non_neutral, na.rm = TRUE) * 100,
    pct_same_non_neutral_direction = mean(same_non_neutral_direction, na.rm = TRUE) * 100
  )

print(pairwise_summary)

write.csv(
  pairwise_summary,
  file.path(outdir, "sgRNA_pairwise_distance_direction_summary.csv"),
  row.names = FALSE
)

# ----------------------------
# 6) Statistical test:
#    are closer pairs more likely to have same direction?
# ----------------------------
pairwise_test <- pairwise_df %>%
  filter(!is.na(distance_bp), !is.na(same_direction)) %>%
  mutate(
    log10_distance = log10(distance_bp + 1),
    same_direction_num = as.integer(same_direction)
  )

m_same <- glm(
  same_direction_num ~ log10_distance,
  data = pairwise_test,
  family = binomial()
)

summary(m_same)

capture.output(
  summary(m_same),
  file = file.path(outdir, "glm_same_direction_vs_distance.txt")
)

# Optional stricter test: exclude neutral-neutral / neutral-any pairs
pairwise_test_non_neutral <- pairwise_df %>%
  filter(both_non_neutral) %>%
  mutate(
    log10_distance = log10(distance_bp + 1),
    same_non_neutral_num = as.integer(same_non_neutral_direction)
  )

if (nrow(pairwise_test_non_neutral) > 10) {
  m_same_non_neutral <- glm(
    same_non_neutral_num ~ log10_distance,
    data = pairwise_test_non_neutral,
    family = binomial()
  )

  capture.output(
    summary(m_same_non_neutral),
    file = file.path(outdir, "glm_same_non_neutral_direction_vs_distance.txt")
  )
}

# ----------------------------
# 7) Visualizations
# ----------------------------

# Panel A: scatter of pairwise distance vs same-direction
p_dist_scatter <- ggplot(pairwise_test, aes(x = log10_distance, y = same_direction_num)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.15, size = 1) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  theme_bw() +
  labs(
    title = "Are nearby sgRNAs more likely to show the same direction?",
    x = "log10(pairwise distance in bp + 1)",
    y = "Same direction (0/1)"
  )

# Panel B: boxplot of distances by concordance
p_dist_box <- ggplot(pairwise_df, aes(x = same_direction, y = distance_bp)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Pairwise genomic distance by direction concordance",
    x = "Same direction",
    y = "Distance (bp, log10 scale)"
  )

# Panel C: stricter non-neutral comparison
p_non_neutral_box <- ggplot(
  pairwise_df %>% filter(both_non_neutral),
  aes(x = same_non_neutral_direction, y = distance_bp)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Distance among non-neutral sgRNA pairs",
    x = "Same non-neutral direction",
    y = "Distance (bp, log10 scale)"
  )

# Panel D: per-CRE fraction concordant
cre_pairwise_summary <- pairwise_df %>%
  group_by(target_gene, target) %>%
  summarise(
    n_pairs = n(),
    median_distance_bp = median(distance_bp, na.rm = TRUE),
    pct_same_direction = mean(same_direction, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_same_direction))

write.csv(
  cre_pairwise_summary,
  file.path(outdir, "CRE_pairwise_concordance_summary.csv"),
  row.names = FALSE
)

p_cre_summary <- ggplot(cre_pairwise_summary, aes(x = median_distance_bp, y = pct_same_direction)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Per-CRE concordance vs median sgRNA spacing",
    x = "Median pairwise distance within CRE (bp, log10 scale)",
    y = "% sgRNA pairs with same direction"
  )

pairwise_panel <- (p_dist_scatter / p_dist_box / p_non_neutral_box / p_cre_summary) +
  plot_layout(heights = c(1, 1, 1, 1))

ggsave(
  file.path(outdir, "sgRNA_distance_vs_direction_panels.pdf"),
  pairwise_panel,
  height = 14,
  width = 10
)

# ----------------------------
# 8) Optional: inspect individual CREs
# ----------------------------
# Plot sgRNAs within a single CRE ordered by genomic position
plot_one_cre <- function(cre_name, tbl = sgRNA_results_geo) {
  df_one <- tbl %>%
    filter(target == cre_name, !is.na(midpoint)) %>%
    arrange(midpoint)

  if (nrow(df_one) == 0) return(NULL)

  ggplot(df_one, aes(x = midpoint, y = logFC_mean, color = direction, label = gRNA)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    theme_bw() +
    labs(
      title = paste0("sgRNAs across CRE: ", cre_name),
      x = "Genomic midpoint (hg38)",
      y = "sgRNA logFC"
    )
}

# Example:
# p_example_cre <- plot_one_cre("APP_enhancer_intron_520")
# print(p_example_cre)
# ggsave(file.path(outdir, "Example_CRE_layout_APP_enhancer_intron_520.pdf"),
#        p_example_cre, width = 10, height = 5)
```

## ------------------------------------------------------------
## STEP 11: Evaluate if closer sgRNAs have same response (CRE-focus)
## ------------------------------------------------------------
```{r}
## ============================================================
## STEP X: Are nearby CREs more likely to show the same trend?
## Uses cre_results_directional_wilcox.csv (or object in memory)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(patchwork)
})

# ----------------------------
# Settings
# ----------------------------
effect_thr <- 0.10
fdr_cut <- 0.10

# Restrict pairwise comparisons:
#   TRUE  = only compare CREs assigned to the same target_gene
#   FALSE = compare all CREs on same chromosome
restrict_same_target_gene <- TRUE

# Whether to analyze all CREs or also a significant-only subset
analyze_significant_subset <- TRUE

# ----------------------------
# 1) Load CRE results
# ----------------------------
# If already in memory, use that. Otherwise read from disk.
if (!exists("cre_results_directional")) {
  cre_results_directional <- read_csv(
    file.path(outdir, "cre_results_directional_wilcox.csv"),
    show_col_types = FALSE
  )
}

# ----------------------------
# 2) Parse genomic coordinates
# ----------------------------
# Expects format like: chr21:26168392-26169982
cre_geo <- cre_results_directional %>%
  filter(!is.na(coordinates), coordinates != "") %>%
  mutate(
    chr = str_extract(coordinates, "^chr[^:]+"),
    coord_start = str_extract(coordinates, "(?<=:)[0-9]+"),
    coord_end   = str_extract(coordinates, "(?<=-)[0-9]+"),
    start = as.numeric(coord_start),
    end   = as.numeric(coord_end),
    midpoint = (start + end) / 2,
    width = end - start + 1,
    direction = case_when(
      logFC_mean <= -effect_thr ~ "down",
      logFC_mean >=  effect_thr ~ "up",
      TRUE ~ "neutral"
    ),
    direction = factor(direction, levels = c("down", "neutral", "up")),
    row_id = row_number()
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(midpoint))

write.csv(
  cre_geo,
  file.path(outdir, "cre_results_directional_with_coordinates_parsed.csv"),
  row.names = FALSE
)

# ----------------------------
# 3) Pairwise CRE comparisons
# ----------------------------
pair_by <- c("chr")
if (restrict_same_target_gene) {
  pair_by <- c(pair_by, "target_gene")
}

cre_pairs <- cre_geo %>%
  select(
    row_id, target_gene, target, target_name, pooled_group,
    chr, start, end, midpoint, width,
    logFC_mean, p_adj, fdr_hit, direction, cre_class
  ) %>%
  inner_join(
    .,
    .,
    by = pair_by,
    suffix = c("_1", "_2")
  ) %>%
  filter(row_id_1 < row_id_2) %>%
  mutate(
    distance_bp = abs(midpoint_1 - midpoint_2),
    same_direction = direction_1 == direction_2,
    both_non_neutral = direction_1 != "neutral" & direction_2 != "neutral",
    same_non_neutral_direction = both_non_neutral & (direction_1 == direction_2),
    opposite_non_neutral_direction = both_non_neutral & (direction_1 != direction_2),
    both_sig = fdr_hit_1 & fdr_hit_2
  )

write.csv(
  cre_pairs,
  file.path(outdir, "CRE_pairwise_distance_direction.csv"),
  row.names = FALSE
)

# ----------------------------
# 4) Summary tables
# ----------------------------
pairwise_summary_all <- cre_pairs %>%
  summarise(
    n_pairs = n(),
    median_distance_bp = median(distance_bp, na.rm = TRUE),
    pct_same_direction = mean(same_direction, na.rm = TRUE) * 100,
    pct_both_non_neutral = mean(both_non_neutral, na.rm = TRUE) * 100,
    pct_same_non_neutral_direction = mean(same_non_neutral_direction, na.rm = TRUE) * 100
  )

print(pairwise_summary_all)

write.csv(
  pairwise_summary_all,
  file.path(outdir, "CRE_pairwise_distance_direction_summary_all.csv"),
  row.names = FALSE
)

if (analyze_significant_subset) {
  pairwise_summary_sig <- cre_pairs %>%
    filter(both_sig) %>%
    summarise(
      n_pairs = n(),
      median_distance_bp = median(distance_bp, na.rm = TRUE),
      pct_same_direction = mean(same_direction, na.rm = TRUE) * 100,
      pct_both_non_neutral = mean(both_non_neutral, na.rm = TRUE) * 100,
      pct_same_non_neutral_direction = mean(same_non_neutral_direction, na.rm = TRUE) * 100
    )

  print(pairwise_summary_sig)

  write.csv(
    pairwise_summary_sig,
    file.path(outdir, "CRE_pairwise_distance_direction_summary_sig.csv"),
    row.names = FALSE
  )
}

# ----------------------------
# 5) Statistical tests:
#    does smaller distance predict same direction?
# ----------------------------
pairwise_test_all <- cre_pairs %>%
  filter(!is.na(distance_bp), !is.na(same_direction)) %>%
  mutate(
    log10_distance = log10(distance_bp + 1),
    same_direction_num = as.integer(same_direction)
  )

m_same_all <- glm(
  same_direction_num ~ log10_distance,
  data = pairwise_test_all,
  family = binomial()
)

capture.output(
  summary(m_same_all),
  file = file.path(outdir, "glm_CRE_same_direction_vs_distance_all.txt")
)

# stricter version: only non-neutral pairs
pairwise_test_non_neutral <- cre_pairs %>%
  filter(both_non_neutral) %>%
  mutate(
    log10_distance = log10(distance_bp + 1),
    same_non_neutral_num = as.integer(same_non_neutral_direction)
  )

if (nrow(pairwise_test_non_neutral) > 10) {
  m_same_non_neutral <- glm(
    same_non_neutral_num ~ log10_distance,
    data = pairwise_test_non_neutral,
    family = binomial()
  )

  capture.output(
    summary(m_same_non_neutral),
    file = file.path(outdir, "glm_CRE_same_non_neutral_direction_vs_distance.txt")
  )
}

# significant-only subset
if (analyze_significant_subset) {
  pairwise_test_sig <- cre_pairs %>%
    filter(both_sig) %>%
    mutate(
      log10_distance = log10(distance_bp + 1),
      same_direction_num = as.integer(same_direction)
    )

  if (nrow(pairwise_test_sig) > 10) {
    m_same_sig <- glm(
      same_direction_num ~ log10_distance,
      data = pairwise_test_sig,
      family = binomial()
    )

    capture.output(
      summary(m_same_sig),
      file = file.path(outdir, "glm_CRE_same_direction_vs_distance_sig_only.txt")
    )
  }
}

# ----------------------------
# 6) Visualization panels
# ----------------------------

# Panel A: all pairs, logistic trend
p_cre_dist_scatter <- ggplot(pairwise_test_all, aes(x = log10_distance, y = same_direction_num)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.12, size = 1) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  theme_bw() +
  labs(
    title = "All CRE pairs: are nearby CREs more likely to show the same direction?",
    x = "log10(pairwise distance in bp + 1)",
    y = "Same direction (0/1)"
  )

# Panel B: boxplot of distance by same vs different direction
p_cre_dist_box <- ggplot(cre_pairs, aes(x = same_direction, y = distance_bp)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "All CRE pairs: genomic distance by direction concordance",
    x = "Same direction",
    y = "Distance (bp, log10 scale)"
  )

# Panel C: stricter non-neutral comparison
p_cre_non_neutral_box <- ggplot(
  cre_pairs %>% filter(both_non_neutral),
  aes(x = same_non_neutral_direction, y = distance_bp)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Non-neutral CRE pairs: distance by concordance",
    x = "Same non-neutral direction",
    y = "Distance (bp, log10 scale)"
  )

# Panel D: per-target-gene summary
cre_pairwise_summary <- cre_pairs %>%
  group_by(target_gene) %>%
  summarise(
    n_pairs = n(),
    median_distance_bp = median(distance_bp, na.rm = TRUE),
    pct_same_direction = mean(same_direction, na.rm = TRUE) * 100,
    pct_same_non_neutral = mean(same_non_neutral_direction, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_same_direction))

write.csv(
  cre_pairwise_summary,
  file.path(outdir, "CRE_pairwise_concordance_by_target_gene.csv"),
  row.names = FALSE
)

p_cre_summary <- ggplot(cre_pairwise_summary, aes(x = median_distance_bp, y = pct_same_direction)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Per target_gene: CRE concordance vs median CRE spacing",
    x = "Median pairwise CRE distance (bp, log10 scale)",
    y = "% CRE pairs with same direction"
  )

cre_pair_panel <- (p_cre_dist_scatter / p_cre_dist_box / p_cre_non_neutral_box / p_cre_summary) +
  plot_layout(heights = c(1, 1, 1, 1))

ggsave(
  file.path(outdir, "CRE_distance_vs_direction_panels.pdf"),
  cre_pair_panel,
  height = 14,
  width = 10
)

# ----------------------------
# 7) Optional significant-only visualizations
# ----------------------------
if (analyze_significant_subset) {

  sig_pairs <- cre_pairs %>% filter(both_sig)

  p_sig_scatter <- ggplot(
    sig_pairs %>% mutate(log10_distance = log10(distance_bp + 1),
                         same_direction_num = as.integer(same_direction)),
    aes(x = log10_distance, y = same_direction_num)
  ) +
    geom_jitter(height = 0.05, width = 0, alpha = 0.12, size = 1) +
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
    theme_bw() +
    labs(
      title = "Significant CRE pairs only: same direction vs distance",
      x = "log10(pairwise distance in bp + 1)",
      y = "Same direction (0/1)"
    )

  p_sig_box <- ggplot(sig_pairs, aes(x = same_direction, y = distance_bp)) +
    geom_boxplot() +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = "Significant CRE pairs only: distance by concordance",
      x = "Same direction",
      y = "Distance (bp, log10 scale)"
    )

  ggsave(
    file.path(outdir, "CRE_distance_vs_direction_sig_only.pdf"),
    p_sig_scatter / p_sig_box,
    height = 8,
    width = 10
  )
}

# ----------------------------
# 8) Optional: inspect local CRE neighborhoods for one target_gene
# ----------------------------
plot_cre_neighborhood <- function(gene_name, tbl = cre_geo) {
  df_one <- tbl %>%
    filter(target_gene == gene_name) %>%
    arrange(chr, midpoint)

  if (nrow(df_one) == 0) return(NULL)

  ggplot(df_one, aes(x = midpoint, y = logFC_mean, color = direction, label = target)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    facet_wrap(~ chr, scales = "free_x") +
    theme_bw() +
    labs(
      title = paste0("CRE neighborhood plot: ", gene_name),
      x = "CRE midpoint",
      y = "Directional CRE logFC"
    )
}

# Example:
# p_example_gene <- plot_cre_neighborhood("APP")
# print(p_example_gene)
# ggsave(file.path(outdir, "CRE_neighborhood_APP.pdf"), p_example_gene, width = 10, height = 6)
```


# ============================
# Simple plot for one CRE vs NTC
# Define CRE of interest by target_name
# ============================
```{r}
# ============================================================
# Plot one CRE in directional groups:
#   1) all sgRNAs vs NTC
#   2) down sgRNAs vs NTC
#   3) up sgRNAs vs NTC
#   4) per-sgRNA distributions
#   5) sgRNA-level effect sizes
# ============================================================

cre_of_interest <- "SYNGAP1_enhancer_3UTR_1068"

# ----------------------------
# 1) All cells for this CRE
# ----------------------------
df_sub_all <- df_cre %>%
  filter(target_name == cre_of_interest)

if (nrow(df_sub_all) == 0) {
  stop("No cells found for CRE: ", cre_of_interest)
}

tg <- unique(df_sub_all$target_gene)
tg <- tg[!is.na(tg) & tg != ""]
if (length(tg) == 0) stop("No target_gene found for CRE: ", cre_of_interest)
tg <- tg[1]

y_ntc <- as.numeric(as.matrix(expr_all_targets[tg, cells_ntc, drop = FALSE]))

# ----------------------------
# 2) Directional pooling info
# ----------------------------
pool_info <- cre_pool_def %>%
  filter(target_name == cre_of_interest)

if (nrow(pool_info) == 0) {
  stop("No directional pooling definition found for CRE: ", cre_of_interest)
}

guides_down <- unique(pool_info$gRNA[pool_info$pooled_group == "down"])
guides_up   <- unique(pool_info$gRNA[pool_info$pooled_group == "up"])

df_sub_down <- df_cre %>% filter(gRNA %in% guides_down)
df_sub_up   <- df_cre %>% filter(gRNA %in% guides_up)

# ----------------------------
# 3) Summary stats
# ----------------------------
lf_all <- log2((mean(df_sub_all$target_gene_expr) + 1e-3) / (mean(y_ntc) + 1e-3))
pct_all <- 100 * (2^lf_all - 1)
p_all <- if (nrow(df_sub_all) >= 10 && length(y_ntc) >= 10) {
  wilcox.test(df_sub_all$target_gene_expr, y_ntc)$p.value
} else {
  NA_real_
}

lf_down <- if (nrow(df_sub_down) > 0) {
  log2((mean(df_sub_down$target_gene_expr) + 1e-3) / (mean(y_ntc) + 1e-3))
} else NA_real_
pct_down <- if (!is.na(lf_down)) 100 * (2^lf_down - 1) else NA_real_
p_down <- if (nrow(df_sub_down) >= 10 && length(y_ntc) >= 10) {
  wilcox.test(df_sub_down$target_gene_expr, y_ntc)$p.value
} else {
  NA_real_
}

lf_up <- if (nrow(df_sub_up) > 0) {
  log2((mean(df_sub_up$target_gene_expr) + 1e-3) / (mean(y_ntc) + 1e-3))
} else NA_real_
pct_up <- if (!is.na(lf_up)) 100 * (2^lf_up - 1) else NA_real_
p_up <- if (nrow(df_sub_up) >= 10 && length(y_ntc) >= 10) {
  wilcox.test(df_sub_up$target_gene_expr, y_ntc)$p.value
} else {
  NA_real_
}

# rows from CRE pooled results for this CRE
cre_row <- cre_results_directional %>%
  filter(target_name == cre_of_interest)

print(cre_row)

# ----------------------------
# 4) Panel A: all sgRNAs vs NTC
# ----------------------------
plot_df_all <- bind_rows(
  data.frame(group = "All sgRNAs", expr = df_sub_all$target_gene_expr, stringsAsFactors = FALSE),
  data.frame(group = "NTC", expr = y_ntc, stringsAsFactors = FALSE)
)

p_cre_all <- ggplot(plot_df_all, aes(x = group, y = expr, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = paste0(cre_of_interest, ": all sgRNAs vs NTC"),
    subtitle = paste0(
      "Target gene: ", tg,
      " | log2FC = ", round(lf_all, 3),
      " | % change = ", round(pct_all, 1), "%",
      " | p = ", signif(p_all, 3),
      " | n CRE cells = ", nrow(df_sub_all)
    ),
    x = "",
    y = paste0(tg, " expression")
  )

# ----------------------------
# 5) Panel B: down sgRNAs vs NTC
# ----------------------------
if (nrow(df_sub_down) > 0) {
  plot_df_down <- bind_rows(
    data.frame(group = "Down sgRNAs", expr = df_sub_down$target_gene_expr, stringsAsFactors = FALSE),
    data.frame(group = "NTC", expr = y_ntc, stringsAsFactors = FALSE)
  )

  p_cre_down <- ggplot(plot_df_down, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = paste0(cre_of_interest, ": down sgRNAs vs NTC"),
      subtitle = paste0(
        "n guides = ", length(guides_down),
        " | log2FC = ", round(lf_down, 3),
        " | % change = ", round(pct_down, 1), "%",
        " | p = ", signif(p_down, 3),
        " | n cells = ", nrow(df_sub_down)
      ),
      x = "",
      y = paste0(tg, " expression")
    )
} else {
  p_cre_down <- ggplot() + theme_void() + labs(title = "No down sgRNAs for this CRE")
}

# ----------------------------
# 6) Panel C: up sgRNAs vs NTC
# ----------------------------
if (nrow(df_sub_up) > 0) {
  plot_df_up <- bind_rows(
    data.frame(group = "Up sgRNAs", expr = df_sub_up$target_gene_expr, stringsAsFactors = FALSE),
    data.frame(group = "NTC", expr = y_ntc, stringsAsFactors = FALSE)
  )

  p_cre_up <- ggplot(plot_df_up, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = paste0(cre_of_interest, ": up sgRNAs vs NTC"),
      subtitle = paste0(
        "n guides = ", length(guides_up),
        " | log2FC = ", round(lf_up, 3),
        " | % change = ", round(pct_up, 1), "%",
        " | p = ", signif(p_up, 3),
        " | n cells = ", nrow(df_sub_up)
      ),
      x = "",
      y = paste0(tg, " expression")
    )
} else {
  p_cre_up <- ggplot() + theme_void() + labs(title = "No up sgRNAs for this CRE")
}

# ----------------------------
# 7) Panel D: per-sgRNA distributions
# ----------------------------
guide_order <- sgRNA_results %>%
  filter(target_name == cre_of_interest) %>%
  arrange(logFC_mean) %>%
  pull(gRNA)

plot_df_guides <- df_cre %>%
  filter(target_name == cre_of_interest) %>%
  mutate(
    gRNA = factor(gRNA, levels = guide_order),
    pool_label = case_when(
      gRNA %in% guides_down ~ "Down pool",
      gRNA %in% guides_up ~ "Up pool",
      TRUE ~ "Not used"
    )
  )

p_cre_guides <- ggplot(plot_df_guides, aes(x = gRNA, y = target_gene_expr, fill = pool_label)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.25, fill = "white") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Down pool" = "blue", "Up pool" = "red", "Not used" = "gray70")) +
  labs(
    title = paste0(cre_of_interest, ": per-sgRNA distributions"),
    subtitle = "Blue = down pool, red = up pool, gray = not used",
    x = "sgRNA",
    y = paste0(tg, " expression"),
    fill = ""
  )

# ----------------------------
# 8) Panel E: sgRNA-level effect sizes + statistics
# ----------------------------
df_sg_effect <- sgRNA_results %>%
  filter(target_name == cre_of_interest) %>%
  mutate(
    gRNA = factor(gRNA, levels = guide_order),
    pool_label = case_when(
      gRNA %in% guides_down ~ "Down pool",
      gRNA %in% guides_up ~ "Up pool",
      TRUE ~ "Not used"
    ),
    stat_label = paste0(
      "p=", signif(p, 2),
      "\nFDR=", signif(p_adj, 2)
    )
  )

p_cre_effect <- ggplot(df_sg_effect, aes(x = gRNA, y = percent_change, color = pool_label)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = stat_label),
    hjust = ifelse(df_sg_effect$percent_change >= 0, -0.1, 1.1),
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("Down pool" = "blue", "Up pool" = "red", "Not used" = "gray70")) +
  labs(
    title = paste0(cre_of_interest, ": sgRNA-level effect sizes"),
    subtitle = "Percent change vs NTC for each sgRNA, with p-value and FDR",
    x = "sgRNA",
    y = "Percent change vs NTC",
    color = ""
  )

print(p_cre_effect)

# ----------------------------
# 9) Combine and save
# ----------------------------
p_cre_compare <- (p_cre_all / p_cre_down / p_cre_up / p_cre_guides / p_cre_effect) +
  plot_layout(heights = c(1, 1, 1, 1.4, 1))

print(p_cre_compare)

ggsave(
  file.path(outdir, paste0(cre_of_interest, "_all_vs_down_vs_up_vs_guides.pdf")),
  p_cre_compare,
  width = 10,
  height = 17
)

# ----------------------------
# 10) Print summary tables
# ----------------------------
print(
  sgRNA_results %>%
    filter(target_name == cre_of_interest) %>%
    arrange(logFC_mean) %>%
    select(target_gene, target_name, gRNA, n_g, mean_expr, mean_ntc, logFC_mean, percent_change, p, p_adj, direction)
)

print(
  cre_results_directional %>%
    filter(target_name == cre_of_interest) %>%
    select(target_gene, target_name, pooled_group, n_cells, mean_expr, mean_ntc, logFC_mean, percent_change, p, p_adj)
)

# ============================================================
# Merge the "up" sgRNAs for one CRE and remake Panel 8
# ============================================================

cre_of_interest <- "SYNGAP1_enhancer_3UTR_1068"

# target gene
df_sub_all <- df_cre %>%
  filter(target_name == cre_of_interest)

if (nrow(df_sub_all) == 0) stop("No cells found for CRE: ", cre_of_interest)

tg <- unique(df_sub_all$target_gene)
tg <- tg[!is.na(tg) & tg != ""]
if (length(tg) == 0) stop("No target_gene found for CRE: ", cre_of_interest)
tg <- tg[1]

# NTC expression
y_ntc <- as.numeric(as.matrix(expr_all_targets[tg, cells_ntc, drop = FALSE]))

# up sgRNAs from your directional pooling definition
guides_up <- cre_pool_def %>%
  filter(target_name == cre_of_interest, pooled_group == "up") %>%
  pull(gRNA) %>%
  unique()

if (length(guides_up) == 0) {
  stop("No up sgRNAs found for CRE: ", cre_of_interest)
}

message("Up sgRNAs used: ", paste(guides_up, collapse = ", "))

# ------------------------------------------------------------
# 1) Individual up sgRNAs
# ------------------------------------------------------------
df_up_individual <- bind_rows(lapply(guides_up, function(gr) {
  x <- df_cre %>%
    filter(target_name == cre_of_interest, gRNA == gr) %>%
    pull(target_gene_expr)

  if (length(x) == 0) return(NULL)

  lf <- log2((mean(x) + 1e-3) / (mean(y_ntc) + 1e-3))
  pct <- 100 * (2^lf - 1)
  pval <- if (length(x) >= 10 && length(y_ntc) >= 10) wilcox.test(x, y_ntc)$p.value else NA_real_

  data.frame(
    label = gr,
    group_type = "Individual up sgRNA",
    n_cells = length(x),
    logFC_mean = lf,
    percent_change = pct,
    p = pval,
    stringsAsFactors = FALSE
  )
}))

# ------------------------------------------------------------
# 2) Merged up pool
# ------------------------------------------------------------
x_up_merged <- df_cre %>%
  filter(target_name == cre_of_interest, gRNA %in% guides_up) %>%
  pull(target_gene_expr)

lf_up_merged <- log2((mean(x_up_merged) + 1e-3) / (mean(y_ntc) + 1e-3))
pct_up_merged <- 100 * (2^lf_up_merged - 1)
p_up_merged <- if (length(x_up_merged) >= 10 && length(y_ntc) >= 10) {
  wilcox.test(x_up_merged, y_ntc)$p.value
} else {
  NA_real_
}

df_up_merged <- data.frame(
  label = paste0("Merged up pool (n guides=", length(guides_up), ")"),
  group_type = "Merged up pool",
  n_cells = length(x_up_merged),
  logFC_mean = lf_up_merged,
  percent_change = pct_up_merged,
  p = p_up_merged,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3) Combine and compute panel-specific FDR
#    NOTE: this FDR is only across the entries shown in this panel.
# ------------------------------------------------------------
df_panel8_up <- bind_rows(df_up_individual, df_up_merged) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    stat_label = paste0(
      "n=", n_cells,
      "\np=", signif(p, 2),
      "\nFDR=", signif(p_adj, 2)
    )
  )

print(df_panel8_up)

# ------------------------------------------------------------
# 4) Order plot: individual guides first, merged pool last
# ------------------------------------------------------------
plot_order <- c(guides_up, paste0("Merged up pool (n guides=", length(guides_up), ")"))

df_panel8_up <- df_panel8_up %>%
  mutate(
    label = factor(label, levels = plot_order),
    group_type = factor(group_type, levels = c("Individual up sgRNA", "Merged up pool"))
  )

# ------------------------------------------------------------
# 5) Recreate Panel 8 with p-value and FDR
# ------------------------------------------------------------
p_panel8_up <- ggplot(df_panel8_up, aes(x = label, y = percent_change, color = group_type)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = stat_label),
    hjust = ifelse(df_panel8_up$percent_change >= 0, -0.1, 1.1),
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("Individual up sgRNA" = "red", "Merged up pool" = "black")) +
  labs(
    title = paste0(cre_of_interest, ": up sgRNAs and merged up pool"),
    subtitle = paste0("Target gene: ", tg, " | compared vs NTC"),
    x = "",
    y = "Percent change vs NTC",
    color = ""
  )

print(p_panel8_up)

ggsave(
  file.path(outdir, paste0(cre_of_interest, "_up_guides_plus_merged_panel8.pdf")),
  p_panel8_up,
  width = 8,
  height = 5
)


# ------------------------------------------------------------
# Recreate Panel 8 with only the merged group
# and add NTC distribution below
# ------------------------------------------------------------

# rename merged label to CRE name
df_panel8_merged <- df_panel8_up %>%
  filter(group_type == "Merged up pool") %>%
  mutate(
    label = cre_of_interest,
    stat_label = paste0(
      "n=", n_cells,
      "\np=", signif(p, 2),
      "\nFDR=", signif(p_adj, 2)
    )
  )

# top panel: merged group only
p_panel8_merged <- ggplot(df_panel8_merged,
                          aes(x = label, y = percent_change)) +
  geom_point(size = 3, color = "black") +
  ylim(0, 50) +
  geom_text(
    aes(label = stat_label),
    hjust = ifelse(df_panel8_merged$percent_change >= 0, -0.1, 1.1),
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() +
  labs(
    title = paste0(cre_of_interest, ": merged up group"),
    subtitle = paste0("Target gene: ", tg, " | compared vs NTC"),
    x = "",
    y = "Percent change vs NTC"
  )

# bottom panel: NTC distribution
df_ntc_plot <- data.frame(
  group = "NTC",
  expr = y_ntc,
  stringsAsFactors = FALSE
)

p_ntc_dist <- ggplot(df_ntc_plot, aes(x = group, y = expr)) +
  geom_violin(fill = "gray70", alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  theme_bw() +
  labs(
    title = "NTC distribution",
    x = "",
    y = paste0(tg, " expression")
  )

# combine
p_panel8_final <- p_panel8_merged / p_ntc_dist +
  plot_layout(heights = c(1, 1.2))

print(p_panel8_final)

ggsave(
  file.path(outdir, paste0(cre_of_interest, "_merged_only_plus_NTC_distribution.pdf")),
  p_panel8_final,
  width = 7,
  height = 8
)
```