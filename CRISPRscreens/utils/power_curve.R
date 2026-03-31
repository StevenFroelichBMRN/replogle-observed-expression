#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(ggplot2)
})

# ----------------------------
# Minimal CLI parser (no optparse)
# Supports: --key value
# ----------------------------
parse_args <- function(args) {
  if (length(args) == 0) return(list())
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }
  keys <- args[seq(1, length(args), by = 2)]
  vals <- args[seq(2, length(args), by = 2)]
  if (any(!grepl("^--", keys))) {
    stop("All keys must start with '--' (e.g., --seurat path.rds).", call. = FALSE)
  }
  keys <- sub("^--", "", keys)
  out <- as.list(vals)
  names(out) <- keys
  out
}

as_int <- function(x, default=NULL) if (is.null(x)) default else as.integer(x)
as_num <- function(x, default=NULL) if (is.null(x)) default else as.numeric(x)
as_chr <- function(x, default=NULL) if (is.null(x)) default else as.character(x)

split_genes <- function(x) {
  genes <- unlist(strsplit(x, ","))
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  unique(genes)
}

safe_mkdir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Layer helpers (Seurat v5 / BPCells-friendly)
# ----------------------------
get_layers_safe <- function(obj, assay) {
  lyr <- NULL

  # SeuratObject v5 Layers()
  lyr <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (!is.null(lyr) && length(lyr) > 0) return(as.character(lyr))

  # @layers slot
  lyr <- tryCatch(names(obj[[assay]]@layers), error = function(e) NULL)
  if (!is.null(lyr) && length(lyr) > 0) return(as.character(lyr))

  # fallback
  lyr <- tryCatch(names(obj[[assay]]), error = function(e) NULL)
  if (!is.null(lyr) && length(lyr) > 0) return(as.character(lyr))

  character(0)
}

choose_layer <- function(available, requested = "auto") {
  if (length(available) == 0) stop("No layers found for the assay.", call. = FALSE)

  if (requested != "auto") {
    if (!requested %in% available) {
      stop(sprintf("Requested layer '%s' not found. Available layers: %s",
                   requested, paste(available, collapse = ", ")),
           call. = FALSE)
    }
    return(requested)
  }

  preferred <- c("counts", "Counts", "raw", "RAW", "umis", "UMIs")
  hit <- intersect(preferred, available)
  if (length(hit) > 0) return(hit[1])

  hit2 <- available[grepl("count|raw|umi", available, ignore.case = TRUE)]
  if (length(hit2) > 0) return(hit2[1])

  if ("data" %in% available) return("data")

  available[1]
}

# ----------------------------
# Power helpers
# ----------------------------
analytic_n_for_effect <- function(pct_change, sd_log2, alpha, target_power) {
  fc <- 1 + pct_change/100
  delta_log2 <- log2(fc)

  p_out <- tryCatch(
    power.t.test(
      delta = delta_log2,
      sd = sd_log2,
      sig.level = alpha,
      power = target_power,
      type = "two.sample",
      alternative = "two.sided"
    ),
    error = function(e) NULL
  )

  if (is.null(p_out)) {
    tibble(
      percent_change = pct_change,
      fold_change = fc,
      delta_log2 = delta_log2,
      n_cells_per_group = NA_real_,
      n_cells_total = NA_real_,
      note = "power.t.test failed"
    )
  } else {
    n_pg <- ceiling(p_out$n)
    tibble(
      percent_change = pct_change,
      fold_change = fc,
      delta_log2 = delta_log2,
      n_cells_per_group = n_pg,
      n_cells_total = 2 * n_pg,
      note = NA_character_
    )
  }
}

estimate_power_sim_multiplicative_cpm <- function(x_ctrl_cpm, fold_change, n_pg, alpha, reps,
                                                 test = c("t", "wilcox")) {
  test <- match.arg(test)

  pvals <- replicate(reps, {
    ctrl_cpm <- sample(x_ctrl_cpm, size = n_pg, replace = TRUE)
    pert_cpm <- sample(x_ctrl_cpm, size = n_pg, replace = TRUE) * fold_change

    ctrl <- log2(ctrl_cpm + 1)
    pert <- log2(pert_cpm + 1)

    tryCatch({
      if (test == "t") {
        t.test(pert, ctrl, alternative = "two.sided")$p.value
      } else {
        wilcox.test(pert, ctrl, alternative = "two.sided")$p.value
      }
    }, error = function(e) NA_real_)
  })

  mean(pvals < alpha, na.rm = TRUE)
}

find_n_sim_multiplicative_cpm <- function(x_ctrl_cpm, fold_change, alpha, target_power, reps,
                                         n_min, n_max, test = c("t", "wilcox")) {
  test <- match.arg(test)

  p_min <- estimate_power_sim_multiplicative_cpm(x_ctrl_cpm, fold_change, n_min, alpha, reps, test)
  if (is.finite(p_min) && p_min >= target_power) {
    return(list(n = n_min, power_at_n = p_min))
  }

  p_max <- estimate_power_sim_multiplicative_cpm(x_ctrl_cpm, fold_change, n_max, alpha, reps, test)
  if (!is.finite(p_max) || p_max < target_power) {
    return(list(n = NA_integer_, power_at_n = p_max))
  }

  lo <- n_min
  hi <- n_max
  best_n <- hi
  best_p <- p_max

  while ((hi - lo) > 1) {
    mid <- floor((lo + hi) / 2)
    p_mid <- estimate_power_sim_multiplicative_cpm(x_ctrl_cpm, fold_change, mid, alpha, reps, test)
    if (is.finite(p_mid) && p_mid >= target_power) {
      best_n <- mid
      best_p <- p_mid
      hi <- mid
    } else {
      lo <- mid
    }
  }

  list(n = best_n, power_at_n = best_p)
}

# ----------------------------
# Read CLI args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(args)

# Required
seurat_rds <- as_chr(opt$seurat, NULL)
genes_str  <- as_chr(opt$genes, NULL)
target_col <- as_chr(opt$target_col, NULL)
ctrl_label <- as_chr(opt$ctrl_label, NULL)
outdir     <- as_chr(opt$outdir, NULL)

if (is.null(seurat_rds) || is.null(genes_str) || is.null(target_col) ||
    is.null(ctrl_label) || is.null(outdir)) {
  stop(
    paste0(
      "Missing required args.\nExample:\n",
      "Rscript power_curve.R --seurat scobj.rds --genes APP,SCN2A ",
      "--target_col target --ctrl_label NTC --outdir out\n\n",
      "Optional:\n",
      "  --assay RNA --layer auto --alpha 0.05 --power 0.8 --min_pct 5 --max_pct 100 --step_pct 5\n",
      "  --sim_reps 500 --sim_test wilcox --n_min 10 --n_max 50000 --sim_seed 1\n"
    ),
    call. = FALSE
  )
}

# Optional
assay        <- as_chr(opt$assay, "RNA")
layer_req    <- as_chr(opt$layer, "auto")  # NEW: --layer auto|counts|data|...
alpha        <- as_num(opt$alpha, 0.05)
target_power <- as_num(opt$power, 0.80)
min_pct      <- as_int(opt$min_pct, 5)
max_pct      <- as_int(opt$max_pct, 100)
step_pct     <- as_int(opt$step_pct, 5)

sim_reps     <- as_int(opt$sim_reps, 200)
n_min        <- as_int(opt$n_min, 10)
n_max        <- as_int(opt$n_max, 50000)
sim_seed     <- as_int(opt$sim_seed, 1)
sim_test     <- as_chr(opt$sim_test, "t")
if (!sim_test %in% c("t", "wilcox")) stop("--sim_test must be 't' or 'wilcox'", call. = FALSE)

genes_of_interest <- split_genes(genes_str)
effect_percents <- seq(min_pct, max_pct, by = step_pct)
if (length(effect_percents) < 2) stop("Bad percent-change range; adjust min/max/step.", call. = FALSE)

# Folder structure
log2_dir <- file.path(outdir, "log2_CPM")
sim_dir  <- file.path(outdir, "simulation_based")
safe_mkdir(outdir)
safe_mkdir(log2_dir)
safe_mkdir(sim_dir)

# ----------------------------
# Load object and subset controls (robust, avoids FetchData eval)
# ----------------------------
obj <- readRDS(seurat_rds)

if (!target_col %in% colnames(obj@meta.data)) {
  stop(sprintf("Metadata column '%s' not found in obj@meta.data", target_col), call. = FALSE)
}

cells_ctrl <- rownames(obj@meta.data)[obj@meta.data[[target_col]] == ctrl_label]
if (length(cells_ctrl) == 0) {
  stop(sprintf("No cells found with %s == '%s'", target_col, ctrl_label), call. = FALSE)
}

obj_ctrl <- subset(obj, cells = cells_ctrl)

message(sprintf("Total cells in object: %s", format(ncol(obj), big.mark=",")))
message(sprintf("Control cells (%s == %s): %s", target_col, ctrl_label, format(ncol(obj_ctrl), big.mark=",")))

DefaultAssay(obj_ctrl) <- assay

# ----------------------------
# Resolve layer to use (auto-detect if needed)
# ----------------------------
available_layers <- get_layers_safe(obj_ctrl, assay)
message("Available layers in assay '", assay, "': ", paste(available_layers, collapse = ", "))

layer_use <- choose_layer(available_layers, requested = layer_req)
message("Using layer: ", layer_use)

# ----------------------------
# Get data layer (sparse) and compute CPM for genes
# ----------------------------
mat_all <- GetAssayData(object = obj_ctrl, assay = assay, layer = layer_use)  # genes x cells

# CPM denominator should be total molecules per cell. This only truly corresponds to CPM if 'layer_use' is raw counts.
# If layer_use is normalized data, CPM is not meaningful; script will still run but interpret cautiously.
col_totals_num <- as.numeric(Matrix::colSums(mat_all))
col_totals_num[col_totals_num == 0] <- NA_real_

genes_present <- intersect(genes_of_interest, rownames(obj_ctrl[[assay]]))
genes_missing <- setdiff(genes_of_interest, genes_present)
if (length(genes_present) == 0) stop("None of the requested genes were found in the selected assay.", call. = FALSE)
if (length(genes_missing) > 0) warning("Skipping missing genes: ", paste(genes_missing, collapse=", "))
message("Genes found (will analyze): ", paste(genes_present, collapse=", "))

mat_sub <- mat_all[genes_present, , drop = FALSE]
mat_sub_dense <- as.matrix(mat_sub)  # genes x cells; ok for modest gene list

# Treat mat_sub_dense as "counts-like" for CPM computation
cpm_mat <- t( t(mat_sub_dense) / col_totals_num ) * 1e6
log2cpm_mat <- log2(cpm_mat + 1)

gene_summaries <- tibble(
  gene = genes_present,
  n_cells = apply(log2cpm_mat, 1, function(x) sum(is.finite(x))),
  mean_log2cpm = apply(log2cpm_mat, 1, function(x) mean(x, na.rm=TRUE)),
  sd_log2cpm   = apply(log2cpm_mat, 1, function(x) sd(x, na.rm=TRUE)),
  mean_cpm     = apply(cpm_mat, 1, function(x) mean(x, na.rm=TRUE)),
  sd_cpm       = apply(cpm_mat, 1, function(x) sd(x, na.rm=TRUE))
)

write_csv(gene_summaries, file.path(outdir, "control_CPM_and_log2CPM_variability_summary.csv"))

# ----------------------------
# Per-gene outputs: analytic + simulation
# ----------------------------
set.seed(sim_seed)

all_analytic <- list()
all_simreq   <- list()

for (g in genes_present) {
  message("\n=== Gene: ", g, " ===")

  g_log2_dir <- file.path(log2_dir, g)
  g_sim_dir  <- file.path(sim_dir, g)
  safe_mkdir(g_log2_dir)
  safe_mkdir(g_sim_dir)

  x_ctrl_cpm  <- as.numeric(cpm_mat[g, ]);     x_ctrl_cpm  <- x_ctrl_cpm[is.finite(x_ctrl_cpm)]
  x_ctrl_log2 <- as.numeric(log2cpm_mat[g, ]); x_ctrl_log2 <- x_ctrl_log2[is.finite(x_ctrl_log2)]

  if (length(x_ctrl_cpm) < 10 || length(x_ctrl_log2) < 10) {
    warning("Skipping ", g, ": too few finite values.")
    next
  }

  sdg <- sd(x_ctrl_log2)
  if (!is.finite(sdg) || sdg == 0) {
    warning("Skipping ", g, ": SD is zero/NA on log2(CPM+1) scale.")
    next
  }

  # 1) Analytic curve
  analytic_curve <- map_dfr(effect_percents, ~analytic_n_for_effect(.x, sdg, alpha, target_power)) %>%
    mutate(
      gene = g,
      sd_log2 = sdg,
      alpha = alpha,
      target_power = target_power,
      control_cells_available = length(x_ctrl_log2),
      layer_used = layer_use
    ) %>%
    select(gene, percent_change, fold_change, delta_log2, sd_log2, alpha, target_power,
           n_cells_per_group, n_cells_total, control_cells_available, layer_used, note)

  write_csv(analytic_curve, file.path(g_log2_dir, paste0("analytic_power_curve_", g, ".csv")))

  p_analytic <- ggplot(analytic_curve, aes(x = percent_change, y = n_cells_per_group)) +
    geom_line() + geom_point() +
    scale_y_log10() +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(g, " — Analytic required cells vs percent change"),
      subtitle = sprintf("power.t.test on log2(CPM+1); power=%.0f%% alpha=%.02f SD=%.3f; layer=%s",
                         target_power * 100, alpha, sdg, layer_use),
      x = "Percent change in mean expression",
      y = "Required cells per group (log10 scale)"
    ) +
    geom_hline(yintercept = length(x_ctrl_log2), linetype = "dashed") +
    annotate("text", x = max(effect_percents), y = length(x_ctrl_log2),
             vjust = -0.5, hjust = 1, size = 3.5,
             label = paste0("Control cells: ", length(x_ctrl_log2)))

  ggsave(file.path(g_log2_dir, paste0("analytic_power_curve_", g, ".png")), p_analytic, width = 8, height = 5, dpi = 300)
  ggsave(file.path(g_log2_dir, paste0("analytic_power_curve_", g, ".pdf")), p_analytic, width = 8, height = 5)

  all_analytic[[g]] <- analytic_curve

  # 2) Simulation-based (multiplicative on CPM -> test on log2(CPM+1))
  sim_required <- map_dfr(effect_percents, function(pct) {
    fc <- 1 + pct / 100
    found <- find_n_sim_multiplicative_cpm(
      x_ctrl_cpm = x_ctrl_cpm,
      fold_change = fc,
      alpha = alpha,
      target_power = target_power,
      reps = sim_reps,
      n_min = n_min,
      n_max = n_max,
      test = sim_test
    )

    tibble(
      gene = g,
      percent_change = pct,
      fold_change = fc,
      delta_log2 = log2(fc),
      alpha = alpha,
      target_power = target_power,
      sim_reps = sim_reps,
      sim_test = sim_test,
      n_min = n_min,
      n_max = n_max,
      required_n_cells_per_group = found$n,
      achieved_power_at_required_n = found$power_at_n,
      control_cells_available = length(x_ctrl_cpm),
      layer_used = layer_use,
      sim_model = "multiplicative_on_CPM_then_test_on_log2CPMplus1"
    )
  })

  write_csv(sim_required, file.path(g_sim_dir, paste0("simulation_required_n_", g, ".csv")))

  p_sim <- ggplot(sim_required, aes(x = percent_change, y = required_n_cells_per_group)) +
    geom_line() + geom_point() +
    scale_y_log10() +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(g, " — Simulation-based required cells vs percent change"),
      subtitle = sprintf("CPM×FC then test log2(CPM+1); reps=%d test=%s; power=%.0f%% alpha=%.02f; layer=%s",
                         sim_reps, sim_test, target_power * 100, alpha, layer_use),
      x = "Percent change in mean expression",
      y = "Required cells per group (log10 scale)"
    ) +
    geom_hline(yintercept = length(x_ctrl_cpm), linetype = "dashed") +
    annotate("text", x = max(effect_percents), y = length(x_ctrl_cpm),
             vjust = -0.5, hjust = 1, size = 3.5,
             label = paste0("Control cells: ", length(x_ctrl_cpm)))

  ggsave(file.path(g_sim_dir, paste0("simulation_required_n_", g, ".png")), p_sim, width = 8, height = 5, dpi = 300)
  ggsave(file.path(g_sim_dir, paste0("simulation_required_n_", g, ".pdf")), p_sim, width = 8, height = 5)

  all_simreq[[g]] <- sim_required
}

# Combined outputs
if (length(all_analytic) > 0) {
  combined_analytic <- bind_rows(all_analytic) %>% arrange(gene, percent_change)
  write_csv(combined_analytic, file.path(outdir, "analytic_power_curves_all_genes.csv"))
}
if (length(all_simreq) > 0) {
  combined_sim <- bind_rows(all_simreq) %>% arrange(gene, percent_change)
  write_csv(combined_sim, file.path(outdir, "simulation_required_n_all_genes.csv"))
}

message("\nDone.")
message("Outputs written to: ", normalizePath(outdir))
message("Analytic per-gene folders: ", normalizePath(log2_dir))
message("Simulation per-gene folders: ", normalizePath(sim_dir))
message("Layer used: ", layer_use)

