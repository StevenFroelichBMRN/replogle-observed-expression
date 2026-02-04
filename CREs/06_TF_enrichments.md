# Search for enrichment for TF binding motifs in cCREs with eQTL hits

## ----------------------------------------------
## STEP 0: load libraries and working directory
## ----------------------------------------------

```{r, warning=FALSE}
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)
library(tidyplots)
library(circlize)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(httr)
library(readr)
library(biomaRt) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2024)
library(RSQLite)
library(SummarizedExperiment)
library(TFBSTools) 
library(motifmatchr) 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(phastCons100way.UCSC.hg38)
library(phastCons30way.UCSC.hg38)

setwd("/home/users/jo912684")
```

## ----------------------------------------------
## STEP 1: Use the JASPAR database to find enrichment of motif-changing eQTLs within our identified cCREs
## ----------------------------------------------

```{r}
# Download JASPAR database:
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))

# Load TF classification from DelRosso et al 2023 (https://www.nature.com/articles/s41586-023-05906-y)
motifs_info <- read.csv("finding_cCREs/data/JASPAR2024_CORE_non-redundant-pfms_info.csv", head=T)

# Get genomic sequences for all regions of interest:
genome <- BSgenome.Hsapiens.UCSC.hg38

# Load cCREs with the same targetting gene as the TargetTouch gene
CREs_overlap_eQTl_TF_gr <- makeGRangesFromDataFrame(overlap_cCREs_eQTL_info[overlap_cCREs_eQTL_info$same_gene == "yes",],
                                                    seqnames.field = "CREs.seqnames",
                                                    start.field = "CREs.start",
                                                    end.field= "CREs.end",
                                                    keep.extra.columns = T)

# Get TF motifs
# See: https://jaspar.elixir.no

motifs24 <- TFBSTools::getMatrixSet(sq24, list(
  species = "Homo sapiens",
  collection = "CORE"))
motifs_pwm <- toPWM(motifs24) # convert from 'position frequency matrix (PFM)' to 'position weight matrix (PWM)'

# Does the variant (eQTL) disrupt the TF binding?
CREs_overlap_eQTL_df <- as.data.frame(CREs_overlap_eQTl_TF_gr)

variants_df <- data.frame(
  chromosome = CREs_overlap_eQTL_df$seqnames,
  start = CREs_overlap_eQTL_df$start,
  end = CREs_overlap_eQTL_df$end,
  variant_id = CREs_overlap_eQTL_df$eQTLs.eqtl_id,
  variant_location = CREs_overlap_eQTL_df$eQTLs.start,
  reference_allele = CREs_overlap_eQTL_df$eQTLs.eqtl_ref,
  alternative_allele = CREs_overlap_eQTL_df$eQTLs.eqtl_alt,
  gene= CREs_overlap_eQTL_df$eQTLs.eqtl_gene
)
variants_df_unique <- variants_df[!duplicated(variants_df), ] # remove duplicated rows

# Get the sequences centered at the eQTL variant
flank <- 10
ref_gr <- GRanges(
  seqnames = variants_df_unique$chromosome,
  ranges = IRanges(
    start = pmax(1, variants_df_unique$variant_location - flank),
    end = variants_df_unique$variant_location + flank
  )
)
ref_seqs <- getSeq(genome, ref_gr)

# Modify sequences to the alternative allele
alt_seqs <- DNAStringSet(ref_seqs)
pos_in_seq <- variants_df_unique$variant_location - start(ref_gr) + 1
for (i in seq_along(alt_seqs)) {
  alt_seqs[[i]] <- replaceAt(alt_seqs[[i]], IRanges(pos_in_seq[i], pos_in_seq[i]), variants_df_unique$alternative_allele[i])
}

# Run matchMotifs in batch to identify TFs that bind
ref_scores <- matchMotifs(motifs_pwm, ref_seqs, genome = genome, out = "scores")
saveRDS(ref_scores, file="finding_cCREs/results/ref_scores.rds")

alt_scores <- matchMotifs(motifs_pwm, alt_seqs, genome = genome, out = "scores")
saveRDS(alt_scores, file="finding_cCREs/results/alt_scores.rds")

# Evaluate binding motifs found in REF/ALT but not in the other

# Convert scores to matrices
ref_mat <- assay(ref_scores)
alt_mat <- assay(alt_scores)

# Calculate score differences
score_diff <- alt_mat - ref_mat

# Define a threshold for meaningful change
threshold <- 0.1

# Create a summary data frame
summary_list <- lapply(seq_len(nrow(score_diff)), function(i) {
  variant_id <- variants_df_unique$variant_id[i]  #
  ref_row <- ref_mat[i, ]
  alt_row <- alt_mat[i, ]
  diff_row <- score_diff[i, ]
  
  changed_motifs <- which(abs(diff_row) >= threshold)
  
  if (length(changed_motifs) == 0) return(NULL)
  
  data.frame(
    variant_id = variant_id,
    motif_id = colnames(score_diff)[changed_motifs],
    ref_score = ref_row[changed_motifs],
    alt_score = alt_row[changed_motifs],
    score_change = diff_row[changed_motifs],
    binding_change = ifelse(diff_row[changed_motifs] > 0, "gain", "loss"),
    stringsAsFactors = FALSE
  )
})

# Combine into one data frame
summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df %>%
  rownames_to_column(var = "motif_base")
summary_df$motif_base <- str_replace(summary_df$motif_base, "\\..*", "")
  
# Save summary
write.csv(summary_df, "finding_cCREs/results/motif_binding_changes_summary.csv", row.names = FALSE)

# Join results to TF classification data:
merged_TF_results <- left_join(summary_df, motifs_info, by= c("motif_base" = "Mat"))

write.csv(merged_TF_results, file="finding_cCREs/results/TFbinding_cCREs_eQTL.csv", row.names = FALSE)

# Count number of observations per combination
merged_TF_results <- as_tibble(merged_TF_results)

# Summarize counts
plot_data <- merged_TF_results %>%
  group_by(Classification, binding_change) %>%
  summarise(n = n(), .groups = "drop")

# Create the bar plot
pdf("finding_cCREs/results/plots/TF_binding_results_summary.pdf")
ggplot(plot_data, aes(x = Classification, y = n, fill = binding_change)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "TF Binding Changes by Classification",
    x = "Classification",
    y = "Number of Observations",
    fill = "Binding Change"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# How are the betas for the eQTLs?
merged_TF_results_eQTLsinfo <- left_join(merged_TF_results, CREs_overlap_eQTL_df, by=c("variant_id" = "eQTLs.eqtl_id"))

# Plot:

pdf("finding_cCREs/results/plots/TF_binding_prediction_eQTLs.pdf")
ggplot(merged_TF_results_eQTLsinfo, aes(x= score_change, y= eQTLs.eqtl_beta)) +
  geom_point() +
  theme_bw() +
  xlab("Match score difference [Alternative - Reference]") +
  ylab("eQTL variant beta") +
  facet_wrap(~binding_change + Classification) +
  geom_hline(yintercept = 0)
dev.off()
```