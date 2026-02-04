# Initial overlap between brain cCREs and neuro genes TSS +/- 100kb

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
library(ggplot2)
library(tidyplots)
library(dplyr)
library(paletteer)
library(ggrepel)
library(patchwork)

setwd("/home/users/jo912684")
```

## ----------------------------------------------
## STEP 1: Initial overlap for the selected genes of interest
## ----------------------------------------------

```{r}
# Expand genomic coordinates from neuro genes TSS by +/- 100 kb:
df_TargetTouch$start_100kb <- df_TargetTouch$transcript_TSS - 100000
df_TargetTouch$end_100kb <- df_TargetTouch$transcript_TSS + 100000
df_TargetTouch_expanded_100kb <- data.frame(chr= df_TargetTouch$chrom,
                                            start= df_TargetTouch$start_100kb,
                                            end= df_TargetTouch$end_100kb,
                                            gene= df_TargetTouch$gene_name)
df_TargetTouch_expanded_100kb_gr <- makeGRangesFromDataFrame(df_TargetTouch_expanded_100kb, keep.extra.columns = T,
                                                             seqnames.field = "chr", 
                                                             start.field="start", 
                                                             end.field = "end")

# Perform initial overlap with brain cCREs:
final_brain_CREs_gr <- final_brain_CREs_gr[final_brain_CREs_gr@seqnames %in% chromosomes,]
initial_overlap_brain_cCREs <- findOverlaps(df_TargetTouch_expanded_100kb_gr, final_brain_CREs_gr)
initial_overlap_brain_cCREs_info <- data.frame(
  TargetTouch = df_TargetTouch_expanded_100kb_gr[queryHits(initial_overlap_brain_cCREs)],
  CREs = final_brain_CREs_gr[subjectHits(initial_overlap_brain_cCREs)]
)
initial_overlap_brain_cCREs_info <- initial_overlap_brain_cCREs_info[initial_overlap_brain_cCREs_info$CREs.classification != "PLS",]
initial_overlap_brain_cCREs_info <- initial_overlap_brain_cCREs_info[!duplicated(initial_overlap_brain_cCREs_info),] # 30,823 cCREs
initial_overlap_brain_cCREs_info <- droplevels(initial_overlap_brain_cCREs_info)
table(initial_overlap_brain_cCREs_info$CREs.classification)
# ATAC_H3K4me1_H3K27ac = 4903
# CA = 3
# CA-CTCF = 239
# CA-H3K4me3 = 225
# CA-TG = 2
# dELS = 9330
# pELS = 3781
# STARRseq = 18
# TF = 6
# ENCODE_REST = 712
# Jayavelu2020 = 5770
# SilencerDB = 1617
# Ss-STARR-seq= 4217

initial_overlap_brain_cCREs_gr <- makeGRangesFromDataFrame(initial_overlap_brain_cCREs_info, 
                                       keep.extra.columns = T,
                                       seqnames.field = "CREs.seqnames",
                                       start.field = "CREs.start", 
                                       end.field = "CREs.end")

# Add context to each cCRE:
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
exons <- exons(txdb)
introns <- unlist(intronsByTranscript(txdb))
promoters <- promoters(genes, upstream= 2000, downstream = 0)
fiveUTR <- fiveUTRsByTranscript(txdb, use.name=TRUE)
threeUTR <- threeUTRsByTranscript(txdb, use.name=TRUE)
intergenic <- gaps(genes)

initial_overlap_brain_cCREs_gr$context <- "intergenic"
initial_overlap_brain_cCREs_gr$context[overlapsAny(initial_overlap_brain_cCREs_gr, promoters)] <- "promoter"
initial_overlap_brain_cCREs_gr$context[overlapsAny(initial_overlap_brain_cCREs_gr, fiveUTR) & initial_overlap_brain_cCREs_gr$context == "intergenic"] <- "5UTR"
initial_overlap_brain_cCREs_gr$context[overlapsAny(initial_overlap_brain_cCREs_gr, threeUTR) & initial_overlap_brain_cCREs_gr$context == "intergenic"] <- "3UTR"
initial_overlap_brain_cCREs_gr$context[overlapsAny(initial_overlap_brain_cCREs_gr, exons) & initial_overlap_brain_cCREs_gr$context == "intergenic"] <- "exon"
initial_overlap_brain_cCREs_gr$context[overlapsAny(initial_overlap_brain_cCREs_gr, introns) & initial_overlap_brain_cCREs_gr$context == "intergenic"] <- "intron"

initial_overlap_brain_cCREs_df <- as.data.frame(initial_overlap_brain_cCREs_gr)

write.csv(initial_overlap_brain_cCREs_df, file="finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", row.names = FALSE)

saveRDS(initial_overlap_brain_cCREs_gr, file="finding_cCREs/results/initial_overlap_brain_cCREs_100kb.rds")

p_brain_cCREs_initial_search <- 
  initial_overlap_brain_cCREs_df %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(TargetTouch.gene)),
                         fill= CREs.cCRE_type)) +
  coord_flip() +
  ylab("N brain cCREs overlapping +/-100kb from TSS") +
  xlab("") +
  scale_fill_manual(values = c("darkred","darkblue"))
ggsave(p_brain_cCREs_initial_search, filename= "finding_cCREs/results/plots/Brain_cCREs_initial_search_neuro_genes_100kb_cCREtype.pdf",
       width = 8, height= 15)

p_brain_cCREs_initial_search2 <- 
  initial_overlap_brain_cCREs_df %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(TargetTouch.gene)),
                         fill= context)) +
  coord_flip() +
  ylab("N brain cCREs overlapping +/-100kb from TSS") +
  xlab("")
ggsave(p_brain_cCREs_initial_search2, filename= "finding_cCREs/results/plots/Brain_cCREs_initial_search_neuro_genes_100kb_cCREtype.pdf",
       width = 8, height= 15)

mean(table(initial_overlap_brain_cCREs_df$TargetTouch.gene)) # 248.5726 cCREs per gene
max(table(initial_overlap_brain_cCREs_df$TargetTouch.gene)) # 624 cCREs per gene
min(table(initial_overlap_brain_cCREs_df$TargetTouch.gene)) # 52 cCREs per gene
```

## ----------------------------------------------
## STEP 2: Summary plots and stats on initial overlap
## ----------------------------------------------

```{r}
df <- read.csv("finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", head=T, stringsAsFactors = T)

all_sum <- as.data.frame(table(df$TargetTouch.gene))
p1 <- ggplot(all_sum, aes(x = "", y = Freq)) + 
  geom_violin(fill = "white", width=0.2) + 
  geom_boxplot(fill= "gray", width=0.1, outlier.shape = NA) +
  geom_jitter(width = 0.01) + 
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  annotate("label", x = 1.15, y = mean(all_sum$Freq), label = paste("Mean:", round(mean(all_sum$Freq), 0)),
           hjust = 0, vjust = 0.5, size = 4, fontface = "bold",
           color = "red") +
  xlab("") + ylab("Number of cCREs") + ggtitle(label="123 neuro genes of interest")

sum_type <- as.data.frame(table(df$CREs.cCRE_type))
sum_type <- sum_type %>%
  mutate(Proportion = Freq / sum(Freq))

p2 <- ggplot(sum_type, aes(x = Var1, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "cCRE type", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values=c("#A80A20", "#336392")) +
  ggtitle(label="Classification of 30,823 cCREs")

sum_context <- as.data.frame(table(df$context))
sum_context <- sum_context %>%
  arrange(desc(-Proportion)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))

p3 <- ggplot(sum_context, aes(x = Var1, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "cCRE context", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values=c(paletteer_d("ggthemes::calc", length(sum_context$Var1)))) +
  ggtitle(label="Context of 30,823 cCREs") + coord_flip()

ggsave("finding_cCREs/results/plots/cCREs_summary_plots.pdf", p1 + p2 + p3, width = 10, height= 6)

df$eQTL_PheWAS <- ifelse(df$CRE_contains_eQTL == "yes" | df$CRE_contains_PheWAS == "yes", "yes", "no")
sum_eQTL <- as.data.frame(table(df$eQTL_PheWAS))
sum_eQTL <- sum_eQTL %>%
  mutate(Proportion = Freq / sum(Freq))

p4 <- ggplot(sum_eQTL, aes(x = Var1, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "cCRE contains eQTL/PheWAS hit", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values=c("darkgray", "gold"))

df_eQTL <- df[df$eQTL_PheWAS == "yes",]
df_eQTL$eQTL_PheWAS_same <- ifelse(df_eQTL$eQTL_gene_same_as_target_gene == "yes" | df_eQTL$PheWAS_gene_same_as_target_gene == "yes", "yes", "no")
sum_eQTL2 <- as.data.frame(table(df_eQTL$eQTL_PheWAS_same))
sum_eQTL2 <- sum_eQTL2 %>%
  mutate(Proportion = Freq / sum(Freq))

p5 <- ggplot(sum_eQTL2, aes(x = Var1, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "eQTL/PheWAS hit to target gene", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values=c("#926C15", "#FFC300"))

sum_cc <- as.data.frame(table(df$chrom_contact_target_gene))
sum_cc <- sum_cc %>%
  mutate(Proportion = Freq / sum(Freq))

p6 <- ggplot(sum_cc, aes(x = Var1, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "cCRE contains contact observation", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values=c("darkgray", "gold"))

ggsave("finding_cCREs/results/plots/cCREs_summary_plots_additionalhits.pdf", p4 + p5 + p6, width = 10, height= 6)
```