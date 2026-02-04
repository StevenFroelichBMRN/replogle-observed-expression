# Use chromatin-contact information to highlight elements in regions with known interactions with TSS of genes of interest.

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
## STEP 1: Use chromatin-contact information from NGN2-derived neurons and overlap with the identified cCREs
## ----------------------------------------------

```{r}
# Load chromatin contact information from NGN2-derived neurons:
df_cc_ngn2 <- read.delim("finding_cCREs/data/Nature2025_Human_accelerated_neurons/hg38.Human.merged.5k.2.sig3Dinteractions.plot.bedpe.sorted.txt", head=T, sep="") # 445,108 contacts

df_cc_neu <- read.delim("finding_cCREs/data/Science2019_brain_celltype_specific_interactome/NeuN.5k_interactions_ucsc_genome_browser.inter.bed", head=F, sep="\t")
colnames(df_cc_neu) <- c("chrom", "chromStart", "chromEnd", "name", "score", "value", "exp", "color", "region1Chr", "region1Start", "region1End", "regionName", "region1Strand", "region2Chrom", "region2Start", "region2End", "region2Name", "region2Strand") # 93,290 contacts

# Combine contact datasets to gather all neuronal-related observed contacts:

df_cc_combined <- data.frame(chr_anchor = c(df_cc_ngn2$chr_anchor, df_cc_neu$region1Chr),
                             start_anchor = c(df_cc_ngn2$start_bin_anchor, df_cc_neu$region1Start),
                             end_anchor = c(df_cc_ngn2$end_bin_anchor, df_cc_neu$region1End),
                             chr_contact = c(df_cc_ngn2$chr_contact, df_cc_neu$region2Chrom),
                             start_contact = c(df_cc_ngn2$start_bin_contact, df_cc_neu$region2Start),
                             end_contact = c(df_cc_ngn2$end_bin_contact, df_cc_neu$region2End)) # 538,398 observed genomic contacts

df_cc_anchors_gr <- makeGRangesFromDataFrame(df_cc_combined, seqnames.field = "chr_anchor", start.field = "start_anchor", end.field = "end_anchor", keep.extra.columns = T)

# Get coordinates of neuro genes of interest:
neuro_genes_TSS_gr <- makeGRangesFromDataFrame(df_TargetTouch[,-10], keep.extra.columns = T, seqnames.field = "chrom", start.field="transcript_TSS", end.field = "transcript_TSS")

# Subset anchors to 5bk bins that are included in the neuro genes TSS:
anchors_overlapping_TSS <- findOverlaps(df_cc_anchors_gr, neuro_genes_TSS_gr)
anchors_overlapping_TSS_info <- data.frame(
  Anchors = df_cc_anchors_gr[queryHits(anchors_overlapping_TSS)],
  TargetTouch = neuro_genes_TSS_gr[subjectHits(anchors_overlapping_TSS)]) # 2,156 observed contacts with bins that include the TSS of the neuro genes

# Grab the bins contacted by the bins that include the TSS of the genes of interest:
bins_contacted_TSS_neuro_gr <- makeGRangesFromDataFrame(anchors_overlapping_TSS_info, seqnames.field = "Anchors.chr_contact", start.field="Anchors.start_contact", end.field="Anchors.end_contact", keep.extra.columns = T)

# Which cCREs are in regions with observed contacts to the 5kb bins that include the neuro genes TSS?

## Obtain the overlap between contacted bins and brain cCREs:
overlap_cCREs_contacts <- findOverlaps(bins_contacted_TSS_neuro_gr, initial_overlap_brain_cCREs_gr)
overlap_cCREs_contacts_info <- data.frame(
  Contacts= bins_contacted_TSS_neuro_gr[queryHits(overlap_cCREs_contacts)],
  CREs = initial_overlap_brain_cCREs_gr[subjectHits(overlap_cCREs_contacts)])
overlap_cCREs_contacts_info <- overlap_cCREs_contacts_info[!duplicated(overlap_cCREs_contacts_info),] # 3,725 observed contacts with CREs found in 5kb bins contacting neuro genes TSS

overlap_cCREs_contacts_gr <- makeGRangesFromDataFrame(overlap_cCREs_contacts_info,
                                                      seqnames.field = "CREs.seqnames",
                                                      start.field = "CREs.start",
                                                      end.field= "CREs.end",
                                                      keep.extra.columns = T)
write.csv(overlap_cCREs_contacts_info, file="finding_cCREs/results/Overlap_brain_CREs_100kb_neuron_chromatin_contacts.csv", row.names = FALSE)
saveRDS(overlap_cCREs_contacts_gr, file="finding_cCREs/results/overlap_cCREs_contacts_gr.rds")

p_brain_cCREs_obs_interaction_neuro_genes <- 
  overlap_cCREs_contacts_info %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(Contacts.TargetTouch.gene_name)))) +
  coord_flip() +
  ylab("N cCREs in bins with contacts to TSS") +
  xlab("")
ggsave(p_brain_cCREs_obs_interaction_neuro_genes, filename= "finding_cCREs/results/plots/Brain_cCREs_interactions_NeuroGenes_bins.pdf",
       width = 8, height= 12)

## Add information of contact to the cCREs database:
initial_overlap_brain_cCREs_gr$chrom_contact_target_gene <- "no"
initial_overlap_brain_cCREs_gr$chrom_contact_target_gene[overlapsAny(initial_overlap_brain_cCREs_gr, overlap_cCREs_contacts_gr)] <- "yes"

## Save:
write.csv(as.data.frame(initial_overlap_brain_cCREs_gr), file="finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", row.names = FALSE)

saveRDS(initial_overlap_brain_cCREs_gr, file="finding_cCREs/results/initial_overlap_brain_cCREs_100kb.rds")
```