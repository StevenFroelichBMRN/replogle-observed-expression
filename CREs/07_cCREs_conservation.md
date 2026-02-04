# Evaluate conservation of the identified cCREs

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
## STEP 1: Evaluate evolutionary conservation
## ----------------------------------------------

```{r}
# 30 mammals:
# Load conservation scores:
phast30 <- getGScores("phastCons30way.UCSC.hg38")

# Get conservation scores across CREs:
scores30 <- score(phast30, CREs_overlap_eQTl_TF_gr)

# Plot:
pdf("finding_cCREs/results/plots/CREs_overlap_eQTL_PhastCons30.pdf")
hist(scores30, breaks = 100, main= "PhastCons30 scores for CREs", xlab= "Conservation score")
dev.off()

# 7 mammals:
# Load conservation scores:
phast7 <- getGScores("phastCons7way.UCSC.hg38")

# Get conservation scores across CREs:
scores7 <- score(phast7, CREs_overlap_eQTl_TF_gr)

# Plot:
pdf("finding_cCREs/results/plots/CREs_overlap_eQTL_PhastCons7.pdf")
hist(scores7, breaks = 100, main= "PhastCons7 scores for CREs", xlab= "Conservation score")
dev.off()

# Add conservation scores to CREs result:
CREs_overlap_eQTL_df$phast7 <- scores7
df_combined <- left_join(disruption_results, CREs_overlap_eQTL_df, by=c("variant_id" ="eQTLs.eqtl_id"))

pdf("finding_cCREs/results/plots/Test_eQTL_motif_analysis_betas_phast7.pdf")
ggplot(df_combined[df_combined$delta != 0,], aes(x= delta, y= beta, size= phast7, color= Classification)) +
  geom_point() +
  theme_bw() +
  xlab("Match score difference [Alternative - Reference]") +
  ylab("eQTL variant beta")
dev.off()
```