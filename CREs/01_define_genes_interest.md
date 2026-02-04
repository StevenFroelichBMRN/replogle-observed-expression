# Basic details on expression patterns of genes of interest

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
## STEP 1: create list of human genes
## ----------------------------------------------

```{r}
# Prepare a genes x ENSEMBL ID file for human GRC38 reference:
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
human_genes <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(human_genes) <- c("GeneName", "ENSEMBLID")
```

## ----------------------------------------------
## STEP 2: define genes of interest for CRISPR screens.
## ----------------------------------------------

We will start with the neuro-related genes from the TargetTouch database (curated by the Genomics team). Note: the canonical isoform of each gene will be used for their genomic coordinates.

```{r}
# Load information of the neuro genes of interest.
df_TargetTouch <- read.csv("finding_cCREs/data/TargetTouch_NeuroGenes_expression.csv", head=T)

# Plot their distribution across the human genome:
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
pdf("finding_cCREs/results/plots/CircusPlot_NeuroGenes_distribution_hg38.pdf")
circos.initializeWithIdeogram(chromosome.index = chromosomes, species= "hg38")
for (i in 1:nrow(df_TargetTouch)){
circos.points(
    x= (df_TargetTouch$chromStart[i] + df_TargetTouch$chromEnd[i]) / 2,
    y= 0,
    sector.index = df_TargetTouch$chrom[i],
    col= "red",
    pch= 16,
    cex= 1)
}
dev.off()
circos.clear()
```

## ----------------------------------------------
## STEP 2: check expression pattern of genes of interest and filter out low/no-expressing genes.
## ----------------------------------------------

### Load the CRISPRbrain scRNAseq and bulkRNAseq as proxies for the expression profile in NGN2-derived neurons.
```{r}
## scRNAseq:

### Create Seurat objects:
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632022_CRISPRi_filtered_feature_bc_matrix_lib1.h5") -> CRISPRi1
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632023_CRISPRi_filtered_feature_bc_matrix_lib2.h5") -> CRISPRi2
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632025_CRISPRi_filtered_feature_bc_matrix_lib3.h5") -> CRISPRi3
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632026_CRISPRi_filtered_feature_bc_matrix_lib4.h5") -> CRISPRi4
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632027_CRISPRa_filtered_feature_bc_matrix_lib1.h5") -> CRISPRa1
Read10X_h5("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSM4632028_CRISPRa_filtered_feature_bc_matrix_lib2.h5") -> CRISPRa2

CreateSeuratObject(counts=CRISPRi1, project = "CRISPRi1", min.cells = 1,min.features=200) -> CRISPRi1
CreateSeuratObject(counts=CRISPRi2, project = "CRISPRi2", min.cells = 1, min.features=200) -> CRISPRi2
CreateSeuratObject(counts=CRISPRi3, project = "CRISPRi3", min.cells = 1, min.features=200) -> CRISPRi3
CreateSeuratObject(counts=CRISPRi4, project = "CRISPRi4",min.cells = 1, min.features=200) -> CRISPRi4
CreateSeuratObject(counts=CRISPRa1, project = "CRISPRa1",min.cells = 1, min.features=200) -> CRISPRa1
CreateSeuratObject(counts=CRISPRa2, project = "CRISPRa2",min.cells = 1, min.features=200) -> CRISPRa2

### QC:
CRISPRi1$percent.MT <- PercentageFeatureSet(CRISPRi1,pattern="^MT-")
CRISPRi2$percent.MT <- PercentageFeatureSet(CRISPRi2,pattern="^MT-")
CRISPRi3$percent.MT <- PercentageFeatureSet(CRISPRi3,pattern="^MT-")
CRISPRi4$percent.MT <- PercentageFeatureSet(CRISPRi4,pattern="^MT-")
CRISPRa1$percent.MT <- PercentageFeatureSet(CRISPRa1,pattern="^MT-")
CRISPRa2$percent.MT <- PercentageFeatureSet(CRISPRa2,pattern="^MT-")

CRISPRi1filt <- CRISPRi1 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRi2filt <- CRISPRi2 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRi3filt <- CRISPRi3 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRi4filt <- CRISPRi4 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRi4filt <- CRISPRi4 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRa1filt <- CRISPRa1 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)
CRISPRa2filt <- CRISPRa2 %>% subset(nFeature_RNA>200 & nFeature_RNA < 5000 & percent.MT < 10)

### Aggregate counts across cells:
CRISPRi1filt_counts <- data.frame(gene= rownames(CRISPRi1filt@assays$RNA$counts), CRISPRi1= rowSums(CRISPRi1filt@assays$RNA$counts))
CRISPRi2filt_counts <- data.frame(gene= rownames(CRISPRi2filt@assays$RNA$counts), CRISPRi2= rowSums(CRISPRi2filt@assays$RNA$counts))
CRISPRi3filt_counts <- data.frame(gene= rownames(CRISPRi3filt@assays$RNA$counts), CRISPRi3= rowSums(CRISPRi3filt@assays$RNA$counts))
CRISPRi4filt_counts <- data.frame(gene= rownames(CRISPRi4filt@assays$RNA$counts), CRISPRi4= rowSums(CRISPRi4filt@assays$RNA$counts))
CRISPRa1filt_counts <- data.frame(gene= rownames(CRISPRa1filt@assays$RNA$counts), CRISPRa1= rowSums(CRISPRa1filt@assays$RNA$counts))
CRISPRa2filt_counts <- data.frame(gene= rownames(CRISPRa2filt@assays$RNA$counts), CRISPRa2= rowSums(CRISPRa2filt@assays$RNA$counts))

### Merge counts across objects:
sc_dfs <- list(CRISPRi1filt_counts, CRISPRi2filt_counts, CRISPRi3filt_counts, CRISPRi4filt_counts, CRISPRa1filt_counts, CRISPRa2filt_counts)
sc_merged_counts <- Reduce(function(x,y) merge(x, y, by="gene", all= TRUE), sc_dfs)

### Obtain "CPM" values:
sums_cols <- colSums(sc_merged_counts[,-1], na.rm = T)  # need library sizes

sc_merged_counts$CRISPRi1_CPM <- (sc_merged_counts$CRISPRi1 / sums_cols[1]) * 1e6
sc_merged_counts$CRISPRi2_CPM <- (sc_merged_counts$CRISPRi2 / sums_cols[2]) * 1e6
sc_merged_counts$CRISPRi3_CPM <- (sc_merged_counts$CRISPRi3 / sums_cols[3]) * 1e6
sc_merged_counts$CRISPRi4_CPM <- (sc_merged_counts$CRISPRi4 / sums_cols[4]) * 1e6
sc_merged_counts$CRISPRa1_CPM <- (sc_merged_counts$CRISPRa1 / sums_cols[5]) * 1e6
sc_merged_counts$CRISPRa2_CPM <- (sc_merged_counts$CRISPRa2 / sums_cols[6]) * 1e6

sc_merged_counts$mean_cpm <- rowMeans(sc_merged_counts[,8:13],na.rm = T) # get mean CPM per gene across all CRISPRi and CRISPRa libraries
write.csv(sc_merged_counts, file= "Analyses/Kampmann_scRNAseq_summarized.csv", row.names = FALSE)

## bulk RNAseq:

### Load:

bulk_14d1 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neuron\/GSM3543612_D14_C1_read_counts.txt", head=F, sep="")
bulk_14d2 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neurons/GSM3543613_D14_C2_read_counts.txt", head=F, sep="")
bulk_21d1 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neurons/GSM3543614_D21_C1_read_counts.txt", head=F, sep="")
bulk_21d2 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neurons/GSM3543615_D21_C2_read_counts.txt", head=F, sep="")
bulk_28d1 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neurons/GSM3543616_D28_C1_read_counts.txt", head=F, sep="")
bulk_28d2 <- read.delim("finding_cCREs/data/Kampmann_lab/CRISPRbrain/GSE124703_bulkRNAseq_NGN2neurons/GSM3543617_D28_C2_read_counts.txt", head=F, sep="")

bulk_df <- data.frame(gene= bulk_14d1$V1,
                      d14_1= bulk_14d1$V2,
                      d14_2= bulk_14d2$V2,
                      d21_1= bulk_21d1$V2,
                      d21_2= bulk_21d2$V2,
                      d28_1= bulk_28d1$V2,
                      d28_2= bulk_28d2$V2)

### Calculate CPM values (Quant-seq method does not return TPM given it's 3' poly-A restricted)

sums_cols2 <- colSums(bulk_df[,-1])

bulk_df$d14_1_cpm <- (bulk_df$d14_1 / sums_cols2[1]) * 1e6
bulk_df$d14_2_cpm <- (bulk_df$d14_2 / sums_cols2[2]) * 1e6
bulk_df$d21_1_cpm <- (bulk_df$d21_1 / sums_cols2[3]) * 1e6
bulk_df$d21_2_cpm <- (bulk_df$d21_2 / sums_cols2[4]) * 1e6
bulk_df$d28_1_cpm <- (bulk_df$d28_1 / sums_cols2[5]) * 1e6
bulk_df$d28_2_cpm <- (bulk_df$d28_2 / sums_cols2[6]) * 1e6

bulk_df$mean_cpm_14d <- rowMeans(bulk_df[, 8:9])
bulk_df$mean_cpm_21d <- rowMeans(bulk_df[, 10:11])
bulk_df$mean_cpm_28d <- rowMeans(bulk_df[, 12:13])

### Add gene name based on their ENSEMBL ID:

bulk_df_final <- merge(x=bulk_df, y=human_genes, by.x= "gene", by.y= "ENSEMBLID")
write.csv(bulk_df_final, file= "finding_cCREs/results/Kampmann_bulkRNAseq_summarized.csv", row.names = FALSE)
```

### Get basic expression metrics for neuro genes of interest.
```{r}
# Correlation between scRNAseq and bulkRNAseq for the genes of interest:
cor_df <- merge(sc_merged_counts, bulk_df_final, by.x= "gene", by.y="GeneName")

## single-cell vs bulk 14 days:
p14dsc <- ggplot(cor_df, aes(x= log2(mean_cpm+1), y= log2(mean_cpm_14d+1))) + 
  theme_bw() +
  geom_point(color= "gray") +
  geom_smooth(method="lm", color= "black", lty=2) +
  geom_point(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(x= log2(mean_cpm+1), y= log2(mean_cpm_14d+1)), color= "red") +
  xlab("log2(Mean single-cell CPM + 1)") + ylab("log2(Mean 14-days bulk CPM +1)") +
  ggtitle("Bulk 14-days vs single-cell") +
  geom_text_repel(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(label = gene), max.overlaps = 30)
ggsave(p14dsc, filename="finding_cCREs/results/plots/Plot_Bulk14d_singlecell.pdf")

cor.test(log2(cor_df$mean_cpm+1), log2(cor_df$mean_cpm_14d+1)) # r= 0.9242911, p < 2.2x10-16

## single-cell vs bulk 21 days:
p21dsc <- ggplot(cor_df, aes(x= log2(mean_cpm+1), y= log2(mean_cpm_21d+1))) + 
  theme_bw() +
  geom_point(color= "gray") +
  geom_smooth(method="lm", color= "black", lty=2) +
  geom_point(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(x= log2(mean_cpm+1), y= log2(mean_cpm_21d+1)), color= "red") +
  xlab("log2(Mean single-cell CPM + 1)") + ylab("log2(Mean 21-days bulk CPM +1)") +
  ggtitle("Bulk 21-days vs single-cell") +
  geom_text_repel(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(label = gene), max.overlaps = 30)
ggsave(p21dsc, filename="finding_cCREs/results/plots/Plot_Bulk21d_singlecell.pdf")

cor.test(log2(cor_df$mean_cpm+1), log2(cor_df$mean_cpm_21d+1)) # r= 0.9106507, p < 2.2x10-16

## single-cell vs bulk 28 days:
p28dsc <- ggplot(cor_df, aes(x= log2(mean_cpm+1), y= log2(mean_cpm_28d+1))) + 
  theme_bw() +
  geom_point(color= "gray") +
  geom_smooth(method="lm", color= "black", lty=2) +
  geom_point(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(x= log2(mean_cpm+1), y= log2(mean_cpm_28d+1)), color= "red") +
  xlab("log2(Mean single-cell CPM + 1)") + ylab("log2(Mean 28-days bulk CPM +1)") +
  ggtitle("Bulk 28-days vs single-cell") +
  geom_text_repel(data = cor_df[cor_df$gene %in% df_TargetTouch$gene_name,], aes(label = gene), max.overlaps = 30)
ggsave(p28dsc, filename="finding_cCREs/results/plots/Plot_Bulk28d_singlecell.pdf")

cor.test(log2(cor_df$mean_cpm+1), log2(cor_df$mean_cpm_28d+1)) # r= 0.9028248, p < 2.2x10-16
```

### Compare the expression of the neuro genes with the genes selected for CRISPRi and CRISPRa from Kampmann lab.
```{r}
# Load summarized results from the CRISPRbrain webpage:
df_CRISPRi <- read.csv("finding_cCREs/data/Kampmann_lab/Tian2020_CRISPRi_results.csv", head=T, stringsAsFactors = T)
df_CRISPRi$TargetFilter <- as.factor(df_CRISPRi$TargetFilter)
df_CRISPRa <- read.csv("finding_cCREs/data/Kampmann_lab/Tian2020_CRISPRa_results.csv", head=T, stringsAsFactors = T)
df_CRISPRa$TargetFilter <- as.factor(df_CRISPRa$TargetFilter)

# What is the expression of the genes selected for CRISPRi/CRISPRa and the neuro genes?
pall <- ggplot(sc_merged_counts, aes(x= log2(mean_cpm+1))) + 
  theme_bw() +
  geom_histogram(na.rm = T, fill="gray", bins = 100) +
  xlab("Mean log2(CPM+1)") +
  xlim(c(-1,15)) +
  ggtitle("Expression all genes (n= 28,839)")


pCRISPRi <- ggplot(sc_merged_counts[sc_merged_counts$gene %in% df_CRISPRi$TargetGene,], aes(x= log2(mean_cpm+1))) + 
  theme_bw() +
  geom_histogram(na.rm = T, fill="#3B78B0", bins = 100) +
  xlab("Mean log2(CPM+1)") +
  xlim(c(-1,15)) +
  ggtitle("Expression CRISPRi genes (n= 129)")

pCRISPRa <- ggplot(sc_merged_counts[sc_merged_counts$gene %in% df_CRISPRa$TargetGene,], aes(x= log2(mean_cpm+1))) + 
  theme_bw() +
  geom_histogram(na.rm = T, fill="#D1352C", bins = 100) +
  xlab("Mean log2(CPM+1)") +
  xlim(c(-1,15)) +
  ggtitle("Expression CRISPRa genes (n= 97)")

pTT <- ggplot(sc_merged_counts[sc_merged_counts$gene %in% df_TargetTouch$gene_name,], aes(x= log2(mean_cpm+1))) + 
  theme_bw() +
  geom_histogram(na.rm = T, fill="black", bins = 100) +
  xlab("Mean log2(CPM+1)") +
  xlim(c(-1,15)) +
  ggtitle("Expression TargetTouch neuro genes (n= 122)")

(pall + pCRISPRi) / (pCRISPRa + pTT)

ggsave(filename = "finding_cCREs/results/plots/Kampmann_singlecell_expression_summaries.pdf",
       (pall + pCRISPRi) / (pCRISPRa + pTT))


# Do more expression results in best identification of increase (CRISPRa) or decrease (CRISPRi) from a CROP-seq experiment?
df_CRISPR_merged <- rbind(df_CRISPRi,df_CRISPRa)
pCRISPRi_CRISPRa_lm <- ggplot(df_CRISPR_merged[df_CRISPR_merged$TargetFilter=="1",], aes(x= Log2CPM, y= Log2FC)) +
  theme_bw() +
  geom_point(color= "black") +
  geom_smooth(method= "lm", color= "black", lty= 2) +
  geom_point(data= df_CRISPR_merged[df_CRISPR_merged$TargetFilter=="1" & df_CRISPR_merged$FDR < 0.05,], color= "#D1352C") +
  xlab("log2(average CPM)") +
  ylab("log2(fold change)") +
  ggtitle("CRISPR genes (n= 192), significant highlighted") +
  geom_text_repel(aes(label = Gene))

ggsave(filename = "finding_cCREs/results/plots/Kampmann_CRISPRi_CRISPR_all_CPMvsFC.pdf",
       pCRISPRi_CRISPRa_lm)

lm_CRISPR_all <- lm(data=df_CRISPR_merged[df_CRISPR_merged$TargetFilter=="1",], Log2CPM ~ Log2FC)
summary(lm_CRISPR_all) # estimate= -0.25689, p= 1.19e-07

pCRISPRi_lm <- ggplot(df_CRISPRi[df_CRISPRi$TargetFilter=="1",], aes(x= Log2CPM, y= Log2FC)) +
  theme_bw() +
  geom_point(color= "black") +
  geom_smooth(method= "lm", color= "black", lty= 2) +
  geom_point(data= df_CRISPRi[df_CRISPRi$TargetFilter=="1" & df_CRISPRi$FDR < 0.05,], color= "#3B78B0") +
  xlab("log2(average CPM)") +
  ylab("log2(fold change)") +
  ggtitle("CRISPRi genes (n= 129), significant highlighted")

lm_CRISPRi <- lm(data=df_CRISPRi[df_CRISPRi$TargetFilter=="1",], Log2CPM ~ Log2FC)
summary(lm_CRISPRi) # estimate= -0.5202, p= 2.87e-05

pCRISPRa_lm <- ggplot(df_CRISPRa[df_CRISPRa$TargetFilter=="1",], aes(x= Log2CPM, y= Log2FC)) +
  theme_bw() +
  geom_point(color= "black") +
  geom_smooth(method= "lm", color= "black", lty= 2) +
  geom_point(data= df_CRISPRa[df_CRISPRa$TargetFilter=="1" & df_CRISPRa$FDR < 0.05,], color= "#D1352C") +
  xlab("log2(average CPM)") +
  ylab("log2(fold change)") +
  ggtitle("CRISPRa genes (n= 97), significant highlighted")

lm_CRISPRa <- lm(data=df_CRISPRa[df_CRISPRa$TargetFilter=="1",], Log2CPM ~ Log2FC)
summary(lm_CRISPRa) # estimate= -0.26858, p= 0.00159

ggsave(filename = "finding_cCREs/results/plots/Kampmann_CRISPRi_CRISPRa_CPMvsFC.pdf",
       pCRISPRi_lm+pCRISPRa_lm)
```