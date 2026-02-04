# Leverage eQTL and PheWAS databases to highlight cCREs of interest.

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
library(httr)
library(jsonlite)
library(dplyr)
library(readr)
library(LDlinkR)
library(GenomicRanges)

setwd("/home/users/jo912684")
```

## ----------------------------------------------
## STEP 1: Grab eQTL/pQTL databases and overlap with identified cCREs
## ----------------------------------------------

```{r}
# Load eQTL data and overlap with initial set of cCREs:
# Note: due to the large size of the eQTL databases, need to overlap each individually.

### ENCODE eQTL:
gtex_brain_eQTL <- read.delim("finding_cCREs/data/GTEx_Analysis_v10_eQTL_updated/Brain_all_combined_eGenes_v10.txt", head=T , sep="") # 310,349 eQTL hits

gtex_bed <- data.frame(
  chromosome= gtex_brain_eQTL$chr,
  pos= as.numeric(gtex_brain_eQTL$variant_pos),
  eqtl_id= gtex_brain_eQTL$rs_id_dbSNP155_GRCh38p13,
  eqtl_ref= gtex_brain_eQTL$ref,
  eqtl_alt = gtex_brain_eQTL$alt,
  eqtl_gene= gtex_brain_eQTL$gene_name,
  eqtl_beta= as.numeric(gtex_brain_eQTL$slope),
  eqtl_pval = as.numeric(gtex_brain_eQTL$pval_perm),
  eqtl_source = "GTEx",
  eqtl_type= "eQTL"
)
gtex_bed <- gtex_bed[gtex_bed$eqtl_pval < 0.05,] # 126,492 significant hits

### ROSMAP pQTL/eQTL:
rosmap_pqtl <- read.delim("finding_cCREs/data/ROSMAP.for_smr.pQTL v2.txt", head=T, sep="") # 929,695 pQTL hits
rosmap_pqtl$chrom <- paste0("chr", rosmap_pqtl$Chr)

rosmap_pqtl_bed <- data.frame(
  chromosome= rosmap_pqtl$chrom,
  pos= as.numeric(rosmap_pqtl$BP),
  eqtl_id= rosmap_pqtl$SNP,
  eqtl_ref= rosmap_pqtl$A1,
  eqtl_alt = rosmap_pqtl$A2,
  eqtl_gene= rosmap_pqtl$GeneSymbol,
  eqtl_beta= as.numeric(rosmap_pqtl$Beta),
  eqtl_pval = as.numeric(rosmap_pqtl$P),
  eqtl_source = "ROSMAP",
  eqtl_type= "pQTL"
)
rosmap_pqtl_bed <- rosmap_pqtl_bed[rosmap_pqtl_bed$eqtl_pval < 0.05,] # 141,857 significant hits

rosmap_eqtl <- read.csv("finding_cCREs/data/SYNAPSE_eqtl/DLPFC_ROSMAP_cis_eQTL_release.csv", head=T) # 58,766,154 eQTL hits
rosmap_eqtl$chrom <- paste0("chr", rosmap_eqtl$chromosome)

rosmap_eqtl_bed <- data.frame(
  chromosome= rosmap_eqtl$chrom,
  pos= as.numeric(rosmap_eqtl$snpLocation),
  eqtl_id= rosmap_eqtl$snpid,
  eqtl_ref= rosmap_eqtl$A1,
  eqtl_alt = rosmap_eqtl$A2,
  eqtl_gene= rosmap_eqtl$geneSymbol,
  eqtl_beta= as.numeric(rosmap_eqtl$beta),
  eqtl_pval = as.numeric(rosmap_eqtl$FDR),
  eqtl_source = "ROSMAP",
  eqtl_type= "eQTL"
)
rosmap_eqtl_bed <- rosmap_eqtl_bed[rosmap_eqtl_bed$eqtl_pval < 0.05,] # 2,353,264 significant hits

### Mayo cohorts:
mayo_CortexMeta <- read.csv("finding_cCREs/data/SYNAPSE_eqtl/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv", head=T) # 100,080,618 eQTL hits
mayo_CortexMeta$chrom <- paste0("chr", mayo_CortexMeta$chromosome)

mayo_CortexMeta_bed <- data.frame(
  chromosome= mayo_CortexMeta$chrom,
  pos= as.numeric(mayo_CortexMeta$snpLocation),
  eqtl_id= mayo_CortexMeta$snpid,
  eqtl_ref= mayo_CortexMeta$A1,
  eqtl_alt = mayo_CortexMeta$A2,
  eqtl_gene= mayo_CortexMeta$geneSymbol,
  eqtl_beta= as.numeric(mayo_CortexMeta$beta),
  eqtl_pval = as.numeric(mayo_CortexMeta$FDR),
  eqtl_source = "Mayo_CortexMeta",
  eqtl_type= "eQTL"
)

mayo_CortexMeta_bed <- mayo_CortexMeta_bed[mayo_CortexMeta_bed$eqtl_pval < 0.05,] # 4,153,598 significant hits

# Merge all overlapping hits across eQTL/pQTL databases

eqtls_combined <- rbind(gtex_bed, rosmap_pqtl_bed, rosmap_eqtl_bed, mayo_CortexMeta_bed)
eqtls_combined <- eqtls_combined[!duplicated(eqtls_combined),]
eqtls_combined <- eqtls_combined[complete.cases(eqtls_combined),] # 6,562,506 eQTL/pQTL hits

eqtls_combined_gr <- makeGRangesFromDataFrame(eqtls_combined, seqnames.field = "chromosome", start.field="pos", end.field = "pos", keep.extra.columns = T)
  
# Overlap with cCREs:

overlap_cCREs_eQTL <- findOverlaps(initial_overlap_brain_cCREs_gr, eqtls_combined_gr)
overlap_cCREs_eQTL_info <- data.frame(
  CREs= initial_overlap_brain_cCREs_gr[queryHits(overlap_cCREs_eQTL)],
  eQTLs = eqtls_combined_gr[subjectHits(overlap_cCREs_eQTL)]) # 46,151 cCREs have an eQTL or more within them

table(overlap_cCREs_eQTL_info$CREs.TargetTouch.gene == overlap_cCREs_eQTL_info$eQTLs.eqtl_gene) # How many eQTL hits are for the same gene as the predicted brain cCREs?
# FALSE  TRUE 
#  44445   1706 

eqtl_overlap_counts <- data.frame(
  Category = c("TargetGene different to eQTL gene", "TargetGene equal to eQTL gene"),
  Count = table(overlap_cCREs_eQTL_info$CREs.TargetTouch.gene == overlap_cCREs_eQTL_info$eQTLs.eqtl_gene))

p_bar_eQTL <- ggplot(eqtl_overlap_counts, aes(x = 1, y = Count.Freq, fill = Category)) +
  geom_bar(stat = "identity", position="stack") +
  labs(title = "Overlap target genes by brain cCRE and eQTL",
       x = "",
       y = "Number of brain cCREs initial search +/- 100kb from TSS") +
  theme_bw() +
  scale_fill_manual(values= c("gray", "navyblue")) +
  theme(axis.text.x= element_blank(),
        axis.ticks.x = element_blank())
ggsave(p_bar_eQTL, filename= "finding_cCREs/results/plots/Intial_overlap_eQTL_genes.pdf", height=4, width=6)

p_brain_eQTL_genes <- 
  overlap_cCREs_eQTL_info[overlap_cCREs_eQTL_info$CREs.TargetTouch.gene == overlap_cCREs_eQTL_info$eQTLs.eqtl_gene,] %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(CREs.TargetTouch.gene)))) +
  coord_flip() +
  ylab("N brain cCREs with eQTL for the same neuro gene") +
  xlab("")
ggsave(p_brain_eQTL_genes, filename= "finding_cCREs/results/plots/Intial_overlap_eQTL_genes_bar.pdf", width = 8, height = 10)

# Add information to database:
overlap_cCREs_eQTL_info$same_gene <- ifelse(overlap_cCREs_eQTL_info$CREs.TargetTouch.gene == overlap_cCREs_eQTL_info$eQTLs.eqtl_gene, "yes", "no")

# Save:
write.csv(overlap_cCREs_eQTL_info, file="finding_cCREs/results/Overlap_brainCREs_eQTLs.csv")

# Add information to cCREs database:
overlap_cCREs_eQTL_info_gr <- makeGRangesFromDataFrame(overlap_cCREs_eQTL_info, seqnames.field = "CREs.seqnames", start.field = "CREs.start", end.field="CREs.end", keep.extra.columns = T)
overlap_cCREs_eQTL_info_same_gene_gr <- overlap_cCREs_eQTL_info_gr[overlap_cCREs_eQTL_info_gr$same_gene == "yes",]

initial_overlap_brain_cCREs_gr$CRE_contains_eQTL <- "no"
initial_overlap_brain_cCREs_gr$CRE_contains_eQTL[overlapsAny(initial_overlap_brain_cCREs_gr, overlap_cCREs_eQTL_info_gr)] <- "yes"

initial_overlap_brain_cCREs_gr$eQTL_gene_same_as_target_gene <- "no"
initial_overlap_brain_cCREs_gr$eQTL_gene_same_as_target_gene[overlapsAny(initial_overlap_brain_cCREs_gr, overlap_cCREs_eQTL_info_same_gene_gr)] <- "yes"

## Save:
write.csv(as.data.frame(initial_overlap_brain_cCREs_gr), file="finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", row.names = FALSE)

saveRDS(initial_overlap_brain_cCREs_gr, file="finding_cCREs/results/initial_overlap_brain_cCREs_100kb.rds")

```

## ----------------------------------------------
## STEP 2: Grab PheWAS databases and overlap with identified cCREs
## ----------------------------------------------

```{r}
## Load PheWAS data and format into BED:
phewas_data <- read_tsv("finding_cCREs/data/gwas_catalog_v1.0-associations_e113_r2024-11-03.tsv")
phewas_data <- phewas_data[phewas_data$`P-VALUE` < 1e-5, ]
phewas_bed <- data.frame(
  chromosome = paste0("chr",phewas_data$CHR_ID),
  start = as.numeric(phewas_data$CHR_POS) - 1,
  end = as.numeric(phewas_data$CHR_POS),
  SNP = phewas_data$SNPS,
  reported_gene= phewas_data$MAPPED_GENE,
  context= phewas_data$CONTEXT,
  disease_trait= phewas_data[,8]
)
phewas_bed <- phewas_bed[complete.cases(phewas_bed), ] # 593,471 PheWAS hits
phewas_gr <- makeGRangesFromDataFrame(phewas_bed, keep.extra.columns = T, seqnames.field = "chromosome", start.field="start", end.field = "end")

## PheWAS:
phewas_overlaps <- findOverlaps(initial_overlap_brain_cCREs_gr, phewas_gr)
phewas_overlap_info <- data.frame(
  CREs = initial_overlap_brain_cCREs_gr[queryHits(phewas_overlaps)],
  PheWAS = phewas_gr[subjectHits(phewas_overlaps)]
) 
phewas_overlap_info <- phewas_overlap_info[!duplicated(phewas_overlap_info),] # 7,033 brain cCREs from initial search have a PheWAS hit

table(phewas_overlap_info$CREs.TargetTouch.gene == phewas_overlap_info$PheWAS.reported_gene) # How many PheWAS hits are for the same gene as the predicted brain cCREs?
# FALSE  TRUE 
#  5312   1721 

PheWAS_overlap_counts <- data.frame(
  Category = c("TargetGene different to PheWAS gene", "TargetGene equal to PheWAS gene"),
  Count = table(phewas_overlap_info$CREs.TargetTouch.gene == phewas_overlap_info$PheWAS.reported_gene))

p_bar_PheWAS <- ggplot(PheWAS_overlap_counts, aes(x = 1, y = Count.Freq, fill = Category)) +
  geom_bar(stat = "identity", position="stack") +
  labs(title = "Overlap target genes by brain cCRE and PheWAS",
       x = "",
       y = "Number of brain cCREs initial search +/- 100kb from TSS") +
  theme_bw() +
  scale_fill_manual(values= c("gray", "navyblue")) +
  theme(axis.text.x= element_blank(),
        axis.ticks.x = element_blank())
ggsave(p_bar_PheWAS, filename= "finding_cCREs/results/plots/Intial_overlap_PheWAS_genes.pdf", height=4, width=6)

p_brain_PheWAS_genes <- 
  phewas_overlap_info[phewas_overlap_info$CREs.TargetTouch.gene == phewas_overlap_info$PheWAS.reported_gene,] %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(CREs.TargetTouch.gene)))) +
  coord_flip() +
  ylab("N brain cCREs with PheWAS hitting the same neuro gene") +
  xlab("")
ggsave(p_brain_PheWAS_genes, filename= "finding_cCREs/results/plots/Intial_overlap_PheWAS_genes_bar.pdf", , width = 8, height = 10)

p_brain_PheWAS_genes_noAPOE <- 
  phewas_overlap_info[phewas_overlap_info$CREs.TargetTouch.gene == phewas_overlap_info$PheWAS.reported_gene &
                        phewas_overlap_info$CREs.TargetTouch.gene != "APOE",] %>%
  ggplot() + 
  theme_bw() +
  geom_bar(mapping = aes(x = fct_rev(fct_infreq(CREs.TargetTouch.gene)))) +
  coord_flip() +
  ylab("N brain cCREs with PheWAS hitting the same neuro gene") +
  xlab("")
ggsave(p_brain_PheWAS_genes_noAPOE, filename= "finding_cCREs/results/plots/Intial_overlap_PheWAS_genes_bar_noAPOE.pdf", , width = 8, height = 10)

# Add information to database:
phewas_overlap_info$same_gene <- ifelse(phewas_overlap_info$CREs.TargetTouch.gene == phewas_overlap_info$PheWAS.reported_gene, "yes", "no")

# Save:
write.csv(phewas_overlap_info, file="finding_cCREs/results/Overlap_brainCREs_PheWAS.csv", row.names = FALSE)

# Add information to cCREs database:
phewas_overlap_info_gr <- makeGRangesFromDataFrame(phewas_overlap_info, seqnames.field = "CREs.seqnames", start.field = "CREs.start", end.field="CREs.end", keep.extra.columns = T)
phewas_overlap_info_same_gene_gr <- phewas_overlap_info_gr[phewas_overlap_info_gr$same_gene == "yes",]

initial_overlap_brain_cCREs_gr$CRE_contains_PheWAS <- "no"
initial_overlap_brain_cCREs_gr$CRE_contains_PheWAS[overlapsAny(initial_overlap_brain_cCREs_gr, phewas_overlap_info_gr)] <- "yes"

initial_overlap_brain_cCREs_gr$PheWAS_gene_same_as_target_gene <- "no"
initial_overlap_brain_cCREs_gr$PheWAS_gene_same_as_target_gene[overlapsAny(initial_overlap_brain_cCREs_gr, phewas_overlap_info_same_gene_gr)] <- "yes"

## Save:
write.csv(as.data.frame(initial_overlap_brain_cCREs_gr), file="finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", row.names = FALSE)

saveRDS(initial_overlap_brain_cCREs_gr, file="finding_cCREs/results/initial_overlap_brain_cCREs_100kb.rds")
```

## ----------------------------------------------
## STEP 3: Find cCREs with variants in LD with the eQTLs found within cCREs (to expand search)
## ----------------------------------------------

```{r}
eqtls <- read.csv("finding_cCREs/results/Overlap_brainCREs_eQTLs.csv", head=T)
phewas <- read.csv("finding_cCREs/results/Overlap_brainCREs_PheWAS.csv", head=T)

# Keep eQTLs that land outside promoter/exon and in which the eQTL gene = target gene
eqtls <- eqtls[eqtls$CREs.context %in% c("intron", "intergenic", "5UTR", "3UTR") & eqtls$same_gene=="yes",]
phewas <- phewas[phewas$CREs.context %in% c("intron", "intergenic", "5UTR", "3UTR") & phewas$same_gene=="yes",]

# Merge eQTL + PheWAS variants:
all_variants <- unique(eqtls$eQTLs.eqtl_id, phewas$PheWAS.SNP) # 631 variants

# Search in LDlink all variants in LD with our variants:
population <- "ALL"
token <- "c6c7da943a7c" # Obtained here: https://ldlink.nih.gov/?tab=apiaccess

# Function to get LD proxies using LDlinkR
get_ld_proxies <- function(rsid, token, pop) {
  result <- tryCatch({
    LDproxy(snp = rsid,
            pop = population,
            r2d = "r2",
            token = token,
            genome_build = "grch38",
            win_size="500000")
  }, error = function(e) {
    warning(paste("Error with", rsid, ":", e$message))
    return(NULL)
  })
  
  if (is.null(result)) return(NULL)
  
  result <- result %>%
    mutate(eQTL = rsid)
  
  return(result)
}

# Run for all eQTLs
ld_results <- lapply(all_variants, get_ld_proxies, token = token, pop = population)

# Check for the structure of each elements to identify empty ones:
for (i in seq_along(ld_results)) {
  entry <- ld_results[[i]]
  
  if (is.null(entry)) {
    cat("Element", i, ": NULL\n")
  } else if (length(entry) == 0) {
    cat("Element", i, ": EMPTY LIST\n")
  } else if (!is.data.frame(entry)) {
    cat("Element", i, ": NOT A DATA FRAME\n")
  } else {
    dims <- dim(entry)
    cat("Element", i, ": rows =", dims[1], ", cols =", dims[2], "\n")
  }
}

# Manually check the dimensions

# Remove problematic elements due to the eQTL not found in the database
ld_results_cleaned <- ld_results[-c(131, 156, 168, 181, 182, 188, 395, 530, 574, 595, 617)]

# Combine all cleaned data frames into one
ld_df <- do.call(rbind, ld_results_cleaned)

# How many are R2 >= 0.5?
ld_df_filt <- ld_df[ld_df$R2 >= 0.5,]
length(unique(ld_df_filt$RS_Number)) # 8210 LD-variants

# How many of the R2 >= 0.5 land within cCREs?
cCREs_coords <- read.csv("finding_cCREs/results/Initial_overlap_brain_cCREs_neuro_genes_100kb.csv", head=T)

# Step 1: Split Coord into chr and location
coord_split <- strsplit(ld_df_filt$Coord, ":")
ld_df_filt$chr <- sapply(coord_split, `[`, 1)
ld_df_filt$location <- as.integer(sapply(coord_split, `[`, 2))

# Step 2: Create GRanges objects
eqtl_gr <- GRanges(
  seqnames = ld_df_filt$chr,
  ranges = IRanges(start = ld_df_filt$location, end = ld_df_filt$location)
)

ccre_gr <- GRanges(
  seqnames = cCREs_coords$seqnames,
  ranges = IRanges(start = cCREs_coords$start, end = cCREs_coords$end)
)

# Step 3: Find overlaps
hits <- findOverlaps(eqtl_gr, ccre_gr)

# Step 4: Extract overlapping entries from both data frames
ld_df_overlap <- ld_df_filt[queryHits(hits), ]
ccres_overlap <- cCREs_coords[subjectHits(hits), ]

# Step 5: Combine both data frames
ld_df_cre <- cbind(ld_df_overlap, ccres_overlap)

# Save:
write.csv(ld_df_cre, file= "finding_cCREs/results/Overlap_brainCREs_eQTLs_PheWAS_LDvariants.csv")
```