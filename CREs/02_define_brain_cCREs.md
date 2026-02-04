# Combine relevant datasets to define a thorough list of brain cCREs

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
## STEP 1: Grab all ENCODE defined cCREs
## ----------------------------------------------

```{r}
# Load all ENCODE annotated cCREs:
all_CREs_bed <- read.table("finding_cCREs/data/ENCODE/ENCODE_GRCh38-cCREs.bed", head=F, sep="\t", stringsAsFactors = F) # 2,348,854 candidate CREs defined across multiple cell types 
colnames(all_CREs_bed) <- c("chr", "start", "end", "ENCODEID", "ENCODEID2", "classification")
all_CREs_gr <- makeGRangesFromDataFrame(all_CREs_bed, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = T)

# Load files for brain-specific cCREs:
brain_CREs_folder <- "finding_cCREs/data/ENCODE/Brain_cCREs/"
brain_CREs_beds <- list.files(brain_CREs_folder, pattern= "\\.bed$", full.names = TRUE) # data from 16 samples gathered

# Format files into GRanges objects:
brain_CREs_list <- lapply(brain_CREs_beds, function(file){
  read.table(file, header=F, sep="\t", stringsAsFactors = F)
})
brain_CREs_list_gr <- lapply(brain_CREs_list, function(file){
  GRanges(seqnames = file$V1, ranges= IRanges(start= file$V2, end= file$V3))
})

# Combine all brain CREs GRanges objects to obtain a list of shared CREs:
brain_CREs_list_gr_combined <- do.call(c,brain_CREs_list_gr)
brain_CREs_shared <- GenomicRanges::reduce(brain_CREs_list_gr_combined) # 1,056,084 brain cCREs
median(brain_CREs_shared@ranges@width) # median size of brain CREs= 289 bp
write.csv(brain_CREs_shared@ranges, file="finding_cCREs/results/ENCODE_brain_cCREs_shared.csv", row.names = FALSE)

# Overlap the brain cCREs from all cCREs to retain their classification (pELS, dELS, PLS, etc):

overlap_cCREs <- findOverlaps(all_CREs_gr, brain_CREs_shared)
overlap_cCREs_info <- data.frame(
  allCREs = all_CREs_gr[queryHits(overlap_cCREs)],
  brainCREs = brain_CREs_shared[subjectHits(overlap_cCREs)]
) # 1,068,618 brain cCREs

# Remove PLS (promoter-like signatures):
table(overlap_cCREs_info$allCREs.classification)
# CA    CA-CTCF CA-H3K4me3      CA-TF       dELS       pELS        PLS         TF 
# 268      30135      21077         45     797614     176362      42829        288

overlap_cCREs_info_filtered <- overlap_cCREs_info[overlap_cCREs_info$allCREs.classification != "PLS",] # 1,025,789 retained

# Obtain sizes of brain CREs:
p_histbrainCREs <- overlap_cCREs_info_filtered %>%
  ggplot(aes(x= allCREs.width)) +
  theme_bw() +
  geom_histogram(bins = 100, color= "black") +
  geom_vline(xintercept = median(overlap_cCREs_info_filtered$allCREs.width), color= "red", lty= 2) +
  xlab("CRE width") +
  ggtitle("ENCODE - shared Brain cCREs (n= 1,025,789 regions)")

ggsave(p_histbrainCREs, filename= "finding_cCREs/results/plots/ENCODE_brain_cCREs_shared_size.pdf")
```

## ----------------------------------------------
## STEP 2: Use published data for NGN2-neurons for additional cCREs
## ----------------------------------------------

```{r}
# First start with open chromatin regions flagged as accessible via ATAC-seq:
atac_seq <- read.delim("finding_cCREs/data/Nature2025_Human_accelerated_neurons/Human_idr.optimal_peak.narrowPeak.bed", head=F, sep="") # 115,619 regions
atac_seq_gr <- GRanges(seqnames = atac_seq$V1, ranges= IRanges(start= atac_seq$V2, end= atac_seq$V3))
median(atac_seq_gr@ranges@width) # median size of brain CREs= 669 bp

p_histATACseq <- as.data.frame(atac_seq_gr) %>%
  ggplot(aes(x= width)) +
  theme_bw() +
  geom_histogram(bins = 100, color= "black") +
  xlab("ATAC-seq peak width") +
  ggtitle("NGN2-neurons ATAC-seq peaks (n= 115,619)")
ggsave(p_histATACseq, filename= "finding_cCREs/results/plots/NGN2neurons_ATACseq_peaks.pdf")

write.csv(atac_seq_gr@ranges, file="finding_cCREs/results/NGN2neurons_ATACseq_peaks.csv", row.names = FALSE)

# Then, overlap with the two typical marks enriched at enhancer regions: H3Kme1 and H3K27ac:

## Prepare object with all identified H3Kme1 available in ENCODE for brain/neuronal tissues:

### Load files:
H3K4me1_folder <- "finding_cCREs/data/ENCODE/Brain_H3K4me1/"
H3K4me1_folder_bed <- list.files(H3K4me1_folder, pattern= "\\.bed$", full.names = TRUE) # data from 20 samples gathered

### Format files into GRanges objects:
H3K4me1_list <- lapply(H3K4me1_folder_bed, function(file){
  read.table(file, header=F, sep="\t", stringsAsFactors = F)
})
chromosomes_to_retain <- c(paste0("chr",1:22), "chrX", "chrY") # Previous QC revealed that additional sequences were included in the original files (e.g., chrUN)
H3K4me1_list_filtered <- lapply(H3K4me1_list, function(file){
  file <- file[file$V1 %in% chromosomes_to_retain,]
})
H3K4me1_list_gr <- lapply(H3K4me1_list_filtered, function(file){
  GRanges(seqnames = file$V1, ranges= IRanges(start= file$V2, end= file$V3))
})
H3K4me1_list_gr_combined <- do.call(c, H3K4me1_list_gr)

### Combine samples to gather all possible H3K4me1 peaks:
H3K4me1_peaks_shared <- GenomicRanges::reduce(H3K4me1_list_gr_combined) # 409,849 regions
median(H3K4me1_peaks_shared@ranges@width) # median size of H3K4me1 peaks= 649 bp

p_histH3K4me1 <- as.data.frame(H3K4me1_peaks_shared) %>%
  ggplot(aes(x= width)) +
  theme_bw() +
  geom_histogram(bins = 100, color= "black") +
  xlab("H3K4me1 peak width") +
  ggtitle("ENCODE H3K4me1 peaks (n= 409,849)")
ggsave(p_histH3K4me1, filename= "finding_cCREs/results/plots/ENCODE_H3K4me1_peaks.pdf")

write.csv(H3K4me1_peaks_shared@ranges, file="finding_cCREs/results/ENCODE_combined_H3K4me1_peaks.csv", row.names = FALSE)

## Prepare object with all identified H3K27ac available in ENCODE for brain/neuronal tissues:

### Load files:
H3K27ac_folder <- "finding_cCREs/data/ENCODE/Brain_H3K27ac/"
H3K27ac_folder_bed <- list.files(H3K27ac_folder, pattern= "\\.bed$", full.names = TRUE) # data from 38 samples gathered

### Format files into GRanges objects:
H3K27ac_list <- lapply(H3K27ac_folder_bed, function(file){
  read.table(file, header=F, sep="\t", stringsAsFactors = F)
})
H3K27ac_list_filtered <- lapply(H3K27ac_list, function(file){
  file <- file[file$V1 %in% chromosomes_to_retain,]
})
H3K27ac_list_gr <- lapply(H3K27ac_list_filtered, function(file){
  GRanges(seqnames = file$V1, ranges= IRanges(start= file$V2, end= file$V3))
})
H3K27ac_list_gr_combined <- do.call(c, H3K27ac_list_gr)

### Combine samples to gather all possible H3K27ac peaks:
H3K27ac_peaks_shared <- GenomicRanges::reduce(H3K27ac_list_gr_combined) # 308,571 regions
median(H3K27ac_peaks_shared@ranges@width) # median size of H3K27ac peaks= 786 bp

p_histH3K27ac <- as.data.frame(H3K27ac_peaks_shared) %>%
  ggplot(aes(x= width)) +
  theme_bw() +
  geom_histogram(bins = 100, color= "black") +
  xlab("H3Kme1 peak width") +
  ggtitle("ENCODE H3K27ac peaks (n= 308,571)")
ggsave(p_histH3K27ac, filename= "Analyses/Plots/ENCODE_H3K27ac_peaks.pdf")

write.csv(H3K27ac_peaks_shared@ranges, file="finding_cCREs/results/ENCODE_combined_H3K27ac_peaks.csv", row.names = FALSE)

# Now, overlap the chromatin marks (H3K4me1 & H3K27ac) and open chromatin marks (ATAC-seq) to find genomic locations flagged as potential enhancers:

overlap_enhancer_marks_atac <- Reduce(subsetByOverlaps, list(H3K4me1_list_gr_combined, H3K27ac_list_gr_combined, atac_seq_gr))
overlap_enhancer_marks_atac_df <- as.data.frame(overlap_enhancer_marks_atac)
overlap_enhancer_marks_atac_df <- overlap_enhancer_marks_atac_df[overlap_enhancer_marks_atac_df$width <= 1000,] # drop sites larger than 1000 bp since they're likely due to errors in the definition of peaks
# 226,692 peaks

```

## ----------------------------------------------
## STEP 3: Additional STARR-seq data for NGN2-neurons
## ----------------------------------------------

```{r}
### Load files:
STARR_folder <- "finding_cCREs/data/ENCODE/NGN2neurons_STARRseq/"
STARR_folder_bed <- list.files(STARR_folder, pattern= "\\.bed$", full.names = TRUE) # data from 4 samples gathered

### Format files into GRanges objects:
STARR_list <- lapply(STARR_folder_bed, function(file){
  read.table(file, header=F, sep="\t", stringsAsFactors = F)
})
chromosomes_to_retain <- c(paste0("chr",1:22), "chrX", "chrY")
STARR_list_filtered <- lapply(STARR_list, function(file){
  file <- file[file$V1 %in% chromosomes_to_retain,]
})
STARR_list_gr <- lapply(STARR_list_filtered, function(file){
  GRanges(seqnames = file$V1, ranges= IRanges(start= file$V2, end= file$V3))
})
STARR_list_gr_combined <- do.call(c, STARR_list_gr)
STARR_list_gr_combined_df <- as.data.frame(STARR_list_gr_combined) # 1,797 STARR-seq peaks
```

## ----------------------------------------------
## STEP 4: Gather candidate silencing regions
## ----------------------------------------------

Overlap chromatin accesibility information with available data on transcription factor REST. Despite not having chromatin marks for silencers, studies have found enrichment for REST transcription factor.

```{r}
# REST Chip-seq data from neuronal cell lines

## Load files:
REST_folder <- "finding_cCREs/data/ENCODE/Brain_REST/"
REST_bed <- list.files(REST_folder, pattern= "\\.bed$", full.names = TRUE) # data from 4 samples gathered

## Format files into GRanges objects:
REST_list <- lapply(REST_bed, function(file){
  read.table(file, header=F, sep="\t", stringsAsFactors = F)
})
REST_list_gr <- lapply(REST_list, function(file){
  GRanges(seqnames = file$V1, ranges= IRanges(start= file$V2, end= file$V3))
})
REST_list_gr_combined <- do.call(c, REST_list_gr)

### Combine samples to gather all possible REST peaks:
REST_list_gr_shared <- GenomicRanges::reduce(REST_list_gr_combined) # 38,508 REST peaks
REST_list_gr_shared_df <- as.data.frame(REST_list_gr_shared)

# Load silencer databases available

## Candidate silencers from Ss-STARR-seq:
sil_db_STARRseq1 <- read.delim("finding_cCREs/data/NatComm2025_Uncovering whole genome silencers using SsSTARRseq/GSE283488_K562_silencer.bed", head=T, sep = "")
sil_db_STARRseq2 <- read.delim("finding_cCREs/data/NatComm2025_Uncovering whole genome silencers using SsSTARRseq/GSE283488_293T_silencer.bed", head=T, sep="")
sil_db_STARRseq3 <- sil_db_STARRseq3 <- read.delim("finding_cCREs/data/NatComm2025_Uncovering whole genome silencers using SsSTARRseq/GSE283488_LNcap_silencer.bed", head=T, sep="")
sil_db_STARRseq_combined <- rbind(sil_db_STARRseq1, sil_db_STARRseq2, sil_db_STARRseq3) # 397,266 total silencers predicted

## Candidate silencers from Jayavelu 2020:
sil_2020 <- read.delim("finding_cCREs/data/NatComm2020_Candidate_silencer_elements_human_mouse_genomes/Candidate_silencers_and_uncharacterized_CREs_human_hg19_ENCODE_cell_types.txt", head=T, sep="\t")
sil_2020 <- sil_2020[sil_2020$SVM_prediction == "Candidate_Silencer",] # 661,417 total silencers predicted

## Human SilencerDB:
sil_db <- read.delim("finding_cCREs/data/SilencerDB/Homo_sapiens.bed", head=F, sep="\t")
colnames(sil_db) <- c("chr", "start", "end", "silencerID", "size", "strand", "Celtype", "Tissue", "TissueGeneral", "Method")
sil_db <- sil_db[sil_db$TissueGeneral=="Brain",] # 206,514 total silencers predicted

# Combine all predicted silencer regions

## Add an identified of the database the prediction comes from:
colnames(REST_list_gr_shared_df) <- c("chr", "start", "end", "width", "strand", "source")
REST_list_gr_shared_df$source <- "ENCODE_REST"
sil_db_STARRseq_combined$source <- "Ss-STARR-seq"
sil_2020$source <- "Jayavelu2020"
sil_db$source <- "SilencerDB"

silencers_predictions <- rbind(REST_list_gr_shared_df[,c(1:3,6)],
                              sil_db_STARRseq_combined[,c(1:3,5)],
                               sil_2020[,c(1:3,9)],
                               sil_db[,c(1:3,11)])

silencers_predictions$start <- as.numeric(silencers_predictions$start)
silencers_predictions$end <- as.numeric(silencers_predictions$end)
silencers_predictions <- silencers_predictions[complete.cases(silencers_predictions),]
silencers_predictions$width <- silencers_predictions$end - silencers_predictions$start
silencers_predictions <- silencers_predictions[silencers_predictions$width > 0 & silencers_predictions$width <= 1000,]
```

## ----------------------------------------------
## STEP 5: Put together a thorough database of potential neuronal cCREs (enhancers + silencers)
## ----------------------------------------------

```{r}
# Enhancers:

## Combine all regions into a large dataframe containing all brain/NGN2 cCREs:

# Need to add a 'classification' column to keep consistency with ENCODE data:
overlap_enhancer_marks_atac_df$classification <- "ATAC_H3K4me1_H3K27ac"
STARR_list_gr_combined_df$classification <- "STARRseq"
ENCODE_brain_cCREs <- data.frame(seqnames= overlap_cCREs_info$allCREs.seqnames,
                                 start=overlap_cCREs_info$allCREs.start,
                                 end= overlap_cCREs_info$allCREs.end,
                                 width= overlap_cCREs_info$allCREs.width,
                                 strand= overlap_cCREs_info$allCREs.strand,
                                 classification= overlap_cCREs_info$allCREs.classification)

final_brain_enhancers <- rbind(ENCODE_brain_cCREs,
                           overlap_enhancer_marks_atac_df,
                           STARR_list_gr_combined_df)
final_brain_enhancers$cCRE_type <- "enhancer"

# Silencers:
silencers_predictions$cCRE_type <- "silencer"
final_brain_silencers <- data.frame(seqnames= silencers_predictions$chr,
                                    start= silencers_predictions$start,
                                    end= silencers_predictions$end,
                                    width= silencers_predictions$width,
                                    strand= "*",
                                    classification= silencers_predictions$source,
                                    cCRE_type= silencers_predictions$cCRE_type)

# Final list of cCREs:
final_brain_cCREs <- rbind(final_brain_enhancers, final_brain_silencers) # 2,589,508 cCREs
# enhancer silencer 
# 1297107  1292401 
  
median(final_brain_cCREs$width) # median size of CREs= 192 bp

p_brainCREs_box <- ggplot(final_brain_cCREs, aes(x=classification, y=width, fill=classification)) +
  theme_bw() +
  geom_boxplot() +
  ylab("CRE width") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
ggsave(p_brainCREs_box, filename= "finding_cCREs/results/plots/Final_brain_cCREs_boxplot.pdf", width = 10, height=6)

p_brainCREs_hist <- ggplot(final_brain_cCREs, aes(x= width)) +
  theme_bw() +
  geom_histogram(bins = 100, color= "black") +
  xlab("CRE width") +
  ggtitle("Final list of neuronal cCREs (n= 2,589,508)")
ggsave(p_brainCREs_hist, filename= "finding_cCREs/results/plots/Final_brain_cCREs_hist.pdf")

final_brain_CREs_gr <- makeGRangesFromDataFrame(final_brain_cCREs, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# Save list of genome-wide brain/neuron cCREs:
write.csv(final_brain_cCREs, file="finding_cCREs/results/Final_brain_cCREs.csv", row.names = FALSE)
saveRDS(final_brain_CREs_gr, file="finding_cCREs/results/final_brain_CREs_gr.rds")
```