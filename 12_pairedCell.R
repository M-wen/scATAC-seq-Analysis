#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="integrated scATAC-seq with scRNA-seq rds file")
parser$add_argument("-o", "--outpath", required=TRUE, help="out put path")

args <- parser$parse_args()

library(Seurat)
library(FigR)
library(optmatch)
library(dplyr)
# 01 using FigR pairing scRNA and snATAC
## cca
# object <- readRDS("${project_path}/results/06_integrated/pbmc_scRNA_scATAC_integrated_harmony.rds")
object <- readRDS(args$input)
CCA_PCs <- object@reductions$harmony@cell.embeddings
isATAC <- grepl("#",rownames(CCA_PCs))
table(isATAC) # ATAC vs RNA

ATACcells <- rownames(CCA_PCs)[isATAC]
RNAcells <- rownames(CCA_PCs)[!isATAC]

ATAC_PCs <- CCA_PCs[isATAC,]
RNA_PCs <- CCA_PCs[!isATAC,]

pairing <- pairCells(ATAC = ATAC_PCs,
                     RNA = RNA_PCs,
                     keepUnique = TRUE
                     )
## deduplicate
if (length(unique(pairing$RNA)) > length(unique(pairing$ATAC))) {
    dupli_ATAC <- unique(pairing$ATAC[duplicated(pairing$ATAC)])
    dupli_pairing <- pairing[which(pairing$ATAC %in% dupli_ATAC),]
    dupli_pairing <- data.frame(ATAC = dupli_pairing$ATAC,RNA = dupli_pairing$RNA, dist = dupli_pairing$dist)
    ## deduplicate
    dedupli_pairing <- dupli_pairing %>% group_by(ATAC) %>% top_n(n=-1, wt=dist)
    unique_pairing <- pairing[which(! pairing$ATAC %in% dupli_ATAC),]
    total_uniq_pairing <- rbind(unique_pairing,dedupli_pairing)
}else {
    dupli_RNA <- unique(pairing$RNA[duplicated(pairing$RNA)])
    dupli_pairing <- pairing[which(pairing$RNA %in% dupli_RNA),]
    dupli_pairing <- data.frame(ATAC = dupli_pairing$ATAC,RNA = dupli_pairing$RNA, dist = dupli_pairing$dist)
    ## deduplicate
    dedupli_pairing <- dupli_pairing %>% group_by(RNA) %>% top_n(n=-1, wt=dist)
    unique_pairing <- pairing[which(! pairing$RNA %in% dupli_RNA),]
    total_uniq_pairing <- rbind(unique_pairing,dedupli_pairing)
}

paired <- data.frame(ATAC = total_uniq_pairing$ATAC, RNA = total_uniq_pairing$RNA, muti = paste0("multi_cell","_",c(1:length(total_uniq_pairing$ATAC))))
saveRDS(paired,paste0(args$outpath,"/results/06_integrated/paired.rds"))
