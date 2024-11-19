
#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-a", "--ATAC", required=TRUE, help="scATAC-seq dataset annotated and add peakMatrix rds file")
parser$add_argument("-p", "--paired", required=TRUE, help="scATAC and scRNA paired rds file")
parser$add_argument("-i", "--input", required=TRUE, help="integrated scATAC-seq with scRNA-seq rds file")
parser$add_argument("-o", "--outpath", required=TRUE, help="out put path")

args <- parser$parse_args()

library(tidyverse)
library(Signac)
library(data.table)
library(GenomicRanges)
library(Seurat)
library(data.table)
library(Pando)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ArchR)
addArchRGenome("hg38")
# proj <- >readRDS("${project_path}/results/04_parsePeak/pbmc_scATAC_addpeakMatrix.rds")
proj <- >readRDS(args$ATAC)
useMatrix <- "PeakMatrix" 
# paired <- readRDS("${project_path}/06_integrated/paired.rds")
paired <- readRDS(args$paired)
peakDF <- ArchR:::.getFeatureDF(getArrowFiles(proj), useMatrix)
peaks <- paste(peakDF$seqnames, peakDF$start,peakDF$end,sep="_")
cells <- paired$ATAC
PeakMatrix <- ArchR:::.getPartialMatrix(
  ArrowFiles = getArrowFiles(proj),
  featureDF = peakDF,
  threads = 10,
  cellNames = cells,
  useMatrix = useMatrix,
  verbose = FALSE
) 
rownames(PeakMatrix) <- peaks
geneAnnotation <- geneAnnoHg38
genes <- geneAnnotation$genes
exons <-  geneAnnotation$exons
genesdf <- data.frame(seqnames = seqnames(genes),start = start(genes),end = end(genes),strand = strand(genes),gene_id = genes$gene_id,symbol = genes$symbol,gene_name = genes$symbol,type = "gene")
exonsdf <- data.frame(seqnames = seqnames(exons),start = start(exons),end = end(exons),strand = strand(exons),gene_id = exons$gene_id,symbol = exons$symbol,gene_name = exons$symbol,type = "exon")
gene_annot_df <- rbind(genesdf,exonsdf)
gene_annot <- makeGRangesFromDataFrame(gene_annot_df,keep.extra.columns= T)
peakassay <- CreateChromatinAssay(PeakMatrix, sep=c("_","_"),annotation = gene_annot)
ATAC <- CreateSeuratObject(
 counts = peakassay,
 assay = "peaks"
)
ATAC <- ATAC[,paired$ATAC]
ATAC <- RenameCells(ATAC,old.names = colnames(ATAC),new.names = paired$muti)
peakassay <- CreateChromatinAssay(PeakMatrix, sep=c("_","_"),annotation = gene_annot)
ATAC <- CreateSeuratObject(
 counts = peakassay,
 assay = "peaks"
)
ATAC <- ATAC[,paired$ATAC]
ATAC <- RenameCells(ATAC,old.names = colnames(ATAC),new.names = paired$muti)
coembed <- readRDS(args$input)
RNA <- coembed[,paired$RNA]
RNA <- RenameCells(RNA,old.names = colnames(RNA),new.names = paired$muti)
multi <- RNA 
multi[["peaks"]] <- ATAC[["peaks"]]
saveRDS(multi,paste0(args$output,"/06_integrated/pbmc_multiomic_object.rds"))

