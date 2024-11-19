#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-a", "--ATAC", required=TRUE, help="scATAC-seq dataset annotated and add peakMatrix rds file")
parser$add_argument("-r", "--RNA", required=TRUE, help="scRNA-seq dataset annotated rds file")
parser$add_argument("-o", "--outpath", required=TRUE, help="out put path")

args <- parser$parse_args()
library(ArchR)
library(Seurat)
library(Signac)
library(future)
library(future.apply)
library(parallel)
library(dplyr)
library(pbmc3k.SeuratData)
# data("pbmc3k")
# pbmc_rna <- pbmc3k
# seRNA@active.assay <- "RNA"
pbmc_rna <- readRDS(args$RNA)
pbmc_rna <- NormalizeData(pbmc_rna)
pbmc_rna <- FindVariableFeatures(pbmc_rna,nfeatures = 3000)
pbmc_rna <- ScaleData(pbmc_rna)
pbmc_rna <- RunPCA(pbmc_rna)
pbmc_rna <- FindNeighbors(pbmc_rna, dims = 1:30)
pbmc_rna <- FindClusters(pbmc_rna)
genesUse <- VariableFeatures(object = pbmc_rna)
pbmc_atac <- readRDS(args$ATAC)
allCells <- pbmc_atac$cellNames
useMatrix <- "GeneScoreMatrix"
geneDF <- ArchR:::.getFeatureDF(getArrowFiles(pbmc_atac), useMatrix)
GeneScoreMatrix <- ArchR:::.getPartialMatrix(
  ArrowFiles = getArrowFiles(pbmc_atac),
  featureDF = geneDF[geneDF$name %in% genesUse,],
  threads = 10,
  cellNames = allCells,
  useMatrix = useMatrix,
  verbose = FALSE
)
rownames(GeneScoreMatrix) <- geneDF[geneDF$name %in% genesUse, "name"]
mat <- log(GeneScoreMatrix + 1)
Seurat_ATAC_pbmc <- CreateSeuratObject(
        counts = GeneScoreMatrix,
        assay = 'GeneScore',
        project = 'ATAC',
        min.cells = 1,
        meta.data = as.data.frame(pbmc_atac@cellColData)
)
rm(list=c("mat","GeneScoreMatrix"))
gc()
Seurat_ATAC_pbmc <- ScaleData(Seurat_ATAC_pbmc, verbose = FALSE)
DefaultAssay(Seurat_ATAC_pbmc) <- 'GeneScore'
Seurat_ATAC_pbmc$Tech <- "ATAC"
pbmc_rna$Tech <- "RNA"
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024^3)

DefaultAssay(pbmc_rna) <- 'RNA'
reference.list <- c(pbmc_rna, Seurat_ATAC_pbmc)
names(reference.list) <- c("RNA", "ATAC")
rna_atac.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
rna_atac_integrated <- IntegrateData(anchorset = rna_atac.anchors, dims = 1:30)
rna_atac_integrated <- ScaleData(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- RunPCA(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- FindNeighbors(object = rna_atac_integrated, dims = 1:30)
rna_atac_integrated <- FindClusters(object = rna_atac_integrated, resolution = 0.5)
library(harmony)
rna_atac_integrated <- RunHarmony(rna_atac_integrated, "Tech")
rna_atac_integrated <- RunUMAP(rna_atac_integrated, reduction =  "harmony",reduction.name = "UMAPHarmony",dims = 1:30)
rna_atac_integrated$Merged_cluster <- ifelse(rna_atac_integrated@meta.data$Tech == "RNA", 
na_atac_integrated@meta.data$seurat_annotations, rna_atac_integrated@meta.data$predictedsubGroup)
saveRDS(rna_atac_integrated, paste0(args$outpath,"/results/06_integrated/pbmc_scRNA_scATAC_integrated_harmony.rds"))
