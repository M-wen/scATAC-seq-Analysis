#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-p", "--preprocess", required=TRUE, help="preprocessed dir")
parser$add_argument("-s", "--singlet", required=TRUE, help="singlet dir")
parser$add_argument("-g", "--genome", default = "hg38", help="used genome [default %(default)s]")
parser$add_argument("--blacklist", help="blacklist in bed gz")
parser$add_argument("--chromSize", help="chrom size")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")

# load library
library(ArchR)
library(ggplot.multistats)
library(readxl)
library(scuttle)
library(cowplot)

addArchRGenome(args$genome)

# 1. make arrow files
raw <- args$preprocess #insert path to raw data
amulet <- args$singlet
ArrowDir <- paste0(args$dir,"/ArrowFiles")
cmd = paste("mkdir -p",ArrowDir)
system(cmd)
setwd(ArrowDir)

samples <- list.files(raw)

getValidBarcodes <- function (csvFiles = NULL, multiplet = NULL, sampleNames = NULL)
{
   ArchR:::.validInput(input = csvFiles, name = "csvFiles", valid = "character")
   ArchR:::.validInput(input = multiplet, name = "multiplet", valid = "character")
   ArchR:::.validInput(input = sampleNames, name = "sampleNames", valid = "character")
   if (length(sampleNames) != length(csvFiles)) {
       stop("csvFiles and sampleNames must exist!")
   }
   if (!all(file.exists(csvFiles))) {
       stop("Not All csvFiles exists!")
   }
     if (!all(file.exists(multiplet))) {
       stop("Not All multipletFiles exists!")
   }
   barcodeList <- lapply(seq_along(csvFiles), function(x) {
       df <- ArchR:::.suppressAll(data.frame(readr::read_csv(csvFiles[x])))
       multi <- read.table(multiplet[x],header=F)
       if ("is__cell_barcode" %ni% colnames(df)) {
           stop("is__cell_barcode not in colnames of 10x singlecell.csv file! Are you sure inut is correct?")
       }
       as.character(df[which(paste0(df$is__cell_barcode) != 0 & !df$barcode %in% multi$V1),
           ]$barcode)
   }) %>% SimpleList
   names(barcodeList) <- sampleNames
   barcodeList
}

samples <- c("10k_pbmc_ATACv1-1_nextgem_Chromium_X","10k_pbmc_ATACv2_nextgem_Chromium_Controller","10k_pbmc_ATACv2_nextgem_Chromium_X")
for(i in 1:length(samples)){

   # use only valid 10x barcodes
   barcodes <- getValidBarcodes(csvFiles = paste0(raw, samples[i],
       "/outs/singlecell.csv"), multiplet = paste0(amulet,samples[i],"/MultipletBarcodes_01.txt"),sampleNames =
       samples[i])

   ArrowFiles <- createArrowFiles(inputFiles = paste0(raw,
       samples[i], "/outs/fragments.tsv.gz"),
       sampleNames = samples[i], minTSS = 0, minFrags = 0,
       validBarcodes = barcodes[[1]], addTileMat = TRUE,
       addGeneScoreMat = TRUE, force=TRUE, offsetPlus = 0, offsetMinus = 0,
       excludeChr = c("chrM", "chrY", "chrX"))
}

# 2. find doublets

ArrowFiles <- paste0(ArrowDir, samples, ".arrow")

set.seed(1)
doubScores <- addDoubletScores(input = ArrowFiles, k = 20,
    knnMethod = "UMAP", LSIMethod = 1, force=TRUE)

# 3. make ArchR project
setwd(args$dir)

sc <- ArchRProject(ArrowFiles = ArrowFiles,
    outputDirectory = args$dir, copyArrows = FALSE)

# 4. quality control plots

metadata <- as.data.frame(sc@cellColData)
colours <- ArchR:::paletteDiscrete( values = unique(metadata$Sample))
system("mkdir ./Plots")

metadata$Sample <- factor(metadata$Sample, levels = unique(metadata$Sample))
P1<- ggplot(metadata, aes(x = Sample, y = TSSEnrichment, fill=Sample)) +
  geom_violin(trim=TRUE,color="white",show.legend = F) +
  geom_boxplot(width=0.1,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,outlier.stroke = 0)+
  scale_fill_manual(values = colours)+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank())+
  ylab("TSSEnrichment")

P2<- ggplot(metadata, aes(x = Sample, y = log10(nFrags), fill=Sample)) +
  geom_violin(trim=TRUE,color="white",show.legend = F) +
  geom_boxplot(width=0.1,position=position_dodge(0.9),show.legend = F,fill="white",outlier.size = 0,outlier.stroke = 0)+
  scale_fill_manual(values = colours)+
  theme_cowplot()+
  scale_y_continuous(limits = c(3, 5))+
  theme(axis.text.x=element_text(colour="black",family="Times",size=10))+
  ylab(expression(Log[10]*paste("(","nFragment",")",sep = "")))+xlab("")
ggsave("./Plots/sample_qc.pdf",do.call(plot_grid,c(list(P1,P2),ncol=1,align = "v")),width=3,height=3)

gg <- ArchR:::ggPoint(
  x = pmin(log10(metadata$nFrags), 5) + rnorm(length(metadata$nFrags), sd = 0.00001),
  y = metadata$TSSEnrichment + rnorm(length(metadata$nFrags), sd = 0.00001), 
  colorDensity = TRUE,
  xlim = c(2.5, 5),
  ylim = c(0, max(metadata$TSSEnrichment) * 0.8),
  baseSize = 6,
  continuousSet = "sambaNight",
  xlabel = "Log 10 (Unique Fragments)",
  ylabel = "TSS Enrichment",
  rastr = TRUE) + 
  geom_hline(yintercept=4, lty = "dashed", size = 0.25) +
  geom_vline(xintercept=log10(1000), lty = "dashed", size = 0.25)
ggsave("./Plots/totalcell_qc.pdf",gg,width=5,height=5)

# 5. filter doublets

sc <- filterDoublets(sc)

#6. filter cells
tss_outliers <- list()
tss_outliers_names <- list()

for(i in 1:length(samples)){

    sample_i <- samples[i]
    tss_enrich <- sc$TSSEnrichment[sc$Sample==sample_i]
    tss_outliers[[sample_i]] <- isOutlier(tss_enrich, nmads=1,  type="lower")
    tss_outliers_names[[sample_i]] <-
        rownames(sc@cellColData)[sc$Sample==sample_i][tss_outliers[[sample_i]]]
}

lower_threshold <- sapply(tss_outliers, function(x) attr(x, "thresholds")[1])
lower_threshold
filter_tss <- sapply(tss_outliers, function(x) sum(x))
filter_tss
tss_outliers_names <- unlist(tss_outliers_names)
sc <- sc[!rownames(sc@cellColData) %in% tss_outliers_names,]
sc <- sc[which(log10(sc$nFrags) > 3.5 ),]

saveRDS(sc,
    file="./scatac_filtered.rds")

# 7. dimension reduction

sc <- addIterativeLSI(ArchRProj = sc, useMatrix = "TileMatrix",
    name = "IterativeLSI", iterations = 4, clusterParams = list(
    resolution = 4, sampleCells = 2000, n.start = 10),
    varFeatures = 50000, dimsToUse = 1:30, force=TRUE, seed=1)
# 8. Assess the correlation between each LSI component and sequencing depth
embed <- sc@reducedDims$IterativeLSI$matSVD
counts <- subset(sc@cellColData,select= "nFrags")
embed <- embed[rownames(x = counts), ]

n = 10 
n <- Signac:::SetIfNull(x = n, y = ncol(x = embed))
embed <- embed[, seq_len(length.out = n)]
counts$nFrags <- as.numeric(counts$nFrags)

depth.cor <- as.data.frame(cor(x = embed, y = counts$nFrags))
depth.cor$Component <- rownames(depth.cor)
colnames(depth.cor)[1] <- "correaltion"
depth.cor$Component <- factor(depth.cor$Component, levels = depth.cor$Component)
p <- ggplot(depth.cor, aes(Component, correaltion)) +
    geom_point() +
    # scale_x_continuous(n.breaks = 10, limits = c(1, 10)) +
    ylab("Correlation") +
    ylim(c(-1, 1)) +
    theme_light() +
    ggtitle("Correlation between depth and reduced dimension components")
ggsave("./Plots/correaltion_with_depth.pdf",p,width=5,height=5)

# 9. add batch correction

sc <- addHarmony(ArchRProj = sc, reducedDims = "IterativeLSI",
    name = "Harmony", groupBy = "Sample", force=TRUE,dimsToUse = 2:30)

# 10. clustering

sc <- addClusters(input = sc, reducedDims = "Harmony", method = "Seurat",
    name = "Clusters", resolution = 0.8, force=TRUE, seed = 1,
    maxClusters = 100,dimsToUse = 2:30)

# 11. UMAP

sc <- addUMAP(ArchRProj = sc, reducedDims = "Harmony", name = "UMAP",
    nNeighbors = 30, minDist = 0.5, metric = "cosine",
    force=TRUE, seed = 1,dimsToUse = 2:30)

p1 <- plotEmbedding(ArchRProj = sc, colorBy = "cellColData", name = "Clusters",
    embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = sc, colorBy = "cellColData", name = "Sample",
    embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = sc,
    addDOC = FALSE, width = 5, height = 5)
saveRDS(sc,
    file="./scatac_final.rds")

# 12. imputation with MAGIC
sc <- addImputeWeights(
  ArchRProj = sc,
  reducedDims = "IterativeLSI",
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  td = 3,
  ka = 4,
  sampleCells = 1000,
  nRep = 2,
  k = 15,
  epsilon = 1,
  useHdf5 = TRUE,
  randomSuffix = FALSE,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  logFile = createLogFile("addImputeWeights")
)

# 13. annotation using marker genes
markerGenes <- c(
    "CD3D","CD4","CD8A",  #T cells
    "CD19", #B cells
    "XBP1", #plasma 
    "CD14", #Monocytes
    "GNLY" # NK cells
  )
p <- plotEmbedding(
        ArchRProj = sc,
        colorBy = "GeneScoreMatrix",
        name = markerGenes,
        embedding = "UMAP",
        imputeWeights = getImputeWeights(sc),
        plotAs = 'points'
    )
P5 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
ggsave("./Plots/MarkerGenes_Feartureplot.pdf",do.call(cowplot::plot_grid,c(list(ncol = 3),P5)))

library(pbmc3k.SeuratData)

seRNA <- data("pbmc3k")
seRNA@active.assay <- "RNA"
sc <- addGeneIntegrationMatrix(
    ArchRProj = sc,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "seurat_annotations",
    nameCell = "predictedsubCell",
    nameGroup = "predictedsubGroup",
    nameScore = "predictedsubScore"
)
P6 <- plotEmbedding(
    proj,
    colorBy = "cellColData",
    name = "predictedsubGroup",
    embedding = "UMAP"
 )
plotPDF(P6, name = "Plot-UMAP-CellTypeAnnotation.pdf", ArchRProj = sc,
addDOC = FALSE, width = 5, height = 5)
cM <- confusionMatrix(sc$Clusters, sc$predictedsubGroup)
cM
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)
p7 <- pheatmap::pheatmap(mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), border_color = "black")

ggsave("./Plots/AnnotationHeatmapPlot.pdf",p7,width=5,height=5)
sc$annotation <- mapLabs[match(sc$seurat_clusters, mapLabs[,1]),2]
colours <- ArchR:::paletteDiscrete( values = unique(sc@meta.data$annotation))

p8 <- DimPlot(sc, reduction = "umap",group.by="annotation",cols = colours) +
     theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave("./Plots/AnnotationHeatmapPlot.pdf",p8,width=5,height=5)

saveRDS(sc, "pbmc_scATAC_annotation.rds")

