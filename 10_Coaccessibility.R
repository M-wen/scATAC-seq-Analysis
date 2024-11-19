#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="scATAC-seq dataset annotated rds file")
args <- parser$parse_args()
library(ArchR)
addArchRGenome("hg38")
proj <- readRDS(args$input)
proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)
markerGenes  <- c(
    "CD3D","CD4","CD8A",  #T cells
    "CD19", #B cells
    "XBP1", #plasma 
    "CD14", #Monocytes
    "GNLY" # NK cells
  )

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "anno", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(proj)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)