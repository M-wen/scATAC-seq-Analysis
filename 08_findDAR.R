#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="scATAC-seq dataset annotated and add peakMatrix rds file")

args <- parser$parse_args()

library(ArchR)
library(data.table)

proj <- readRDS(args$input)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "predictedsubGroup",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
