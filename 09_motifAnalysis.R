#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="scATAC-seq dataset annotated rds file")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

library(ArchR)
set.seed(10)
addArchRGenome("hg38")

proj <- readRDS(args$input)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
 ArchRProj = proj,
 peakAnnotation = "Motif",
 force = TRUE
)
saveRDS(proj,psate0(args$dir,"/",args$output,"_addMotifMatrix.rds"))

