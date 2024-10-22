#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="scATAC-seq dataset annotated rds file")
parser$add_argument("-p", "--peakset", required=TRUE,help="union peakset file")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

ibrary(ArchR)
library(data.table)
set.seed(10)
addArchRGenome("hg38")
proj <- readRDS(args$input)
peakset <- as.data.frame(fread(args$peakset,header=T)),header=T))
peakSet_gr <- GRanges(
    peakset$seqnames,
    IRanges(peakset$start, peakset$end)
  )
mcols(peakSet_gr)$score <- peakset$score
mcols(peakSet_gr)$name <- peakset$name
mcols(peakSet_gr)$Group <- peakset$label
mcols(peakSet_gr)$spm <- peakset$spm
genomeAnnotation = getGenomeAnnotation(proj)
BSgenome <- eval(parse(text = genomeAnnotation$genome))
BSgenome <- validBSgenome(BSgenome)
geneAnnotation = getGeneAnnotation(proj)
promoterRegion = c(2000, 100)
peakSet_gr <- ArchR:::.fastAnnoPeaks(peakSet_gr, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion)
proj <- addPeakSet(
    ArchRProj = proj,
    peakSet = peakSet_gr,
    force = TRUE
)
proj <- addPeakMatrix(proj)
saveRDS(proj,psate0(args$dir,"/",args$output,"_addPeakMatrix.rds"))