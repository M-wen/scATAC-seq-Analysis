#!/usr/bin/env Rscript
# create parser object
parser <- argparse::ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="list of summit bed")
parser$add_argument("-g", "--genome", default = "mm10", help="used genome [default %(default)s]")
parser$add_argument("--chromSize", default = "/projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes", help="chrom size")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("BSgenome"))
library(tictoc)

inF = args$input
chromF = args$chromSize
outDir = args$dir
outF = args$output

source("./bin/function_peakmergeing.R")
################################
# working on peaks

# load summit
summitF <- read.table(inF, sep="\t", header=F)
label.lst <- as.character(summitF$V1)
file.lst <- as.character(summitF$V2)


tic("parse summit peak set")
peak.list = lapply(seq(file.lst), function(i){
  message("working on summit set for... ", label.lst[i])
  p.gr <- read2gr(file.lst[i], label=label.lst[i])
  p.gr <- extendSummit(p.gr, size=500)
  p.gr <- filter4chrom(p.gr, chromF)
  p.gr <- filter4N(p.gr, genome=genome)
  p.gr <- nonOverlappingGR(p.gr, by = "score", decreasing = TRUE)
  p.gr <- norm2spm(p.gr, by="score")
  p.gr
})
toc()


tic("write filtered & fixed peak set")
for(i in 1:length(label.lst)){
  message("write fixed & filtered peak set to: ", outDir, label.lst[i], ".filteredNfixed.peakSet")
  outPeak <- as.data.frame(peak.list[[i]])
  outFname <- paste(outDir, label.lst[i], ".filterNfixed.peakset", sep="")
  fwrite(outPeak, file=outFname, sep="\t", quote = F, col.names = T, row.names = F)
}
toc()


tic("merge to union peak list")
# merge
merged.gr <- do.call(c, peak.list)
merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)
toc()

tic("filter union peak list")
# filter reproducible peaks by choosing a spm cut-off of 2 
merged.filtered.gr <- merged.gr[which(mcols(merged.gr)$spm >2),]


tic("save union peak set")
outUnion <- as.data.frame(merged.filtered.gr)
outfname = paste(outDir, outF, ".filteredNfixed.union.peakSet",sep="")
fwrite(outUnion, file=outfname, sep="\t", quote = F, col.names = T, row.names = F)
toc()
