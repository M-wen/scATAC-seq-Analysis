#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="annotated rds file")
parser$add_argument("--blacklist", help="blacklist file")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")

options(scipen=999)

library(ArchR)
library(rtracklayer)
library(parallel)
library(readxl)
library(data.table)
library(GenomicRanges)
addArchRGenome("hg38")
addArchRThreads(threads = 1)
set.seed(10)

source("./bin/function_peakcalling.R")

# 1. read in data and blacklist
setwd(args$dir)
sc <- readRDS(args$input)
blacklist <- import(args$blacklist)
sc$predictedsubGroup <- gsub(" ","_",sc$predictedsubGroup)
sc$predictedsubGroup <- gsub("\\+","",sc$predictedsubGroup)

ArrowFiles <- getArrowFiles(sc)
Groups <- getCellColData(ArchRProj = sc, select = "predictedsubGroup", drop = TRUE)
Cells <- sc$cellNames
cellGroups <- split(Cells, Groups)
availableChr <- ArchR:::.availableSeqnames(head(getArrowFiles(sc)))
chromLengths <- getChromLengths(sc)
chromSizes <- getChromSizes(sc)
cell_types <- names(cellGroups)
input <- lapply(1:length(cellGroups), function(x) lapply(names(ArrowFiles),
     function(y) {
        if(sum(grepl(paste0(y, "#"), cellGroups[[x]]))>=40)
        c(names(cellGroups)[x], y)
      }
))
input <- unlist(input, recursive = FALSE)           
input <- input[!sapply(input, is.null)]

system("mkdir -p 03_peakcalling/bedfiles_pseudo")
system("mkdir -p 03_peakcalling/bedfiles")
system("mkdir -p 03_peakcalling/peaks")
              
# get fragments from ArrowFiles for every sample and every celltype and split into two pseudo-pseudobulks of equal size
mcmapply(function(X,Y) make_beds_pseudo(X,Y, cellGroups=cellGroups),
    X=ArrowFiles, Y=names(ArrowFiles),
    mc.cores=10)

bed_files <- list.files("03_peakcalling/bedfiles_pseudo/","_blacklistrm.bed")
celltypes <- unique(sapply(seq_along(input),function(x){y=input[[x]][1] 
                                   return(y)}))

# make bedfiles and call peaks for pseudo-pseudobulks for cell types
mclapply(celltypes, function(x) call_pseudo_cluster(x, bed_files), mc.cores=1) 

mclapply(function(X,Y) make_beds(X, Y, cellGroups),
  X=ArrowFiles, Y=names(ArrowFiles),
  mc.cores=10)
# make bedfiles for pseudobulk for celltype and sample
mclapply(input, function(x) call_peaks(x[1], x[2]), mc.cores=10)
bed_files <- list.files("03_peakcalling/bedfiles")
# make bedfiles and call peaks for pseudobulks for each celltype
mclapply(celltypes, function(x) call_cluster(x, bed_files), mc.cores=10)

