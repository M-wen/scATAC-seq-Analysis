
# read bed to gr
read2gr <- function(bedF, label){
  df <- fread(bedF, sep="\t", header=F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  );
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}


# extend summit to 500 bp
extendSummit <- function(gr, size=500){
  gr <- resize(gr, width=size, fix="center")
  return(gr)
}

# filter non-chromosome
filter4chrom <- function(gr, chromF){
  chrom = read.table(chromF);
  chrom.gr = GRanges(
    chrom[,1],
    IRanges(0, chrom[,2])
  );
  idx = queryHits(
    findOverlaps(gr, chrom.gr, type="within")
  );
  if(length(idx) > 0){
    gr = gr[idx]
  }
  return(gr)
}


# filter N containing regions
filter4N <- function(gr, genome=genome){
  genome <- getBSgenome(genome)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(gr)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  gr[which(mcols(gr)$N < 0.001)] #Remove N Containing Peaks
  return(gr)
}

# get non-overlapped regions
#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
  ){
  
  stopifnot(by %in% colnames(mcols(gr)))

  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)

    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
    
    }else{
      grAll <- c(grAll, grSelect)
    } 

  }
  message(sprintf("Converged after %s iterations!", i))

  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))

  return(grAll)

}


# normlize to score per million

norm2spm <- function(gr, by = "score"){
  mlogp = mcols(gr)[,by]
  normmlogp = 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}

