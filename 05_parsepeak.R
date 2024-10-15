library(ArchR)
library(parallel)

overlap_peak <- function(pooled_gr, rep_gr) {  
    overlap_pooled_rep <- findOverlaps(pooled_gr, rep_gr)  
    overlap_pooled <- pooled_gr[queryHits(overlap_pooled_rep)]
    overlap_rep1gr <- rep_gr[subjectHits(overlap_pooled_rep)]
    overlapgr <- pintersect(overlap_pooled,overlap_rep1gr)
    filtered_gr <- overlap_pooled[width(overlapgr)/width(overlap_pooled) >= 0.5 | width(overlapgr)/width(overlap_rep1gr) >= 0.5]  
    # filtered_gr <- unique(filtered_gr)
    peak_name <- unique(filtered_gr$name)
    return(peak_name)
}

filter_peaks <- function(celltype){
Peaks <- list.files("03_peakcalling/bedfiles",".narrowPeak")
reps <- paste0("03_peakcalling/bedfiles/",Peaks[grep(celltype,Peaks)])
pooled_peak <- paste0("03_peakcalling/bedfiles/",celltype,"_peaks.narrowPeak")
reps <- reps[which(!reps %in% pooled_peak)]

pooled_gr <- rtracklayer::import(pooled_peak)
## Find pooled peaks that overlap Reps where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
pooledInReps <- lapply(reps,function(x){
    # pooled_gr <- rtracklayer::import(pooled_peak)
    rep_gr <- rtracklayer::import(x)
    pooledInRep <- overlap_peak(pooled_gr,rep_gr)
    return(pooledInRep)
})

##  peaks in both individual replicates were kept
pooledInRepsStat <- table(unlist(pooledInReps))
pooledInBothRepsPeak <- names(pooledInRepsStat[pooledInRepsStat > 1])

## Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
# psreps <- c(pseudo1_peak,pseudo2_peak)
pseudoPeaks <- list.files("03_peakcalling/bedfiles_pseudo",".narrowPeak")
psreps <- paste0("03_peakcalling/bedfiles_pseudo/",pseudoPeaks[grep(celltype,pseudoPeaks)])
pooledInPsReps <- lapply(psreps,function(x){
    # pooled_gr <- rtracklayer::import(pooled_peak)
    psrep_gr <- rtracklayer::import(x)
    pooledInpsRep <- overlap_peak(pooled_gr,psrep_gr)
    return(pooledInpsRep)
})
##  peaks in both pseudo-replicates were kept
# pooledInBothPsReps <- overlap_peak_gr(pooledInPsReps[[1]],pooledInPsReps[[2]])
pooledInPsRepsStat <- table(unlist(pooledInPsReps))
pooledInBothPsRepsPeak <- names(pooledInPsRepsStat[pooledInPsRepsStat > 1])
rm(pooledInPsReps)
rm(pooled_gr)
# pooledInBothPsReps_df <- unique(as.data.frame(pooledInBothPsReps))

## Combine peak lists
final_peak <- c(pooledInBothRepsPeak,pooledInBothPsRepsPeak)
pooled_peak_df <- fread(pooled_peak,header = F)
# final_df <- unique(rbind(pooledInBothReps_df,pooledInBothPsReps_df))
finle_df <- pooled_peak_df[which(pooled_peak_df$V4 %in% final_peak),]
final_df$V5[final_df$V5 > 1000]  <- 1000
chr_pattern <- '^chr[0-9XY]+$'
final_df <- final_df[grepl(chr_pattern, final_df$V1),]
# final_df <- final_df[,-4]
fwrite(final_df,paste0("04_parsePeak/",celltype,".naivePeakList.narrowPeak"),col.names = F,sep="\t")                             
# rtracklayer::export(final_df, con=paste0("04_parsePeak/",celltype,".naivePeakList.narrowPeak"), format="bed")

## Get summit
pooled_summit <- paste0("03_peakcalling/bedfiles/",celltype,"_summits.bed")
pooled_summit_df <- as.data.frame(fread(pooled_summit,header = F))
pooled_summit_df <- pooled_summit_df[which(pooled_summit_df$V4 %in% final_peak),]
fwrite(pooled_summit_df,paste0("04_parsePeak/",celltype,".naiveSummitList.bed"),col.names = F,sep="\t")
# rtracklayer::export(pooled_summit_df, con=paste0("04_parsePeak/",celltype,".naiveSummitList.bed"), format="bed")

path=paste0(getwd(),"/04_parsePeak/",celltype,".naiveSummitList.bed")
cmd <- paste0("echo -e \"",celltype,"\t",path,"\"" ,">> pbmc.naiveSummitList.list" )
system(cmd)

}
setwd("/hwfssz1/ST_SUPERCELLS/P21Z10200N0090/mawen3/project/03.snATAC-seq_analysis_pipeline/results")
system("mkdir 04_parsePeak")

sc <- readRDS("02_clustering/pbmc_scATAC_annotation.rds")
sc$predictedsubGroup <- gsub(" ","_",sc$predictedsubGroup)
sc$predictedsubGroup <- gsub("\\+","",sc$predictedsubGroup)
ArrowFiles <- getArrowFiles(sc)
Groups <- getCellColData(ArchRProj = sc, select = "predictedsubGroup", drop = TRUE)
Cells <- sc$cellNames
cellGroups <- split(Cells, Groups)
input <- lapply(1:length(cellGroups), function(x) lapply(names(ArrowFiles),
     function(y) {
        if(sum(grepl(paste0(y, "#"), cellGroups[[x]]))>=40)
        c(names(cellGroups)[x], y)
      }
))
input <- unlist(input, recursive = FALSE)           
input <- input[!sapply(input, is.null)]
celltypes <- unique(sapply(seq_along(input),function(x){y=input[[x]][1] 
                                   return(y)}))
celltypes  

#cores <- length(celltypes)
mclapply(celltypes,function(x) filter_peaks(x),mc.cores = cores)