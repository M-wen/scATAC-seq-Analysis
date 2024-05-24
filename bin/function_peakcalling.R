## Functions for peak calling
getFrags <- function(ArrowFile, cellNames){
    covList <- lapply(availableChr, function(kk){
        ArchR:::.getFragsFromArrow(ArrowFile,
             chr = kk, out = "GRanges",
             cellNames = cellNames)
    })
    covList <- do.call(c, covList)
    covList <- c(GRanges(seqnames = seqnames(covList),
        ranges = IRanges(start = start(covList), end = start(covList))),
        GRanges(seqnames = seqnames(covList),
        ranges = IRanges(start = end(covList), end = end(covList))))
    covList <- covList + 25
    covList <- sort(covList)
}
concateBeds <- function(x){
    all_beds <- lapply(x, function(y) rtracklayer::import(y))
    all_beds <- as(all_beds, "GRangesList")
    all_beds <- unlist(all_beds)
    return(all_beds)
}
macsPeaks <- function(input, output){
system(paste0("$software_path/ minconda3/envs/scATACPipeline/bin/macs2 callpeak -t ", input,
       " -f BED -n ", output,
       " -g 2871842584 -p 0.01 ",
       "--nomodel --keep-dup all --call-summits --nolambda"))
}
# get fragments from ArrowFiles for every sample and every celltype
# split into two pseudo-pseudobulks of equal size

make_beds_pseudo <- function(ArrowFile, ArrowFileName,
                             cellGroups_new=cellGroups_new){
    for(i in 1:length(cellGroups_new)){

        cellGroupi <- cellGroups_new[[i]]
        cellGroupi <- cellGroupi[grepl(paste0(ArrowFileName, "#"),
                        cellGroupi)]
        if(length(cellGroupi)<40) next

        ind_groupi <- sample(1:length(cellGroupi), floor(length(cellGroupi)/2))
        ind_groupi <- list(pseudo_1=ind_groupi,
            pseudo_2=setdiff(1:length(cellGroupi), ind_groupi))

        for(k in 1:2){

            cellNames <- cellGroupi[ind_groupi[[k]]]

            covList <- getFrags(ArrowFile, cellNames)

            # remove blacklist and set strand to positive
            covList <- subsetByOverlaps(covList, blacklist, invert=TRUE)
            strand(covList) <- "+"

            rtracklayer::export(covList, con=paste0("snATACseq/bedfiles_pseudo/",
                names(cellGroups_new)[i],"_", ArrowFileName, "_pseudo", k,
                "_blacklistrm.bed"), format="bed")
        }
    }
}
call_pseudo <- function(cellGroup, ArrowFile, rep=2){

   for(k in 1:rep){
       macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cellGroup, "_", ArrowFile,
       "_", "pseudo", k, "_blacklistrm.bed"),
       paste0("snATACseq/bedfiles_pseudo/", cellGroup, "_", ArrowFile,
       "_pseudo", k))
   }
}

call_pseudo_cluster <- function(cell_type, bed_files){

    pseudo_files <- bed_files[grepl(cell_type, bed_files)]
    pseudo_files <- pseudo_files[grepl(cell_type, pseudo_files)]
    pseudo_files <- pseudo_files[grepl("_blacklistrm.bed", pseudo_files)]

    pseudo_files1 <- paste0("snATACseq/bedfiles_pseudo/",
        pseudo_files[grepl("_pseudo1", pseudo_files)])
    pseudo_files2 <- paste0("snATACseq/bedfiles_pseudo/",
        pseudo_files[grepl("_pseudo2", pseudo_files)])

    all_pseudo_files_1 <- concateBeds(pseudo_files1)
    export(all_pseudo_files_1, con=paste0("snATACseq/bedfiles_pseudo/",
        cell_type, "_pseudo1.bed"), format="bed")
    all_pseudo_files_2 <- concateBeds(pseudo_files2)
    export(all_pseudo_files_2, con=paste0("snATACseq/bedfiles_pseudo/",
        cell_type, "_pseudo2.bed"), format="bed")

    macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo1.bed"),
        paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo1"))
    macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo2.bed"),
        paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo2"))

}
make_beds<- function(ArrowFile, ArrowFileName,
                             cellGroups_new=cellGroups_new){
    for(i in 1:length(cellGroups_new)){
        cellGroupi <- cellGroups_new[[i]]
        cellGroupi <- cellGroupi[grepl(paste0(ArrowFileName, "#"),
                                   cellGroupi)]
        if(length(cellGroupi)<40) next

        covList <- getFrags(ArrowFile, cellGroupi)
        covList <- subsetByOverlaps(covList, blacklist, invert=TRUE)

        strand(covList) <- "+"

        rtracklayer::export(covList, con=paste0("snATACseq/bedfiles/",
            names(cellGroups_new)[i], "_", ArrowFileName,
            "_blacklistrm.bed"), format="bed")
    }
}
call_cluster <- function(cell_type, bed_files){

    pseudo_files <- bed_files[grepl(cell_type, bed_files)]
    pseudo_files <- pseudo_files[grepl(cell_type, pseudo_files)]
    pseudo_files <- pseudo_files[grepl("_blacklistrm.bed", pseudo_files)]

    all_pseudo_files <- concateBeds(paste0("snATACseq/bedfiles/", pseudo_files))
    export(all_pseudo_files, con=paste0("snATACseq/bedfiles/",
        cell_type, "_blacklistrm.bed"), format="bed")

    macsPeaks(paste0("snATACseq/bedfiles/", cell_type, "_blacklistrm.bed"),
        paste0("snATACseq/bedfiles/", cell_type))

}

filter_peaks <- function(celltype){
Peaks <- list.files("03_peakcalling/bedfiles",".narrowPeak")
reps <- paste0("03_peakcalling/bedfiles/",Peaks[grep(celltype,Peaks)])
pooled_peak <- paste0("03_peakcalling/bedfiles/",celltype,"_peaks.narrowPeak")

## Find pooled peaks that overlap Reps where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs > 0.5
pooledInReps <- mclapply(reps,function(x){
    pooled_gr <- rtracklayer::import(pooled_peak)
    rep_gr <- rtracklayer::import(x)
    pooledInRep <- overlap_peak_gr(pooled_gr,rep_gr)
},mc.cores = 3)

##  peaks in both individual replicates were kept
pairs <- combinat::combn(seq_along(pooledInReps),2,simplify = FALSE)
pooledInBothReps <- mclapply(pairs,function(x) overlap_peak_gr(pooledInReps[[x[1]]],pooledInReps[[x[2]]]),mc.cores = 3)

pooledInBothReps_df <- mclapply(seq_along(pooledInBothReps), function(x){
    df <- as.data.frame(pooledInBothReps[[x]]) 
    # df <- df[order(df$seqnames,df$start),]
    df <- unique(df)
},mc.cores = 3) %>% do.call(rbind,.) 

pooledInBothReps_df <- unique(pooledInBothReps_df)

## Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs > 0.5
# psreps <- c(pseudo1_peak,pseudo2_peak)
pseudoPeaks <- list.files("03_peakcalling/bedfiles_pseudo",".narrowPeak")
psreps <- paste0("03_peakcalling/bedfiles_pseudo/",pseudoPeaks[grep(celltype,pseudoPeaks)])
pooledInPsReps <- mclapply(psreps,function(x){
    pooled_gr <- rtracklayer::import(pooled_peak)
    psrep_gr <- rtracklayer::import(x)
    pooledInpsRep <- overlap_peak_gr(pooled_gr,psrep_gr)
},mc.cores = 2)
##  peaks in both pseudo-replicates were kept
pooledInBothPsReps <- overlap_peak_gr(pooledInPsReps[[1]],pooledInPsReps[[2]])
pooledInBothPsReps_df <- unique(as.data.frame(pooledInBothPsReps))

## Combine peak lists
final_df <- unique(rbind(pooledInBothReps_df,pooledInBothPsReps_df))
final_df$score[final_df$score > 1000]  <- 1000
chr_pattern <- '^chr[0-9XY]+$'
final_df <- final_df[grepl(chr_pattern, final_df$seqnames),]
final_df <- final_df[,-4]
fwrite(final_df,paste0("04_parsePeak/",celltype,".naivePeakList.narrowPeak"),col.names = F,sep="\t")                             
# rtracklayer::export(final_df, con=paste0("04_parsePeak/",celltype,".naivePeakList.narrowPeak"), format="bed")

## Get summit
pooled_summit <- paste0("03_peakcalling/bedfiles/",celltype,"_summits.bed")
pooled_summit_df <- as.data.frame(fread(pooled_summit,header = F))
pooled_summit_df <- pooled_summit_df[which(pooled_summit_df$V4 %in% final_df$name),]
fwrite(pooled_summit_df,paste0("04_parsePeak/",celltype,".naiveSummitList.bed"),col.names = F,sep="\t")
# rtracklayer::export(pooled_summit_df, con=paste0("04_parsePeak/",celltype,".naiveSummitList.bed"), format="bed")

path=paste0(getwd(),"/04_parsePeak/",celltype,".naiveSummitList.bed")
cmd <- paste0("echo -e \"",celltypes[1],"\t",path,"\"" ,">> pbmc.naiveSummitList.list" )
system(cmd)

}
                       

