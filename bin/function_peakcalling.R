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


