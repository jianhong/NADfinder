#' Perform overlap queries between reads and genome by sliding windows
#' Count reads over sliding windows.
#' @param reads An object that represents the names and path of the bam files  to be counted.
#' If reads are more than 1 bam files,
#' it should be a vector of character with full path. This function now works for paired end reads
#'
#' @param fragment.length integer(1). An integer scalar or a list of two integer scalars/vectors,
#'          containing the average length(s) of the sequenced fragments in each libary.
#' @param windowSize numeric(1) or integer(1). Size of the windows.
#' @param step numeric(1) or integer(1). Step of generating silding windows.
#' @param pe a character string indicating whether paired-end data is present; set to "none", "both", "first" or "second"
#' @param restrict restrict to a set of chromosomes, default to mouse chromosomes.
#' @param filter default to 0 without filtering. An integer scalar for the minimum count sum across libraries for each window
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment} object with chromosome-level depth
#' The assays slot holds the counts, rowRanges holds the annotation from the
#' sliding widows of genome.
#' metadata contains lib.size.chrom for holding chromosome-level sequence depth
#' @importFrom csaw windowCounts
#' @import SummarizedExperiment
#' @import rbamtools
#' @import Rsamtools
#' @import GenomicAlignments
#' @export
#' @author Jun Yu,Hervé Pagès and Julie Zhu
#' @examples
#' if (interactive())
#' {
#'     fls <- list.files(system.file("extdata", package="NADfinder"),
#'     recursive=FALSE, pattern="*bam$", full=TRUE)
#'     names(fls) <- basename(fls)
#'    
#'     se <- tileCount2(reads = fls,
#'                     windowSize=50000, step=10000)
#' }
#'
#'
 
 

computeLibSizeChrom <- function(aln_list)
     {
         stopifnot(is.list(aln_list))
         lib_size_list <- lapply(aln_list,
             function(aln) {
                 qname <- names(aln)
                 if (is.null(qname))
                     stop(wmsg("Some of the GAlignments or
GAlignmentsList ",
                               "objects in 'aln_list' don't have names. ",
                               "Did you use 'use.names=TRUE' when loading ",
                               "them with readGAlignments() or ",
                               "readGAlignmentsList()?"))
                 if (is(aln, "GAlignmentsList")) {
                     rname <- seqnames(unlist(aln, use.names=FALSE))
                     qname <- rep.int(qname, lengths(aln, use.names=FALSE))
                 } else {
                     rname <- seqnames(aln)
                 }
                 lengths(unique(split(qname, rname)))
             })
         rnames <- unique(names(unlist(unname(lib_size_list))))
         lib_size_list <- lapply(lib_size_list,
             function(lib_size) {
                 lib_size2 <- setNames(integer(length(rnames)), rnames)
                 lib_size2[names(lib_size)] <- lib_size
                 lib_size2
             })
         ans <- do.call(cbind, unname(lib_size_list))
         colnames(ans) <- names(lib_size_list)
         ans
    }

 
tileCount2 <- function(reads,
                     fragment.length = 100,
                     windowSize = 50000,
                     restrict = paste0("chr", c(1:19, "X", "Y")),
                     step = 1000,
                     filter = 0,
                     pe = "both") {
    rse <- windowCounts(bam.files = reads, ext = fragment.length, 
               param=readParam(pe=pe, restrict = restrict), 
               filter = filter, width=windowSize, spacing = step)
 
    seqlevels(rse, pruning.mode = "coarse") <- restrict
    rse <- sort(rse)
    colnames(rse) <- basename(reads)
    counts.l <- split(assay(rse), f = seqnames(rowRanges(rse)))

#### bamCountAll does not count properly with paired-end data  
    
    if (!testPairedEndBam(reads[1]))
    {
        count <- lapply(reads, function(thisBam) {
    	    reader <- bamReader(thisBam, idx = TRUE) 
            temp <- bamCountAll(reader, verbose = FALSE)
            temp1 <- cbind(rownames(temp), temp$nAligns)
            colnames(temp1) <- c("chr", basename(thisBam))
            temp1
       })
   
       merge.all <- function(x, y) {
           merge(x, y, all=TRUE, by="chr")
       }

       lib.size.chrom <- Reduce(merge.all, count)
       rownames(lib.size.chrom) <- lib.size.chrom[,1]
       lib.size.chrom <- lib.size.chrom[, -1]
       lib.size.chrom <- lib.size.chrom[rownames(lib.size.chrom) %in% restrict, drop = FALSE, ]
    }
    else
    { 
        if (is.null(names(reads)))
          names(reads) <- basename(reads)
        aln_list <- lapply(reads,
           function(file) {
             isPE <- testPairedEndBam(file)
             if (isPE)
                 readGAlignmentsList(file, use.names=TRUE)
             else
                 readGAlignments(file, use.names=TRUE)

         })

       ## aln is a list of list of lists
       lib.size.chrom <- computeLibSizeChrom(aln_list)
    
       ## remove excludeChr from lib.size.chrom
       lib.size.chrom <- lib.size.chrom[rownames(lib.size.chrom) %in% restrict, drop = FALSE, ]
    } 
    metadata(rse)$lib.size.chrom <- lib.size.chrom 
    colData(rse)$records <- colData(rse)$total
    colData(rse)$object <- colnames(rse)
   
    chrs <- seqnames(rowRanges(rse))
    runValue(chrs) <- seq_along(runValue(chrs))
    rowRanges(rse)$oid <- chrs
    return(rse)
}
