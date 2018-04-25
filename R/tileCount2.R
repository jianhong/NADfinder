#' Perform overlap queries between reads and genome by sliding windows
#' Count reads over sliding windows.
#' @param reads An object that represents the names and path of the bam files  to be counted.
#' If reads are more than 1 bam files,
#' it should be a vector of character with full path. This funciton now if for paired end reads
#' @export
#' @return return a summarized experiment object with chromosome-level depth
#' information for each input sample as metadata.
#'
#'
#' @param windowSize numeric(1) or integer(1). Size of the windows.
#' @param step numeric(1) or integer(1). Step of generating silding windows.
#' @param pe a character string indicating whether paired-end data is present; set to "none", "both", "first" or "second"
#' @param restrict restrict to a set of chromosomes, default to mouse chromosomes.
#' @param filter default to 0 without filtering
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
#' The assays slot holds the counts, rowRanges holds the annotation from the
#' sliding widows of genome.
#' metadata contains lib.size.chrom for holding chromosome-level sequence depth
#' @importFrom csaw windowCounts
#' @import SummarizedExperiment
#' @export
#' @author Jun Yu
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
 
 
 
tileCount2 <- function(reads,
                     windowSize = 50000,
                     restrict = paste0("chr", c(1:19, "X", "Y")),
                     step = 1000,
                     filter = 0,
                     pe = "both") {
    rse <- windowCounts(bam.files = reads, 
               param=readParam(pe=pe, restrict = restrict), 
               filter = filter, width=windowSize, spacing = step)
 
    seqlevels(rse, pruning.mode = "coarse") <- restrict
    rse <- sort(rse)
    colnames(rse) <- basename(reads)
    counts.l <- split(assay(rse), f = seqnames(rowRanges(rse)))
   
    metadata(rse)$lib.size.chrom <- round(t(sapply(counts.l, colSums))*windowSize/step, digits = 0)
    colData(rse)$records <- colData(rse)$totals
    colData(rse)$object <- colnames(rse)
   
    chrs <- seqnames(rowRanges(rse))
    runValue(chrs) <- seq_along(runValue(chrs))
    rowRanges(rse)$oid <- chrs
    return(rse)
}
