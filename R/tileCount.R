#' Slide windows on a given \link[GenomicRanges]{GRanges} object
#'
#' slidingWindows returns a set of genomic regions by sliding the windows in
#' a given step. Each window is called a "tile". A complementary function for
#' \link[GenomicRanges]{slidingWindows}.
#'
#' @param targetRegions A \link[GenomicRanges]{GRanges} object of
#' genomic regions of interest.
#' @param windowSize numeric(1) or integer(1). Size of windows.
#' @param step numeric(1) or integer(1). Step of windows.
#' @param keepPartialWindow logical(1). Keep last partial window or not.
#' @param ... Not used.
#'
#' @return A \link[GenomicRanges]{GRanges} object.
#'
#' @import GenomicRanges
#' @import S4Vectors
#' @export
#' @examples
#' genes <- GRanges(
#' seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
#' ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600,
#'                    4000, 7500, 5000, 5400),
#'                    width=c(rep(500, 3), 600, 900, 500, 300, 900,
#'                       300, 500, 500),
#'                   names=letters[seq_len(11)]))
#' tg <- slidingWindows(genes, windowSize=50, step=10)
slidingWindows <- function(targetRegions, windowSize,
                        step, keepPartialWindow=FALSE, ...){
    if(any(width(targetRegions)<windowSize) | any(width(targetRegions)<step)){
        warning("Some of chromosomes are smaller than windowSize or step.",
                "They will be removed.")
    }
    targetRegions <- targetRegions[width(targetRegions)>=windowSize]
    targetRegions <- targetRegions[width(targetRegions)>=step]
    tileTargetRanges <- tile(x=ranges(targetRegions), width=step)
    nt <- elementNROWS(tileTargetRanges)
    tileTargetRanges.end <- rep(end(targetRegions), nt)
    tileTargetRanges <- unlist(tileTargetRanges)
    width(tileTargetRanges) <- windowSize
    tileTargetRegions <- GRanges(rep(seqnames(targetRegions), nt),
                                 tileTargetRanges,
                                 rep(strand(targetRegions), nt),
                                 oid=rep(seq_along(targetRegions), nt))
    id <- end(tileTargetRanges) > tileTargetRanges.end
    if(keepPartialWindow){
        end(tileTargetRegions[id]) <- tileTargetRanges.end[id]
    }else{
        tileTargetRegions <-
            tileTargetRegions[!id]
    }
    return(tileTargetRegions)
}

#' Count overlapping genomic ranges
#'
#' Count the reads in a given feature. This function does not work for parallel.
#'
#' @param features A object of \link[GenomicRanges]{GRanges} represents the 
#' feature regions to be counted.
#' @param reads object that represents the data to be counted. See 
#' \link[GenomicAlignments]{summarizeOverlaps}.
#' @param ignore.strand logical(1). ignore strand?
#' @param inter.feature not used. This parameter is required by 
#' \link[GenomicAlignments]{summarizeOverlaps}.
#'
#' @return return a vector of counts the same length as features.
#'
countByOverlaps <- function(features, reads,  ignore.strand, inter.feature) {
    ## NOT work for parallel
    countOverlaps(features, reads, ignore.strand=ignore.strand)
}

#' Perform overlap queries between reads and genome by windows
#'
#' tileCount extends \link[GenomicAlignments]{summarizeOverlaps} by providing
#' fixed window size and step to split whole genome into windows and then do
#' queries. It will return counts in each window.
#'
#' @param reads A \link[GenomicRanges]{GRanges},
#' \link[GenomicRanges]{GRangesList} (should be one read per list element),
#' \link[GenomicAlignments]{GAlignments},
#' \link[GenomicAlignments]{GAlignmentsList},
#' \link[GenomicAlignments]{GAlignmentPairs} or
#' \link[Rsamtools]{BamFileList} object that represents the data to be
#' counted by \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#' @param genome The object from/on which to get/set the sequence information.
#' @param windowSize numeric(1) or integer(1). Size of windows.
#' @param step numeric(1) or integer(1). Step of windows.
#' @param keepPartialWindow logical(1). Keep last partial window or not.
#' @param mode mode can be one of the pre-defined count methods.
#' see \link[GenomicAlignments]{summarizeOverlaps}.
#' default is countByOverlaps, alia of
#' countOverlaps(features, reads, ignore.strand=ignore.strand)
#' @param dataOverSamples logical(1). Data over several samples when use 
#' \link[GenomicRanges]{GRangesList} as input?
#' @param ... Additional arguments passed to
#' \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment} object. 
#' The assays slot holds the counts, rowRanges holds the annotation from 
#' sliding widows of genome.
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import SummarizedExperiment
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics lengths
#' @importFrom methods as
#' @export
#' @examples
#' \dontrun{
#' fls <- list.files(system.file("extdata", package="GenomicAlignments"),
#' recursive=TRUE, pattern="*bam$", full=TRUE)
#' names(fls) <- basename(fls)
#' genes <- GRanges(seqlengths = c(chr2L=7000, chr2R=10000))
#' se <- tileCount(fls, genes, windowSize=1000, step=500)
#' }
#'
#' ##
#' genome <- GRanges("chr1", IRanges(1, 1))
#' seqlengths(genome) <- c(chr1=1000)
#' reads <- GRanges("chr1", IRanges((seq_len(90))*10, width=10))
#' tileCount(reads, genome, windowSize=100, step=50)
#'

tileCount <- function(reads, genome, windowSize=1e5, step=1e4,
                      keepPartialWindow=FALSE,
                      mode=countByOverlaps, 
                      dataOverSamples=FALSE, ...){
    targetRegions <- as(seqinfo(genome), "GRanges")
    tileTargetRegions <- slidingWindows(targetRegions, windowSize,
                                     step, keepPartialWindow)
    if(inherits(reads, "GRangesList")&&dataOverSamples){
        se <- do.call(cbind, lapply(reads, summarizeOverlaps, 
                                    features=tileTargetRegions,
                                    mode=mode, ...))
        if(length(names(reads))==ncol(se)) colnames(se) <- names(reads)
    }else{
        if(inherits(reads, "GRangesList")){
            ol <- lapply(head(reads, 100), findOverlaps, 
                         drop.redundant=TRUE, drop.self=TRUE)
            if(any(lengths(ol))){
                stop("reads are GRangesList with overlaps in elements.",
                     "Please try to set dataOverSamples as TRUE.")
            }
        }
        se <- summarizeOverlaps(features=tileTargetRegions, reads=reads,
                                mode=mode, ...)
    }
    
    names(assays(se)) <- "counts"
    se
}
