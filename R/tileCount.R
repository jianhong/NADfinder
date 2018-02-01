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
#' @return return a summarized experiment object with chromosome-level depth 
#' information for each input sample as metadata.
#'
IntersectionNotStrict <- function(features, reads,  ignore.strand = TRUE, inter.feature = FALSE) {
    ## NOT work for parallel
    ov <- findOverlaps(reads, features, type = "within", ignore.strand = ignore.strand)
    countSubjectHits(ov)
}

#' Perform overlap queries between reads and genome by windows
#'
#' tileCount extends \link[GenomicAlignments]{summarizeOverlaps} by finding coverage for
#' each fixed window in the whole genome
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
#' The assays slot holds the counts, rowRanges holds the annotation from the
#' sliding widows of genome.
#' metadata contains lib.size.chrom for holding chromosome-level sequence depth
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import SummarizedExperiment
#' @importFrom IRanges IRanges elementNROWS
#' @importFrom BiocGenerics lengths table
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
#' genome <- GRanges(c("chr1","chr2"), IRanges(c(1,1), c(1000,1000)))
#' seqlengths(genome) <- c(chr1=1000, chr2=1000)
#' reads <- GRanges("chr1", IRanges((seq_len(90))*10, width=10))
#' tileCount(reads, genome, windowSize=100, step=50)
#' reads.2 <- GRangesList(GRanges("chr2", IRanges((seq_len(90))*10, width=10)), reads)
#' tileCount(reads.2, genome, windowSize=100, step=50)
#' @author Jianhong Ou and Julie Zhu

tileCount <- function(reads, genome, windowSize=1e5L, step=1e4L,
                      mode=IntersectionNotStrict, 
                      dataOverSamples=FALSE, ...){
    targetRegions <- as(seqinfo(genome), "GRanges")
    tileTargetRegions <- 
        slidingWindows(x=targetRegions, 
                       width=windowSize,
                       step=step)
    tileTargetRegionsLen <- elementNROWS(tileTargetRegions)
    tileTargetRegions <- unlist(tileTargetRegions)
    mcols(tileTargetRegions)$oid <- rep(seq_along(targetRegions), 
                                        tileTargetRegionsLen)
    if(inherits(reads, "GRangesList")&&dataOverSamples){
        se <- do.call(cbind, lapply(reads, summarizeOverlaps, 
                                    features=tileTargetRegions,
                                    inter.feature = FALSE,
                                    mode=mode, ...))
        if(length(names(reads))==ncol(se)) colnames(se) <- names(reads)
        lib.size.chrom <- t(BiocGenerics::table(seqnames(reads))) 
        lib.size.chrom <- cbind(rownames(lib.size.chrom), lib.size.chrom)
        colnames(lib.size.chrom) <- c("chrom", paste("lib.size", 1:length(reads), sep="."))
    }else{
        if(inherits(reads, "GRangesList")){
            ol <- lapply(head(reads, 100), findOverlaps, 
                         drop.redundant=TRUE, drop.self=TRUE)
            if(any(lengths(ol))){
                stop("reads are GRangesList with overlaps in elements.",
                     "Please try to set dataOverSamples as TRUE.")
            }
            lib.size.chrom <- as.data.frame(rowSums(t(BiocGenerics::table(seqnames(reads)))))
            lib.size.chrom <- cbind(rownames(lib.size.chrom), lib.size.chrom)
            colnames(lib.size.chrom) <- c("chrom", "lib.size")
        }
        else if (inherits(reads, "GRanges")) {
            lib.size.chrom <- as.data.frame(BiocGenerics::table(seqnames(reads))) 
            colnames(lib.size.chrom) <- c("chrom", "lib.size")
        }
        else {
            aln <- scanBam(reads)
            lib.size.chrom <- as.data.frame(BiocGenerics::table(aln[[1]]$rname))
            if (length(aln) > 1)
            {
               lib.size.chrom <- cbind(lib.size.chrom, do.call(cbind, lapply(2:length(aln),function(i){
                 as.data.frame(BiocGenerics::table(aln[[i]]$rname))[,2] } )))
            }
        }
        se <- summarizeOverlaps(features=tileTargetRegions, reads=reads,
                                mode=mode, inter.feature = FALSE, ...)
    }
     
    names(assays(se)) <- "counts"
    metadata(se)$lib.size.chrom <- lib.size.chrom
    se
}
