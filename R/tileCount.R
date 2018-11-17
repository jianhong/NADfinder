#' Count reads overlapping genomic ranges
#'
#' Count reads overlapping a set of genimc features represented as
#' genomic ranges. This function does not work for parallel.
#' @param features A object of \link[GenomicRanges:GRanges-class]{GRanges} representing the
#' feature regions to be counted.
#' @param reads An object that represents the data to be counted. See
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory. For paired end reads, 
#' @param ignore.strand logical(1). ignore strand?
#' @param inter.feature not used. This parameter is required by
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
#' @export
#' @importFrom stats setNames
#' @return return a summarized experiment object with chromosome-level depth
#' information for each input sample as metadata.
#'
IntersectionNotStrict <-function(features,
                                 reads,
                                 ignore.strand = TRUE,
                                 inter.feature = FALSE) 
{
    ## NOT work for parallel
    ov <- findOverlaps(reads,
                       features,
                       type = "within",
                       ignore.strand = ignore.strand)
    countSubjectHits(ov)
}

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


#' Perform overlap queries between reads and genome by windows
#'
#' tileCount extends \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps} by finding coverage for
#' each fixed window in the whole genome
#'
#' @param reads A \link[GenomicRanges:GRanges-class]{GRanges},
#' \link[GenomicRanges:GRangesList-class]{GRangesList} (should be one read per list element),
#' \link[GenomicAlignments:GAlignments-class]{GAlignments},
#' \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList},
#' \link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs} or
#' \link[Rsamtools:BamFile-class]{BamFileList} object that represents the data to be
#' counted by \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory.
#' @param genome A BSgenome object from/on which to get/set the sequence and metadata information.
#' @param windowSize numeric(1) or integer(1). Size of the windows.
#' @param step numeric(1) or integer(1). Step of generating silding windows.
#' @param mode One of the pre-defined count methods.
#' @param excludeChrs A vector of string: chromosomes/scaffolds of no interest for NAD analysis.
#' see \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
#' default is countByOverlaps, alia of countOverlaps(features, reads, ignore.strand=ignore.strand)
#' @param dataOverSamples logical(1). Data over several samples when use
#' \link[GenomicRanges:GRangesList-class]{GRangesList} as input.
#' @param ... Additional arguments passed to
#' \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
#'
#' @return A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} object.
#' The assays slot holds the counts, rowRanges holds the annotation from the
#' sliding widows of genome.
#' metadata contains lib.size.chrom for holding chromosome-level sequence depth
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import SummarizedExperiment
#' @importFrom IRanges IRanges elementNROWS
#' @importFrom BiocGenerics table
#' @importFrom methods as
#' @importFrom ATACseqQC readBamFile
#' @export
#' @author Jianhong Ou, Haibo Liu, Herve Pages and Julie Zhu
#' @examples
#' if (interactive())
#' {
#'     fls <- list.files(system.file("extdata", package="NADfinder"),
#'     recursive=FALSE, pattern="*bam$", full=TRUE)
#'     names(fls) <- basename(fls)
#'     if (!require(BSgenome.Mmusculus.UCSC.mm10))
#'     {
#'         if (!requireNamespace("BiocManager", quietly=TRUE))
#'         install.packages("BiocManager")
#'         BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#'         library(BSgenome.Mmusculus.UCSC.mm10)
#'     }
#'     se <- tileCount(reads = fls, 
#'                     genome = Mmusculus,
#'                     excludeChrs = c("chrM", paste0("chr", c(1:17,19)), 
#'                                     "chrX", "chrY"), 
#'                     windowSize=50000, step=10000)
#' }
#'
#' 



tileCount<- function(reads,
                      genome,
                      excludeChrs = c("chrM", "M", "Mt", "MT"),
                      windowSize = 50000,
                      step = 10000,
                      mode = IntersectionNotStrict,
                      dataOverSamples = FALSE,
                      ...) 
{
    stopifnot(all(grepl(".bam$", reads)), all(file.exists(paste0(reads, ".bai"))))
    stopifnot(class(genome) == "BSgenome")
    stopifnot(windowSize %% 1 == 0, step %% 1 ==0, windowSize > 0, step > 0, step < windowSize)

    targetRegions <- as(seqinfo(genome), "GRanges")

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
    lib.size.chrom <- lib.size.chrom[!rownames(lib.size.chrom) %in% excludeChrs, drop = FALSE, ]
    
    ## filtering GRanges to keep only those chromosomal scaffolds that are in the BAM file
    targetRegions <- targetRegions[seqnames(targetRegions) %in% rownames(lib.size.chrom)]
    tileTargetRegions <-slidingWindows(x = targetRegions,
                                       width = windowSize,
                                       step = step)
    tileTargetRegionsLen <- elementNROWS(tileTargetRegions)
    tileTargetRegions <- unlist(tileTargetRegions)
    mcols(tileTargetRegions)$oid <- rep(seq_along(targetRegions),tileTargetRegionsLen)
    
     rse_list <- lapply(aln_list,
         function(aln) summarizeOverlaps(features=tileTargetRegions,
                                         reads=aln,
                                         mode=mode,
                                         inter.feature=FALSE))
     rse <- do.call(cbind, unname(rse_list))
     colnames(rse) <- names(rse_list)

    names(assays(rse)) <- "counts"
    metadata(rse)$lib.size.chrom <- lib.size.chrom
    rse
}
