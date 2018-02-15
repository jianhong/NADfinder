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
#'
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
#' \link[Rsamtools:BamFileList-class]{BamFileList} object that represents the data to be
#' counted by \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory.
#' @param genome A BSgenome object from/on which to get/set the sequence and metadata information.
#' @param windowSize numeric(1) or integer(1). Size of the windows.
#' @param step numeric(1) or integer(1). Step of generating silding windows.
#' @param mode One of the pre-defined count methods.
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
#' @importFrom BiocGenerics lengths table
#' @importFrom methods as
#' @importFrom ATACseqQC readBamFile
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
#' @author Jianhong Ou, Haibo Liu and Julie Zhu

tileCount <- function(reads,
                      genome,
                      windowSize = 1e5L,
                      step = 1e4L,
                      mode = IntersectionNotStrict,
                      dataOverSamples = FALSE,
                      ...) 
{
  ## get sliding windows with size = windowSize and step = step
  targetRegions <- as(seqinfo(genome), "GRanges")
  tileTargetRegions <-slidingWindows(x = targetRegions,
                                     width = windowSize,
                                     step = step)
  tileTargetRegionsLen <- elementNROWS(tileTargetRegions)
  tileTargetRegions <- unlist(tileTargetRegions)
  mcols(tileTargetRegions)$oid <- rep(seq_along(targetRegions),tileTargetRegionsLen)
  
  ## read in BAM files in paired  or single end mode if the input is paired end reads or SE, 
  ## respectively using the readBamFile function in the ATACseqQC package. The output will be
  isPE <- testPairedEndBam(reads[1])
  
  ## get fragment counts per chromosome
  countByChr <- function(algns)
  {
    fragments <- lapply(algns, function (pairedAlignmentItems) {unique(pairedAlignmentItems@seqnames)}) 
    countFreqByChr <- as.data.frame(table(unlist(fragments)))
    countFreqByChr
  }
  
  ## add names and metadata to a RangedSummarizedExperiment object
  addNamesAndMetadata <- function(rse, lib.size.chrom)
  {
    names(assays(rse)) <- "counts"
    metadata(rse)$lib.size.chrom <- lib.size.chrom
    rse
  }
  
  if (isPE)
  {
    ## get count of fragments for each sliding window: a list 
    pe <- lapply(reads, function(x) 
    {
      peCount <- summarizeOverlaps(features = tileTargetRegions,
                                   reads = x,
                                   mode = mode,
                                   singleEnd = FALSE,
                                   fragments = TRUE,
                                   inter.feature = FALSE, ...)
      GAlignmentList <- readBamFile(bamFile = x, what=c("seq"), asMates = TRUE)
      lib.size.chrom <- countByChr(GAlignmentList)
      peCount <- addNamesAndMetadata(peCount, lib.size.chrom)
      peCount
    })
  } else
  {
    ## get count of fragments for each sliding window
    se <- lapply(reads, function(x)
    {
      seCount <- summarizeOverlaps(features = tileTargetRegions,
                                   reads = x,
                                   mode = mode,
                                   singleEnd = TRUE,
                                   inter.feature = FALSE, ...)
      GAlignmentList <- readBamFile(bamFile = x, what=c("seq"), asMates = FALSE)
      lib.size.chrom <- countByChr(GAlignmentList)
      
      seCount <- addNamesAndMetadata(seCount, lib.size.chrom)
      seCount
    })
  }
}