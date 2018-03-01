#' Count reads overlapping genomic ranges
#'
#' Count reads overlapping a set of genimc features represented as
#' genomic ranges. This function does not work for parallel.
#' @param features A object of \link[GenomicRanges]{GRanges} representing the
#' feature regions to be counted.
#' @param reads An object that represents the data to be counted. See
#' \link[GenomicAlignments]{summarizeOverlaps}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory. For paired end reads, 
#' @param ignore.strand logical(1). ignore strand?
#' @param inter.feature not used. This parameter is required by
#' \link[GenomicAlignments]{summarizeOverlaps}.
#' @export
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
#' tileCount extends \link[GenomicAlignments]{summarizeOverlaps} by finding coverage for
#' each fixed window in the whole genome
#'
#' @param reads A \link[GenomicRanges]{GRanges},
#' \link[GenomicRanges]{GRangesList} (should be one read per list element),
#' \link[GenomicAlignments]{GAlignments},
#' \link[GenomicAlignments]{GAlignmentsList},
#' \link[GenomicAlignments]{GAlignmentPairs} or
#' \link[Rsamtools]{BamFileList} object that represents the data to be
#' counted by \code{\link[GenomicAlignments]{summarizeOverlaps}}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory.
#' @param genome A BSgenome object from/on which to get/set the sequence and metadata information.
#' @param windowSize numeric(1) or integer(1). Size of the windows.
#' @param step numeric(1) or integer(1). Step of generating silding windows.
#' @param mode One of the pre-defined count methods.
#' @param excludeChrs A vector of string: chromosomes/scaffolds of no interest for NAD analysis.
#' see \link[GenomicAlignments]{summarizeOverlaps}.
#' default is countByOverlaps, alia of countOverlaps(features, reads, ignore.strand=ignore.strand)
#' @param dataOverSamples logical(1). Data over several samples when use
#' \link[GenomicRanges]{GRangesList} as input.
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
#' @importFrom ATACseqQC readBamFile
#' @export
#' @examples
#' \dontrun{
#' fls <- list.files(system.file("extdata", package="NADfinder"),
#' recursive=FALSE, pattern="*bam$", full=TRUE)
#' names(fls) <- basename(fls)
#' if (!require(BSgenome.Mmusculus.UCSC.mm10))
#' {
#'     source("https://bioconductor.org/biocLite.R")
#'     biocLite("BSgenome.Mmusculus.UCSC.mm10"))
#'     library(BSgenome.Mmusculus.UCSC.mm10)
#' }
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


tileCount<- function(reads,
                      genome,
                      excludeChrs = c("chrM", "M", "Mt", "MT"),
                      windowSize = 1e5L,
                      step = 1e4L,
                      mode = IntersectionNotStrict,
                      dataOverSamples = FALSE,
                      ...) 
{
    stopifnot(all(grepl(".bam$", reads)), all(file.exists(paste0(reads, ".bai"))))
    stopifnot(class(genome) == "BSgenome")
    stopifnot(windowSize %% 1 == 0, step %% 1 ==0, windowSize > 0, step > 0, step < windowSize)

    
    ## the scanbam function canonly read a single bam file not a list of bamfiles at a time
    ## so for a list bamfiles, I modify the code as follows
    ## aln is a list of list of lists
    isPE <- testPairedEndBam(reads[1])
    
    
    targetRegions <- as(seqinfo(genome), "GRanges")

    param <- ScanBamParam(what=c("rname", "qname"))
    aln <- lapply(reads, scanBam, param = param)
    lib.size.chrom <- do.call(cbind, lapply(1:length(aln),function(i)
    {
            ref_qnames <- unique(data.frame(rnames = aln[[i]][[1]]$rname, qnames = aln[[i]][[1]]$qname))
            if(i==1)
            {
                countByChr = as.data.frame(BiocGenerics::table(ref_qnames$rname))
                rownames(countByChr) <- countByChr[,1]
                countByChr <- countByChr[, 2, drop =FALSE]
                
            } else 
            {
                countByChr = as.data.frame(BiocGenerics::table(ref_qnames$rname))[, 2, drop =FALSE]  
            }
    
            ## using the bamfile names to name the column derived from a give bamfiles
            colnames(countByChr) <- reads[i]
            countByChr
    }))
    
    ## remove excludeChr from lib.size.chrom
    lib.size.chrom <- lib.size.chrom[!rownames(lib.size.chrom) %in% excludeChrs, ]
    
    ## filtering GRanges to keep only those chromosomal scaffolds that are in the BAM file
    targetRegions <- targetRegions[seqnames(targetRegions) %in% rownames(lib.size.chrom)]
    tileTargetRegions <-slidingWindows(x = targetRegions,
                                       width = windowSize,
                                       step = step)
    tileTargetRegionsLen <- elementNROWS(tileTargetRegions)
    tileTargetRegions <- unlist(tileTargetRegions)
    mcols(tileTargetRegions)$oid <- rep(seq_along(targetRegions),tileTargetRegionsLen)
    
    if (isPE)
    {
        rse <- summarizeOverlaps(features = tileTargetRegions,
                                reads = reads,
                                mode = mode,
                                singleEnd = FALSE,
                                fragments = TRUE,
                                inter.feature = FALSE)
        
    } else
    {
        rse <- summarizeOverlaps(features = tileTargetRegions,
                                reads = reads,
                                mode = mode,
                                singleEnd = TRUE,
                                inter.feature = FALSE)
    }
    names(assays(rse)) <- "counts"
    metadata(rse)$lib.size.chrom <- lib.size.chrom
    rse
}