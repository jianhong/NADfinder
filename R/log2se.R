#' Calculate the log2-transformed ratios of signal stored in a SummarizedExperiment object
#'
#' @description Calculate the log2 transformed ratios of signals for nucleoleus-associated 
#' DNA vs the whole genome DNA reference. Pseudo-count will be used to avoid zero division or log(0).
#'
#' @param se A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} object.
#' The output of \link{tileCount}. Columnnames of counts data (fragment counts in each tile) 
#' and metadata (total fragment counts per chromosome) are bam file names.
#' @param nucleoluesCols,genomeCols Column names of counts for nucleoleus-associated DNA
#' and the whole genome DNA. They should be the column names in the assays slot of 
#' an RangedSummarizedExperiment object. Ratios will be calculated as log2 (transformed-
#' nucleoluesCols/transformed genomeCols).
#' @param pseudocount Default to 1, pseudo-count used to aviod zero division or log(0).
#' @param transformation Singal transformation method.
#' @param chrom.level.lib Logical (1) indicating whether calculating CPM or odds using 
#' sequence depth of the whole genome or per chromosome.
#' @export
#' @import SummarizedExperiment
#' @return A RangedSummarizedExperiment object with log2 transformed ratios. 
#' Assays will be named as nucleoleus, genome and ratios.
#' @author Jianhong Ou, Haibo Liu and Julie Zhu
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays=list(counts=DataFrame(A=seq_len(3),
#'        B=rep(1, 3), C=rep(4, 3), D=rep(2, 3))),              
#'                   rowRanges=GRanges("chr1",
#'                       IRanges(c(1, 10, 20),
#'                             width=9)))
#' metadata(se)$lib.size.chrom <- data.frame(c("chr1"), c(1000), c(2000), c(200), c(300))
#' as.data.frame(assays(log2se(se, nucleoleusCols = c("A", "C"), genomeCols = c("B", "D"), 
#'     transformation = "log2Ratio")))
#' as.data.frame(assays(log2se(se, nucleoleusCols = c("A", "C"), genomeCols = c("B", "D"),
#'     transformation = "log2CPMRatio")))
#' as.data.frame(assays(log2se(se, nucleoleusCols = c("A", "C"), genomeCols = c("B", "D"),
#'     transformation = "log2OddsRatio")))


log2se <- function(se,
                   nucleoleusCols,
                   genomeCols,
                   pseudocount = 1L,
                   transformation = c("log2OddsRatio", "log2CPMRatio", "log2Ratio"),
                   chrom.level.lib = TRUE)
{
    stopifnot(inherits(se, "RangedSummarizedExperiment"))
    stopifnot("counts" %in% names(assays(se)))
    stopifnot(length(nucleoleusCols) == length(genomeCols))
    stopifnot(length(nucleoleusCols) > 0)
    stopifnot(all(c(nucleoleusCols, genomeCols) %in% colnames(assays(se)$counts)))
    transformation <- match.arg(transformation)
    
    asy <- assays(se)$counts
    nucleoleusCols.ind <- which(colnames(asy) %in% nucleoleusCols)
    genomeCols.ind <- which(colnames(asy) %in% genomeCols)
    
    ## nA will be used as the column names of log2 transformed ratios
    nA <- make.names(nucleoleusCols, unique = TRUE)
    nucleoleus = data.frame(asy[, nucleoleusCols, drop = FALSE])
    genome = data.frame(asy[, genomeCols, drop = FALSE])
    chrNames <- rownames(metadata(se)$lib.size.chrom)
    
    args <- list()
    if (!missing(pseudocount)) {
        args$pseudocount <- pseudocount
    }
    args$transformation <- transformation
    args$chrom.level.lib <- chrom.level.lib
    args$seqnames.A <-  args$seqnames.B <- seqnames(se)
   
    if (transformation != "log2Ratio")
    {
        stopifnot("lib.size.chrom" %in% names(metadata(se)))
        lib.size <- metadata(se)$lib.size.chrom
        
        log2ratios <- do.call(cbind, mapply(
                        function(a, b, lib.size.a, lib.size.b) {
                        args$A <- a
                        args$B <- b
                        args$lib.size.A <- as.data.frame(lib.size.a, row.names = chrNames)
                        args$lib.size.B <- as.data.frame(lib.size.b, row.names = chrNames)
                        
                        do.call(transformData, args = args)
                    },
                    nucleoleus, genome,
                    
                    ## my lib.size.chrom is a data frame with chromosome names as rownames, 
                    ## and sample names as colnames so +1 for the index is not necessary
                    data.frame(lib.size[, nucleoleusCols.ind, drop = FALSE], stringsAsFactors = FALSE),
                    data.frame(lib.size[, genomeCols.ind, drop = FALSE], stringsAsFactors = FALSE),
                    SIMPLIFY = FALSE))
    } else
    {
        log2ratios <- do.call(cbind, mapply(
                        function(a, b) {
                                    args$A <- a
                                    args$B <- b
                                    do.call(transformData, args = args)
                                  },
                        nucleoleus, genome,
                        SIMPLIFY = FALSE))
    }
    colnames(nucleoleus) <- colnames(genome) <- colnames(log2ratios) <- nA
    SummarizedExperiment(assays = list(nucleoleus = nucleoleus, genome = genome, 
                                       ratio = log2ratios), rowRanges = rowRanges(se))
}
