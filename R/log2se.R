#' calculate the log2 transformed ratios for SummarizedExperiment class
#'
#' @description Calculate the log2 transformed ratios for nucleoleus vs genome.
#' pseudo-count will be used to avoid x/0 or log(0).
#'
#' @param se A \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
#' The output of \link{tileCount}.
#' @param nucleolusCols,genomeCols column Names of counts for nucleoleus 
#' and genome. They should be the column names in the assays of se. 
#' Ratios will be calculated as log2(transformed nucleolusCols/transformed genomeCols).
#' @param pseudocount default to 1, pseudo-count used to aviod x/0 or log(0).
#' @param transformation transformation type
#' @param chrom.level.lib indicating whether calculating CPM or odds using 
#' sequence depth of the whole genome or the corresponding chromosome
#' @export
#' @import SummarizedExperiment
#' @return A RangedSummarizedExperiment object with log2 transformed ratios. 
#' Assays will be named as nucleoleus, genome and ratio.
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays=list(counts=DataFrame(A=seq_len(3),
#'        B=rep(1, 3), C=rep(4, 3), D=rep(2, 3))),              
#'                   rowRanges=GRanges(c("chr1","chr1", "chr2"),
#'                       IRanges(c(1, 10, 20),
#'                             width=9)))
#' metadata(se)$lib.size.chrom <- data.frame( c(1000, 1000), c(2000, 2000), c(200,200), c(300,300))
#' colnames(metadata(se)$lib.size.chrom) <- c("A", "B", "C", "D")
#' rownames(metadata(se)$lib.size.chrom) <- c("chr1", "chr2")
#' as.data.frame(assays(log2se(se, nucleolusCols = c("A", "C"), genomeCols = c("B", "D"), 
#'     transformation = "log2Ratio")))
#' as.data.frame(assays(log2se(se, nucleolusCols = c("A", "C"), genomeCols = c("B", "D"),
#'     transformation = "log2CPMRatio")))
#' as.data.frame(assays(log2se(se, nucleolusCols = c("A", "C"), genomeCols = c("B", "D"),
#'     transformation = "log2OddsRatio")))
#' @author Jianhong Ou and Julie Zhu

log2se <- function(se, nucleolusCols, genomeCols, pseudocount = 1L, 
  transformation = c("log2OddsRatio", "log2CPMRatio", "log2Ratio"), chrom.level.lib = TRUE){
  stopifnot(inherits(se, "RangedSummarizedExperiment"))
  stopifnot("counts" %in% names(assays(se)))
  stopifnot(length(nucleolusCols)==length(genomeCols))
  stopifnot(length(nucleolusCols)>0)
  stopifnot(all(c(nucleolusCols, genomeCols) %in% 
                    colnames(assays(se)$counts)))
  transformation <- match.arg(transformation)
  asy <- assays(se)$counts
  if (transformation != "log2Ratio")
  {     stopifnot("lib.size.chrom" %in% names(metadata(se)))
        lib.size <- metadata(se)$lib.size.chrom 
        nucleolusCols.ind <- which(colnames(asy) %in% nucleolusCols)
        genomeCols.ind <- which(colnames(asy) %in% genomeCols)
  }

  nA <- names(nucleolusCols)
  ## nA will be used as the column names of log2 transformed ratios
  if(length(nA)==0){
    nA <- make.names(nucleolusCols, unique = TRUE)
  }
  
  log2ratios <- asy[, nucleolusCols]
  args <- list()

  if(!missing(pseudocount)){
    args$pseudo.count <- pseudocount
  }
  args$transformation <- transformation
  args$chrom.level.lib <- chrom.level.lib
  args$seqnames.A <- seqnames(se) 
  args$seqnames.B <-  args$seqnames.A


  if (transformation != "log2Ratio")
  {
      log2ratios <- do.call(cbind, mapply(function(a, b, lib.size.a, lib.size.b){
          args$A <- a
          args$B <- b
          args$lib.size.A <- cbind(rownames(lib.size), lib.size.a)
          args$lib.size.B <- cbind(rownames(lib.size), lib.size.b)
          #print(args$lib.size.A)
          do.call(transformData, args = args)
      }, data.frame(asy[, nucleolusCols, drop=FALSE]), 
     data.frame(asy[, genomeCols, drop=FALSE]), 
     data.frame(lib.size[, nucleolusCols.ind, drop=FALSE], stringsAsFactors = FALSE),
     data.frame(lib.size[, genomeCols.ind, drop=FALSE], stringsAsFactors = FALSE),
     SIMPLIFY = FALSE))
  }
  else
  {
      log2ratios <- do.call(cbind, mapply(function(a, b){
          args$A <- a
          args$B <- b
          do.call(transformData, args = args)
      }, data.frame(asy[, nucleolusCols, drop=FALSE]),
     data.frame(asy[, genomeCols, drop=FALSE]),
     SIMPLIFY = FALSE))
  }
  nucleoleus=asy[, nucleolusCols, drop=FALSE]
  genome=asy[, genomeCols, drop=FALSE]
  colnames(nucleoleus) <- colnames(genome) <- colnames(log2ratios) <- nA
  SummarizedExperiment(assays=list(nucleoleus=nucleoleus,
                                   genome=genome,
                                   ratio=log2ratios), 
                       rowRanges=rowRanges(se))
}
