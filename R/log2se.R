#' calculate the log2 transformed ratios for SummarizedExperiment class
#'
#' @description Calculate the log2 transformed ratios for nucleosome vs genome.
#' pseudo-count will be used to avoid x/0.
#'
#' @param se A \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
#' The output of \link{tileCount}.
#' @param nucleosomeCols,genomeCols column Names of counts for nucleosome 
#' and genome. They should be the column names in the assays of se. 
#' Ratios will be calculated as log2(nucleosomeCols/genomeCols).
#' @param pseudocount pseudo-count will be used to aviod x/0 by x/pseudocount.
#' If it is not set, pseudo-count will be the minimal count except 0 of inputs.
#' @export
#' @import SummarizedExperiment
#' @return A RangedSummarizedExperiment object with log2 transformed ratios. 
#' Assays will be named as nucleosome, genome and ratio.
#' @examples
#' library(SummarizedExperiment)
#' se <- 
#'  SummarizedExperiment(assays=list(counts=DataFrame(A=seq_len(3), 
#'                                                    B=rep(1, 3))), 
#'                       rowRanges=GRanges("chr1", 
#'                                          IRanges(c(1, 10, 20), 
#'                                                  width=9)))
#' log2se(se, "A", "B")
#'
log2se <- function(se, nucleosomeCols, genomeCols, pseudocount){
  stopifnot(inherits(se, "RangedSummarizedExperiment"))
  stopifnot("counts" %in% names(assays(se)))
  stopifnot(length(nucleosomeCols)==length(genomeCols))
  stopifnot(length(nucleosomeCols)>0)
  stopifnot(all(c(nucleosomeCols, genomeCols) %in% 
                    colnames(assays(se)$counts)))
  asy <- assays(se)$counts
  nA <- names(nucleosomeCols) 
  ## nA will be used as the column names of log2 transform ratios
  if(length(nA)==0){
    nA <- make.names(nucleosomeCols, unique = TRUE)
  }
  log2ratios <- asy[, nucleosomeCols]
  args <- list()
  if(!missing(pseudocount)){
    args$pseudocount <- pseudocount
  }
  log2ratios <- do.call(cbind, mapply(function(a, b){
    args$A <- a
    args$B <- b
    do.call(log2ratio, args = args)
  }, data.frame(asy[, nucleosomeCols, drop=FALSE]), 
     data.frame(asy[, genomeCols, drop=FALSE]), 
     SIMPLIFY = FALSE))
  nucleosome=asy[, nucleosomeCols, drop=FALSE]
  genome=asy[, genomeCols, drop=FALSE]
  colnames(nucleosome) <- colnames(genome) <- colnames(log2ratios) <- nA
  SummarizedExperiment(assays=list(nucleosome=nucleosome,
                                   genome=genome,
                                   ratio=log2ratios), 
                       rowRanges=rowRanges(se))
}