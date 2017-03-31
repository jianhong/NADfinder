#' calculate the log2 transformed ratios
#'
#' @description calculate the log2 transformed ratios for nucleosome vs genome.
#' pseudo-count will be used to avoid x/0.
#'
#' @param A,B counts for nucleosome and genome. They should be
#' numeric vectors with identical length.
#' @param pseudocount pseudo-count will be used to aviod x/0 by x/pseudocount.
#' If it is not set, pseudo-count will be the minimal count except 0 of inputs.
#' @export
#' @return A vector of numeric of log2 transformed ratios.
#' @examples
#' log2ratio(seq_len(10), 10:1)
#'

log2ratio <- function(A, B, pseudocount){
    stopifnot(length(A)==length(B))
    stopifnot(inherits(A, c("numeric", "integer")))
    stopifnot(inherits(B, c("numeric", "integer")))
    stopifnot(all(A>=0))
    stopifnot(all(B>=0))
    if(missing(pseudocount)){
        cnts <- c(A, B)
        cnts <- cnts[cnts>0]
        pseudocount <- min(cnts, na.rm = TRUE)
    }
    r <- log2(A) - log2(B)
    r[is.infinite(r)] <- log2(A[is.infinite(r)]+pseudocount) -
        log2(B[is.infinite(r)]+pseudocount)
    r[is.na(r)] <- 0
    r
}
