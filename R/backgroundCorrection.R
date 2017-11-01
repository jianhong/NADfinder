#' @title Correct ratios for background
#'
#' @description Background correct ratios of read counts for each window.
#'
#' @details This function implements the backgound correction methods of
#' algorithm for polynomial fitting. See details via
#' \code{\link[baseline]{baseline.modpolyfit}}. This function expects 
#' the trendency of decreasing of the ratios from 5' end to 3' end.
#'
#' @param ratios A vector of numeric.
#' It is the ratios of counts for each window.
#' @param degree Degree of polynomial. default 3.
#' @param ... parameters could be passed to
#'  \link[baseline]{baseline.modpolyfit}.
#'
#' @return A vector of numeric.
#' It is the background corrected ratios.
#'
#' @importFrom baseline baseline getCorrected
#' @export
#'
#' @examples
#' x <- runif(200)
#' background <- rep(c(20:1)/100, each=10)
#' backgroundCorrection(x)
#'
backgroundCorrection <- function(ratios, degree=3, ...){
    stopifnot(inherits(ratios, c("numeric", "integer")))
    if(any(is.na(ratios))){
        ratios[is.na(ratios)] <- 0
    }
    if(any(is.infinite(ratios))){
        ratios[is.infinite(ratios)] <- 0
    }
    if(length(unique(ratios))>1){
        bc <- baseline(t(ratios), method="modpolyfit", degree = degree, ...)
        as.numeric(getCorrected(bc))
    }else{
        ratios
    }
}
