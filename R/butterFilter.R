#' Low pass filter on ratios by butterworth filter
#'
#' The Butterworth filter is a type of signal processing filter designed to
#' have as flat a frequency response as possible in the passband.
#'
#' @param ratios A vector of numeric. It is the ratios of counts in each window.
#' @param N numeric(1) or integer(1). 
#' Critical frequencies of the low pass filter will be set as 1/N.
#' 1/N is a cutoff at 1/N-th of the Nyquist frequency.
#' Default suppose there are about 200 peaks in the inputs.
#'
#' @return A vector of numeric with same length of input ratios. 
#' The vector indicates smoothed ratios.
#'
#' @export
#' @importFrom signal butter filter
#' @importFrom stats loess.smooth
#' @examples
#' ratios <- runif(20000)
#' butterFilter(ratios)
#'
butterFilter <- function(ratios, N=ceiling(length(ratios)/200)){
    stopifnot(inherits(ratios, c("numeric", "integer")))
    stopifnot(length(N)==1)
    bf <- butter(2, 1/N, type="low")
    r2 <- as.numeric(filter(bf, ratios))
    W1 <- floor(N/2)+1
    W2 <- N-W1
    if(length(r2)>N && W2>0){
        y <- ratios[(length(r2)-W2):length(r2)]
        y[is.na(y)] <- min(y, na.rm=TRUE)
        c(r2[W1:length(r2)],
          loess.smooth(seq.int(W2+1),
                       y,
                       evaluation=W2+1)$y)[seq_along(r2)]
    }else{
        r2
    }
}
