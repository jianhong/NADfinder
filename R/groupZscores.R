#' Calculate z-scores for each peak
#'
#' Detect peaks and calculate z-scores for each peak
#'
#' @param zscore A vector of numeric.
#' It is the z-scores of ratios for each window.
#' @return A data.frame with column names as "zscore", "group", "grp.zscore",
#' and "pvalue".
#'
#' @export
#' @importFrom stats pnorm
#' @examples
#' x <- seq_len(500)
#' a <- 2 * 2*pi/length(x)
#' y <- 20 * sin(x*a)
#' noise1 <- 20 * 1/10 * sin(x*a*10)
#' zscore <- y+noise1
#' groupZscores(zscore)
#'
groupZscores <- function(zscore)
{
    stopifnot(is.numeric(zscore) || is.integer(zscore))
    ## fix zscore == NA <- mean(i-1, i+1)
    if (any(is.na(zscore)))
    {
        zscore <- c(0, zscore, 0)
        i <- which(is.na(zscore))
        for (j in i) {
            ## avoide NA for continues NAs
            zscore[j] <- rowMeans(cbind(zscore[j - 1], zscore[j + 1]),
            na.rm = TRUE)
        }
        zscore <- zscore[-c(1, length(zscore))]
    }
    peaks <- peakdet(zscore)
    if (length(peaks$peakpos) == 0) {
        return(cbind(
            zscore,
            group = 1,
            grp.zscore = 0,
            pvalue = 1))
    }
    x <- seq.int(length(peaks$peakpos))
    times <- diff(c(0, peaks$valleypos, length(zscore)))
    if (length(times) != length(x))
    {
        x <- c(1, x + 1)
        peaks$peakpos <- c(peaks$peakpos, length(zscore))
    }
    group <- rep(x, times)
    grp.zscore <- zscore[peaks$peakpos[group]]
    pvalue <- 2 * pnorm(-abs(grp.zscore))
    cbind(zscore, group, grp.zscore, pvalue)
}
