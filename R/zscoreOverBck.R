#' Z-scores over the background
#'
#' Calculate the z-scores over the background distribution.
#'
#' @param ratios A numeric vector containing the transformed, background corrected and smoothed ratios in each window.
#' @param backgroundPercentage numeric(1). Low percentile for background distribution.
#' @return A vector of numeric. Z-scores.
#' @export
#' @examples
#' r <- runif(200)
#' zscoreOverBck(r)
#' @author Jianhong Ou and Julie Zhu

zscoreOverBck <- function(ratios, backgroundPercentage=0.25){
    r <- quantile(ratios, probs = backgroundPercentage)
    bg.ratios <- ratios[ratios<=r]
    pop.sd <- sd(bg.ratios) * sqrt(length(bg.ratios) - 1)/sqrt(length(bg.ratios)) 
    pop.mean <- mean(bg.ratios)
    z <- (ratios - pop.mean)/pop.sd
    z
}
