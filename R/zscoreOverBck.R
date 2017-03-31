#' Z-scores over the background
#'
#' Calculate the z-scores over the lower percentage values.
#'
#' @param ratios A vector of numeric. It is the ratios of counts in each window.
#' @param backgroundPercentage numeric(1). Percentage of value for background.
#' @return A vector of numeric. Z-scores.
#' @export
#' @examples
#' r <- runif(200)
#' zscoreOverBck(r)
#'
zscoreOverBck <- function(ratios, backgroundPercentage=0.25){
    r <- quantile(ratios, probs = c(0, backgroundPercentage, 1))
    r <- ratios[ratios<=r[2]]
    pop_sd <- sd(r) * sqrt((length(r)-1)/length(r))
    pop_mean <- mean(r)
    z <- (ratios - pop_mean)/pop_sd
    z
}
