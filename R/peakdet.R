#' Detect peak positions
#'
#' Detect the peak positions and valley positions leveraging
#' github::dgromer/peakdet
#'
#' @param y A numeric vector for searching peaks
#' @param delta A numeric vector of length 1, defining the minimum absolute changes required  
#' for local maximum or minimum detection when slope sign changes. If it is set to 0, 
#' the delta will be set to 1/10 of the range of y.
#' @param silence logical(1). 
#' If false, echo the delta value when delta is set as 0.
#'
#' @return A list with peakpos and valleypos. Both peakpos and valleypos are
#' numeric vectors storing the positions of peaks or valleys.
#'
#' @export
#'
#' @examples
#' y <- runif(200)
#' peakdet(y)
#' y <- sin(seq(0,20))
#' peakdet(y)

peakdet <- function(y, delta = 0, silence = TRUE)
{
    peakpos <- NULL
    valleypos <- NULL
    stopifnot(delta >= 0)
    if (delta == 0) 
    {
        ry <- quantile(y, c(0, .25, .75, 1))[2:3]
        delta <- diff(ry) / 10
        if (!silence) 
        {
            message("auto set delta = ", prettyNum(delta))
        }
    }
    minpos <- 1
    maxpos <- 1
    lookformax <- TRUE
    
    for (i in seq_along(y)) 
    {
        if (y[i] > y[maxpos]) 
        {
            maxpos <- i
        }
        if (y[i] < y[minpos]) 
        {
            minpos <- i
        }
        if (lookformax) 
        {
            if (y[i] < y[maxpos] - delta) 
            {
                peakpos <- c(peakpos, maxpos)
                minpos <- i
                lookformax <- FALSE
            }
        } else
        {
            if (y[i] > y[minpos] + delta) 
            {
                valleypos <- c(valleypos, minpos)
                maxpos <- i
                lookformax <- TRUE
            }
        }
    }
    
    list(peakpos = peakpos, valleypos = valleypos)
}
