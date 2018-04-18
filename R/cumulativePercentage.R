#' Plot the cumulative percentage of tag allocation
#' 
#' Plot the difference between the cumulative percentage of tag allocation in 
#' paired samples.
#' 
#' @param se An object of 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' with assays of raw counts, transfomred ratios, background correct ratios,
#' smoothed ratios and z-scores. It should be an element of the output of 
#' \link{smoothRatiosByChromosome}.
#' @param binWidth numeric(1) or integer(1). The width of each bin.
#' @param backgroundCorrectedAssay character(1). Assays names 
#' for background correction ratios.
#' @param ... Parameter not used.
#' @import SummarizedExperiment
#' @import GenomicAlignments
#' @importFrom IRanges Views viewApply tile viewMeans
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics which.min
#' @importFrom graphics abline axis legend matlines par plot
#' @export
#' @return A list of data.frame with the cumulative percentages.
#' @references Normalization, bias correction, and peak calling for ChIP-seq
#' Aaron Diaz, Kiyoub Park, Daniel A. Lim, Jun S. Song
#' Stat Appl Genet Mol Biol. Author manuscript; 
#' available in PMC 2012 May 3.Published in final edited form as: 
#' Stat Appl Genet Mol Biol. 2012 Mar 31; 11(3): 10.1515/1544-6115.1750
#'  /j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml. 
#'  Published online 2012 Mar 31.  doi: 10.1515/1544-6115.1750
#'  PMCID: PMC3342857
#' @examples 

#' library(SummarizedExperiment)
#' data(triplicate.count)
#' se <- triplicate.count
#' se <- log2se(se, transformation = "log2CPMRatio",
#'              nucleolusCols = c("N18.subsampled.srt-2.bam",
#'              "N18.subsampled.srt-3.bam",
#'              "N18.subsampled.srt.bam"),
#'              genomeCols = c("G18.subsampled.srt-2.bam",
#'              "G18.subsampled.srt-3.bam",
#'              "G18.subsampled.srt.bam"))
#' se <- smoothRatiosByChromosome(se, chr="chr18")
#' cumulativePercentage(se[["chr18"]])

cumulativePercentage <- function(se,
                                 binWidth = 1e5,
                                 backgroundCorrectedAssay = "bcRatio",...) 
{
    stopifnot(is(se, "RangedSummarizedExperiment"))
    assayName <-
        c("nucleolus", "genome", backgroundCorrectedAssay)
    if (any(!assayName %in% names(assays(se)))) 
    {
        stop(
            "nucleolus",
            "genome",
            backgroundCorrectedAssay,
            "should be the assays of se.")
    }
    ## resample
    sampleName <-
        unique(do.call(rbind, lapply(assays(se), colnames)))
    if (nrow(sampleName) != 1) 
    {
        stop("The column names of assays in se are not identical.")
    }
    sampleName <- sampleName[1, , drop = TRUE]
    seqL <- ranges(range(rowRanges(se)))
    if (length(seqL) != 1) 
    {
        stop("One chromosome only.")
    }
    features <- tile(seqL, width = binWidth)[[1]]
    sigBin <- lapply(sampleName, function(.n) {
        sig <- sapply(assayName[-3],
                      exportSignals,
                      dat = se,
                      colName = .n)
        sig <- lapply(sig, function(.ele) {
            .ele <- .ele[[1]]
            v <- Views(.ele, features)
            viewMeans(v, na.rm = TRUE)
        })
        sig <- do.call(cbind, sig)
        sig[order(sig[, "nucleolus"]),]
    })
    sigCumsum <- lapply(sigBin, function(.ele) {
        .ele <- apply(.ele, 2, cumsum)
        .ele <- cbind(Rank = seq_len(nrow(.ele)), .ele)
        sweep(.ele,
              MARGIN = 2,
              STATS = .ele[nrow(.ele), ],
              FUN = `/`)
    })
    sigEnrichStart <- sapply(sigBin, function(.ele) {
        .r <- (.ele[, "nucleolus"] + 1) / (.ele[, "genome"] + 1)
        ## split .r into two parts, background and enriched
        .x <- cumsum(.r)
        .y <- sum(.r) - .x
        .l <- length(.x)
        .mx <- .x / seq_len(.l)
        .my <- .y / (.l - seq_len(.l) + 1)
        .R <- numeric(.l)
        for (i in seq_len(.l - 1)) 
        {
            .R[i] <- sum((.r[seq_len(i)] - .mx[i]) ^ 2) +
                sum((.r[(i + 1):.l] - .my[i]) ^ 2)
        }
        which.min(.R)[1] / .l})
    pin <- par("pin")
    if (pin[2] > 0) 
    {
        ratio <- 2 ^ round(diff(log2(pin)))
        n <- length(sampleName)
        ncol <- ceiling(sqrt(n / ratio))
        nrow <- ceiling(n / ncol)
        op <- par(mfrow = c(nrow, ncol), pty = "s")
        on.exit(par(op))
        for (i in seq_len(n)) 
        {
            ## plot
            plot(
                c(0, 1),
                c(0, 1),
                type = "n",
                xlab = "% of bins",
                ylab = "% of tags",
                main = sampleName[i]
            )
            matlines(x = sigCumsum[[i]][, 1],
                     y = sigCumsum[[i]][, -1])
            zero <- which(sigCumsum[[i]][, "nucleolus"] > 1 / binWidth)
            if (length(zero) > 0) 
            {
                x.tick <- sigCumsum[[i]][zero[1], 1]
                abline(v = x.tick,
                       col = "yellowgreen",
                       lty = 3)
                if (sigEnrichStart[i] < 1) 
                {
                    abline(v = sigEnrichStart[i],
                           col = "violetred",
                           lty = 3)
                    x.tick <- c(x.tick, sigEnrichStart[i])
                }
                axis(3,
                     at = x.tick,
                     labels = formatC(x.tick, digits = 2))
            }
            legend(
                "topleft",
                legend = colnames(sigCumsum[[i]])[-1],
                col = seq_len(6),
                lty = seq_len(5),
                pch = NA,
                box.col = NA
            )
        }
    }
    return(invisible(sigCumsum))
}
