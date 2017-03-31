#' call peaks for ratios of repeats
#'
#' Use limma to call peaks for ratios of repeats
#'
#' @param se An object of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}
#' with assays of raw counts, ratios, background correct ratios,
#' smoothed ratios and z-scores. It should be an element of output of 
#' \link{smoothRatiosByChromosome}
#' @param backgroundCorrectionAssay character(1). Assays names 
#' for background correction ratios
#' @param normlizationMethod character(1) specifying the normalization
#' method to be used. Choices  are "none", "scale", "quantile" or "cyclicloess".
#' See \link[limma]{normalizeBetweenArrays} for details.
#' @param N numeric(1) or integer(1).
#' Critical frequencies of the low pass filter will be set as 1/N.
#' 1/N is a cutoff at 1/N-th of the Nyquist frequency.
#' Default 100.
#' @param cutoffAdjPvalue numeric(1). Cutoff adjust p-value.
#' @param countFilter numeric(1). Cutoff value for mean of raw reads count 
#' in each window.
#' @param ... Parameter not used.
#'
#' @import limma
#' @export
#' @return An object of GRanges of peak list with metadata "AveSig", "P.Value", 
#' and "adj.P.Val", where "AveSig" means average signals.
#' @examples
#'
#' data(triplicates.counts)
#' se <- triplicates.counts
#' gps <- c("26", "28", "29")
#' se <- log2se(se, 
#'              nucleosomeCols = paste0("N", gps, ".bam"),
#'              genomeCols = paste0("G", gps, ".bam"))
#' se<- smoothRatiosByChromosome(se, chr="chr18")
#' peaks <- callPeaks(se[[1]][10000:15000, ], 
#'                 cutoffAdjPvalue=0.05, countFilter=1000)
#'
callPeaks <- function(se, backgroundCorrectionAssay="bcRatio",
                      normlizationMethod="quantile",
                      N=100, cutoffAdjPvalue=0.05,
                      countFilter=1000, ...){
    normlizationMethod <-
        match.arg(normlizationMethod,
                  c("none", "scale", "quantile", "cyclicloess"))
    stopifnot(class(se)=="RangedSummarizedExperiment")
    if(any(!c("nucleosome", "genome", backgroundCorrectionAssay) %in% 
                      names(assays(se)))){
        stop("nucleosome", "genome", backgroundCorrectionAssay, 
             "should be the assays of se.")
    }
    stopifnot(ncol(assays(se)[[backgroundCorrectionAssay]])>=2)
    gr <- rowRanges(se)
    ## fixed window size
    stopifnot(all(width(gr)==width(gr)[1]))
    windowSize <- width(gr)[1]
    ### normalization among ratios
    bc <- as.data.frame(assays(se)[[backgroundCorrectionAssay]])
    bc.norm <- normalizeBetweenArrays(bc, method=normlizationMethod)
    ### call peaks
    bc.norm.rowMeans <- rowMeans(bc.norm)
    bc.smoothed <- butterFilter(bc.norm.rowMeans, N=N)
    peaks <- peakdet(bc.smoothed)
    if(length(peaks$peakpos)==0){
        peaks$peakpos <- ceiling(length(bc.smoothed)/2)
    }
    ## split the signals by peaks
    x <- seq.int(length(peaks$peakpos))
    times <- diff(c(0, peaks$valleypos, length(bc.smoothed)))
    if(length(times)!=length(x)){
        x <- c(1, x+1)
        peaks$peakpos <- c(peaks$peakpos, length(bc.smoothed))
    }
    group <- rep(x, times)
    if(length(group)!=length(bc.smoothed)){
        stop("The length of group is not identical with that of signals.",
             "Please report this bug.")
    }
    fit <- lmFit(bc.norm)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, number=nrow(fit2), sort.by="none")
    res$group <- group

    mcols(gr) <- DataFrame(res)
    keep <- gr$adj.P.Val < cutoffAdjPvalue
    gr <- gr[keep]
    se <- se[keep, ]
    gr <-
        gr[rowMeans(cbind(assays(se)[["nucleosome"]], 
                          assays(se)[["genome"]])) > countFilter]

    if(length(gr)==0){
        mcols(gr) <- DataFrame(AveSig=numeric(0),
                               P.Value=numeric(0),
                               adj.P.Val=numeric(0))
        colnames(mcols(gr)) <- c("AveSig", "P.Value", "adj.P.Val")
        return(gr)
    }

    gr <- split(gr, gr$group)
    gr.rd <- lapply(gr, function(.e) {
        ## peak summit
        if(length(.e)>10){
            sig <- loess.smooth(x=seq_along(.e), y=.e$AveExpr,
                                evaluation = length(.e))$y
        }else{
            sig <- .e$AveExpr
        }
        .idx <- which.max(sig)[1]
        if(!is.na(.idx)){
            .leftmin <- min(mcols(.e)[seq_len(.idx), "AveExpr"],
                            na.rm=TRUE)
            .rightmin <- min(mcols(.e)[.idx:length(.e), "AveExpr"],
                             na.rm=TRUE)
            peakHeight <- max(.e$AveExpr, na.rm=TRUE) -
                min(.leftmin, .rightmin, na.rm=TRUE)
            .z <- (.e$AveExpr - min(.leftmin, .rightmin, na.rm=TRUE)) >
                peakHeight/2
            .e <- .e[.z & .e$AveExpr>max(.leftmin, .rightmin, na.rm=TRUE)]
        }
        ra <- range(.e)
        if(length(.e)>0){
            ra$AveSig <- quantile(.e$AveExpr, probs=c(0, .75, 1), na.rm=TRUE)[2]
            ra$P.value <- quantile(.e$P.Value, probs=c(0, .25, 1),
                                   na.rm=TRUE)[2]
            ra$adj.P.Val <- quantile(.e$adj.P.Val, probs=c(0, .25, 1),
                                     na.rm=TRUE)[2]
        }else{
            mcols(ra) <- DataFrame(AveSig=numeric(0),
                                   P.Value=numeric(0),
                                   adj.P.Val=numeric(0))
            colnames(mcols(ra)) <- c("AveSig", "P.Value", "adj.P.Val")
        }
        ra
    })
    gr <- unlist(GRangesList(gr.rd))
    ## avoid overlaps
    wh <- ceiling(windowSize/2)
    gr <- gr[width(gr) > 2*wh]
    start(gr) <- start(gr) + wh
    end(gr) <- end(gr) - wh
    gr <- sort(gr)
    gr
}
