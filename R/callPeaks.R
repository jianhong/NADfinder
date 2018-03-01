#' 
#' Call peaks using transformed, background corrected, and smoothed ratios with biological replicates
#'
#' Use limma to calculate p-values for NADs
#'
#' By default, use the mean smoothed ratio for each peak region to calculate p-values 
#'
#' @param se An object of 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' with assays of raw counts, tranformed ratios, background corrected ratios,
#' smoothed ratios and z-scores. It should be an element of output of 
#' \link{smoothRatiosByChromosome}
#' @param backgroundCorrectedAssay character(1). Assays names 
#' for background corrected log2-transformed ratios, CPMRatios or OddRatios.
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
#' @param method.combineP A method used to combine P-values. Default meansig
#' @param ... Parameter not used.
#'
#' @import limma
#' @import S4Vectors
#' @import EmpiricalBrownsMethod  
#' @import metap
#' @export
#' @return An object of GRanges of peak list with metadata "AveSig", "P.Value", 
#' and "adj.P.Val", where "AveSig" means average signal such as average log2OddsRatio, log2CPMRatio or log2Ratio.
#' @author Jianhong Ou, Haibo Liu and Julie Zhu

#' @examples
#'
#' data(triplicate.count)
#' se <- triplicate.count
#' se <- log2se(se, transformation = "log2CPMRatio",
#'              nucleoleusCols = c("N18.subsampled.srt-2.bam",
#'              "N18.subsampled.srt-3.bam",
#'              "N18.subsampled.srt.bam"),
#'              genomeCols = c("G18.subsampled.srt-2.bam",
#'              "G18.subsampled.srt-3.bam",
#'              "G18.subsampled.srt.bam"))
#' se<- smoothRatiosByChromosome(se, chr="chr18")
#' peaks <- callPeaks(se[[1]][10000:15000, ], 
#'                 cutoffAdjPvalue=0.001, countFilter=1000)
#'
callPeaks <- function(se,
                      backgroundCorrectedAssay = "bcRatio",
                      normlizationMethod = "quantile",
                      N = 100,
                      cutoffAdjPvalue = 0.0001,
                      countFilter = 1000,
                      method.combineP = "meansig",
                      ...) {
    stopifnot(class(se) == "RangedSummarizedExperiment")
    stopifnot(ncol(assays(se)[[backgroundCorrectedAssay]]) >= 2)
    if (any(!c("nucleoleus", "genome", backgroundCorrectedAssay) %in%
            names(assays(se))))
    {
        stop(
            "nucleoleus",
            "genome",
            backgroundCorrectedAssay,
            "should be the assays of se."
        )
    }
    
    normlizationMethod <-
        match.arg(normlizationMethod,
                  c("none", "scale", "quantile", "cyclicloess"))
    method.combineP <-
        match.arg(method.combineP,
                  c(
                      "Browns",
                      "minimump",
                      "logitp",
                      "Fishers",
                      "sumz",
                      "meansig"
                  ))
    
    gr <- rowRanges(se)
    ## fixed window size
    stopifnot(all(width(gr) == width(gr)[1]))
    windowSize <- width(gr)[1]
    ## normalization among ratios
    bc <- as.data.frame(assays(se)[[backgroundCorrectedAssay]])
    bc.norm <- normalizeBetweenArrays(bc, method = normlizationMethod)
    bc.norm.rowMeans <- rowMeans(bc.norm)
    
    
    
    #### Why are so many smoothing steps needed? (Haibo comments on it)
    ## Jianhong, do you want to add a parameter as smooth.method (loess, buttfilter)
    bc.smoothed <- butterFilter(bc.norm.rowMeans, N = N)
    peaks <- peakdet(bc.smoothed)
    if (length(peaks$peakpos) == 0)
    {
        peaks$peakpos <- which(bc.smoothed == max(bc.smoothed))
    }
    ## split the signals by peaks
    x <- seq.int(length(peaks$peakpos))
    times <- diff(c(0, peaks$valleypos, length(bc.smoothed)))
    if (length(times) != length(x))
    {
        x <- c(1, x + 1)
        peaks$peakpos <- c(peaks$peakpos, length(bc.smoothed))
    }
    ## assume points between previous valley and this valley belong to the group of this valley
    group <- rep(x, times)
    if (length(group) != length(bc.smoothed))
    {
        stop(
            "The length of group is not identical with that of signals.",
            "Please report this bug.")
    }
    
    fit <- lmFit(bc.norm)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, number = nrow(fit2), sort.by = "none")
    res <- cbind(bc.norm, res)
    
    #### group windows by valley plus points after previous valley
    res$group <- group
    
    mcols(gr) <- DataFrame(res)
    gr.all <- gr
    
    keep <- gr$adj.P.Val < cutoffAdjPvalue
    gr <- gr[keep]
    se <- se[keep,]
    ### add cpm filter or count filter to retain windows with at least half of the samples
    ### the following filter favors dataset with more replicates
    
    #### This combined data.frame has already exist in tileCount() output
    gr <- gr[rowMeans(cbind(assays(se)[["nucleoleus"]],
                            assays(se)[["genome"]])) > countFilter]
    
    if (length(gr) == 0)
    {
        mcols(gr) <- DataFrame(
            AveSig = numeric(0),
            P.Value = numeric(0),
            adj.P.Val = numeric(0)
        )
        colnames(mcols(gr)) <- c("AveSig", "P.Value", "adj.P.Val")
        return(gr)
    }
    
    gr <- split(gr, gr$group)
    
    
    gr.rd <- lapply(gr, function(.e) {
        ## peak summit
        ## Jianhong, do we need to smooth the data again? (you used butterFilter
        ## to smooth the data previously)
        if (length(.e) > 10)
        {
            sig <- loess.smooth(
                x = seq_along(.e),
                y = .e$AveExpr,
                evaluation = length(.e)
            )$y
        } else
        {
            sig <- .e$AveExpr
        }
        .idx <- which.max(sig)[1]
        if (!is.na(.idx))
        {
            .leftmin <- min(mcols(.e)[seq_len(.idx), "AveExpr"],
                            na.rm = TRUE)
            .rightmin <- min(mcols(.e)[.idx:length(.e), "AveExpr"],
                             na.rm = TRUE)
            peakHeight <- max(.e$AveExpr, na.rm = TRUE) -
                min(.leftmin, .rightmin, na.rm = TRUE)
            .z <-
                (.e$AveExpr - min(.leftmin, .rightmin, na.rm = TRUE)) >
                peakHeight / 2
            .e <-
                .e[.z & .e$AveExpr > max(.leftmin, .rightmin, na.rm = TRUE)]
        }
        ra <- range(.e)
        
        all.windows.in.ra <-
            gr.all[seqnames(gr.all) == seqnames(ra) &
                       (start(gr.all) + windowSize / 2) >= start(ra) &
                       (end(gr.all) - windowSize / 2) <= end(ra),]
        
        if (length(.e) > 0)
        {
            ra$AveSig <- mean(all.windows.in.ra$AveExpr)
            #ra$AveSig <- quantile(all.windows.in.ra$AveExpr, probs=c(0, .75, 1), na.rm=TRUE)[2]
            if (method.combineP  == "Browns")
            {
                ## sample size too small here for using this methods?
                ra$P.value <- empiricalBrownsMethod(
                    data_matrix =
                        mcols(all.windows.in.ra)[, 1:dim(bc.norm)[2]],
                    p_values = all.windows.in.ra$P.Value,
                    extra_info = FALSE)
                ra$adj.P.Val <-
                    empiricalBrownsMethod(
                        data_matrix =
                            mcols(all.windows.in.ra)[, 1:dim(bc.norm)[2]],
                        p_values = all.windows.in.ra$adj.P.Val,
                        extra_info = FALSE)
            } else if (method.combineP == "Fishers")
            {
                ra$P.value <- sumlog(all.windows.in.ra$P.Value)$p
                ra$adj.P.Val <- sumlog(all.windows.in.ra$adj.P.Val)$p
            } else if (method.combineP == "logitp")
            {
                ra$P.value <- logitp(all.windows.in.ra$P.Value)$p
                ra$adj.P.Val <- logitp(all.windows.in.ra$adj.P.Val)$p
            } else if (method.combineP == "minimump")
            {
                ra$P.value <- minimump(all.windows.in.ra$P.Value)$p
                ra$adj.P.Val <- minimump(all.windows.in.ra$adj.P.Val)$p
            } else if (method.combineP == "sumz")
            {
                print("OK")
                ra$P.value <- sumz(all.windows.in.ra$P.Value)$p
                ra$adj.P.Val <- sumz(all.windows.in.ra$adj.P.Val)$p
            }
            
            ## This is common to all combineP methods, so remove redundancy
            mcols(ra) <- cbind(DataFrame(t(colMeans(as.matrix(
                mcols(all.windows.in.ra)[, 1:dim(bc.norm)[2]])))), mcols(ra))
        } else
        {
            mcols(ra) <- DataFrame(
                AveSig = numeric(0),
                P.Value = numeric(1),
                adj.P.Val = numeric(1)
            )
            colnames(mcols(ra)) <-
                c("AveSig", "P.Value", "adj.P.Val")
            mcols(ra) <- cbind(t(colMeans(as.matrix(
                mcols(all.windows.in.ra)[, 1:dim(bc.norm)[2]]
            ))), mcols(ra))
        }
        ra
    })
    
    gr <- unlist(GRangesList(gr.rd))
    
    if (method.combineP == "meansig")
    {
        fit <- lmFit(mcols(gr)[, 1:dim(bc.norm)[2]])
        fit2 <- eBayes(fit)
        res <- topTable(fit2, number = nrow(fit2), sort.by = "none")
        mcols(gr)$P.Value =  res$P.Value
        mcols(gr)$adj.P.Val = res$adj.P.Val
        mcols(gr)$AveSig  = res$logFC
    }
    
    
    ## avoid overlaps
    wh <- ceiling(windowSize / 2)
    #gr <- gr[width(gr) > 2*wh]
    start(gr[width(gr) > wh]) <- start(gr[width(gr) > wh]) + wh
    end(gr[width(gr) > wh]) <- end(gr[width(gr) > wh]) - wh
    gr <- sort(gr)
    gr
}
