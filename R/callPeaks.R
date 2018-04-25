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
#' @param normalization.method character(1) specifying the normalization
#' method to be used. Choices  are "none", "scale", "quantile" or "cyclicloess".
#' See \link[limma]{normalizeBetweenArrays} for details.
#' @param N numeric(1) or integer(1).
#' The number of  neighboring windows used for loess smoothing or the inverse of
#' the critical frequencies of the low pass filter for butterworth filter.
#' 1/N is a cutoff at 1/N-th of the Nyquist frequency.
#' Default 100.
#' @param cutoffAdjPvalue numeric(1). Cutoff adjust p-value.
#' @param countFilter numeric(1). Cutoff value for mean of raw reads count
#' in each window.
#' @param combineP.method A method used to combine P-values. Default meansig
#' @param smooth.method A method used to smooth the ratios. Choices are "loess",
#' "none" and "butterworthfilter".
#' @param lfc the minimum log2-fold-change that is considered scientifically meaningful
#' @param ... Parameter not used.
#' 
#' @import SummarizedExperiment
#' @import limma
#' @import S4Vectors
#' @import EmpiricalBrownsMethod 
#' @import metap
#' @importFrom stats loess median predict
#' @export
#' @return An object of GRanges of peak list with metadata "AveSig", "P.Value",
#' and "adj.P.Val", where "AveSig" means average signal such as average log2OddsRatio, log2CPMRatio or log2Ratio.
#' @author Jianhong Ou, Haibo Liu and Julie Zhu

#' @examples
#'
#' data(triplicate.count)
#' se <- triplicate.count
#' se <- log2se(se, transformation = "log2CPMRatio",
#'              nucleolusCols = c("N18.subsampled.srt-2.bam",
#'              "N18.subsampled.srt-3.bam",
#'              "N18.subsampled.srt.bam"),
#'              genomeCols = c("G18.subsampled.srt-2.bam",
#'              "G18.subsampled.srt-3.bam",
#'              "G18.subsampled.srt.bam"))
#' se<- smoothRatiosByChromosome(se, chr="chr18")
#' #add some variability to the data since the triplicate.count data was created using one sample only
#' assays(se[[1]])$bcRatio[,2] <- assays(se[[1]])$bcRatio[,2] + 0.3
#' assays(se[[1]])$bcRatio[,3] <- assays(se[[1]])$bcRatio[,3] - 0.3
#' peaks <- callPeaks(se[[1]],
#'                 cutoffAdjPvalue=0.001, countFilter=10)
#'


#### helper function to merge the continuous bins in the peaks

callPeaks <- function(se,
                      backgroundCorrectedAssay = "bcRatio",
                      normalization.method = "quantile",
                      N = 100,
                      cutoffAdjPvalue = 0.0001,
                      countFilter = 1000,
                      combineP.method = "meansig",
                      smooth.method = "loess",
                      lfc = log2(1.5),
                      ...) {
  stopifnot(is(se, "RangedSummarizedExperiment"))
  stopifnot(ncol(assays(se)[[backgroundCorrectedAssay]]) >= 2)
  if (any(!c("nucleolus", "genome", backgroundCorrectedAssay) %in%
          names(assays(se))))
  {
    stop(
      "nucleolus ",
      "genome ",
      backgroundCorrectedAssay,
      " should be the assays of se."
    )
  }
  
  normalization.method <-
    match.arg(normalization.method,
              c("none", "scale", "quantile", "cyclicloess"))


  combineP.method <-
    match.arg(combineP.method,
              c(
                "Browns",
                "minimump",
                "logitp",
                "Fishers",
                "sumz",
                "meansig"
              ))
  smooth.method <-
    match.arg(smooth.method, 
        c("loess", "none", "butterworthfilter")) 
  gr <- rowRanges(se)
  windowSize <- median(width(gr))
  ## normalization among ratios
  bc <- as.data.frame(assays(se)[[backgroundCorrectedAssay]])
  bc.norm <- normalizeBetweenArrays(bc, method = normalization.method)
  bc.norm.rowMeans <- rowMeans(bc.norm)
 
  if (smooth.method == "loess") 
  {
      positions <- 1:length(bc.norm.rowMeans)
      span <- N/length(positions)
      loess.fit <- loess(bc.norm.rowMeans ~ positions, span = span, ...)
      bc.smoothed <- predict(loess.fit, positions, se = FALSE)
  }
  else if (smooth.method == "butterworthfilter")
  {
      bc.smoothed <- butterFilter(bc.norm.rowMeans, N = N)
  }
  else
  {
     bc.smoothed = bc.norm.rowMeans
  }
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
  ## assume points between previous valley and next valley belong to the group of this peak
  group <- rep(x, times)
  if (length(group) != length(bc.smoothed))
  {
    stop(
      "The length of group is not identical with that of signals.",
      "Please report this bug.")
  }
  
  fit <- lmFit(bc.norm)
  fit2 <- treat(fit, lfc=lfc, trend = TRUE, ...)
  res <- topTable(fit2, number = nrow(fit2), sort.by = "none")
  res <- cbind(bc.norm, res)
  
  #### group windows by valley plus points after previous valley
  res$group <- group
  
  mcols(gr) <- DataFrame(res)
  gr.all <- gr
  
  keep <- gr$adj.P.Val < cutoffAdjPvalue
  gr <- gr[keep]
  se <- se[keep,]
  
  gr <- gr[rowMeans(cbind(assays(se)[["nucleolus"]],
                          assays(se)[["genome"]])) >= countFilter]
  se <- se[rowMeans(cbind(assays(se)[["nucleolus"]],
                          assays(se)[["genome"]])) >= countFilter] 
  if (length(gr) > 0) {
    
    gr <- split(gr, gr$group)
    
    gr.rd <- endoapply(gr, function(.e) {
      ## find the peak summit for each groupped window
      if (length(.e) > 10)
      {
        sig <- loess.smooth(
          x = seq_along(.e),
          y = .e$AveExpr,
          degree = 2,
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
          .e[.z & .e$AveExpr >= max(.leftmin, .rightmin, na.rm = TRUE)]
      }
      ra <- range(.e)
      if (length(.e) > 0)
      {
        all.windows.in.ra <-
                 gr.all[seqnames(gr.all) == seqnames(ra) &
                 start(gr.all)  >= start(ra) &
                 end(gr.all)  <= end(ra),]
      
        ra$AveSig <- mean(all.windows.in.ra$AveExpr)
        #ra$AveSig <- quantile(.e$AveExpr, probs=c(0, .75, 1), na.rm=TRUE)[2]
        if (combineP.method  == "Browns")
        {
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
        } else if (combineP.method == "Fishers")
        {
          ra$P.value <- sumlog(all.windows.in.ra$P.Value)$p
          ra$adj.P.Val <- sumlog(all.windows.in.ra$adj.P.Val)$p
        } else if (combineP.method == "logitp")
        {
          ra$P.value <- logitp(all.windows.in.ra$P.Value)$p
          ra$adj.P.Val <- logitp(all.windows.in.ra$adj.P.Val)$p
        } else if (combineP.method == "minimump")
        {
          ra$P.value <- minimump(all.windows.in.ra$P.Value)$p
          ra$adj.P.Val <- minimump(all.windows.in.ra$adj.P.Val)$p
        } else if (combineP.method == "sumz")
        {
          print("OK")
          ra$P.value <- sumz(all.windows.in.ra$P.Value)$p
          ra$adj.P.Val <- sumz(all.windows.in.ra$adj.P.Val)$p
        }
        mcols(ra) <- cbind(DataFrame(t(colMeans(as.matrix(
          mcols(all.windows.in.ra)[, 1:dim(bc.norm)[2]])))), mcols(ra))
      }
      return(ra)
    }
    )
    gr <- unlist((gr.rd))
    
    if (combineP.method == "meansig")
    {
      fit <- lmFit(mcols(gr)[, 1:dim(bc.norm)[2]])
      fit2 <- treat(fit, trend = TRUE, lfc = lfc, ...)
      res <- topTable(fit2, number = nrow(fit2), sort.by = "none")
      mcols(gr)$P.Value =  res$P.Value
      mcols(gr)$adj.P.Val = res$adj.P.Val
      mcols(gr)$AveSig  = res$logFC
    }
    gr <- unique(sort(gr))
    gr <- gr[mcols(gr)$adj.P.Val < cutoffAdjPvalue,]
  }
    ## avoid overlaps
    wh <- ceiling(windowSize/2)
    gr <- gr[width(gr) > 2*wh]
    start(gr) <- start(gr) + wh
    end(gr) <- end(gr) - wh
    gr
}
