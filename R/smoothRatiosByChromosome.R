#' Backgound correction and signal smoothing per chromosome
#'
#' Split the ratios by chromosome and do background correction and signal
#' smoothing. The future_lapply() function is used for multithreading.
#'
#' @param se An [SummarizedExperiment::RangedSummarizedExperiment-class] object.
#' with log2-transformed ratios, CPMRatios or OddRatios. Output of \link{log2se}.
#' @param chr A character vector, used to filter out seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay The name of assay in se, which store the values
#' (log2-transformed ratios,
#' CPMRatios or OddRatios) to be smoothed.
#' @param backgroundCorrectedAssay,smoothedRatioAssay,zscoreAssay character(1).
#' Assays names for background corrected ratios, smoothed ratios and
#' z-scores based on background corrected ratios.
#' @param backgroundPercentage numeric(1). Percentage of values for background,
#' see \link{zscoreOverBck}. The percentage of values lower than this threshold
#' will be treated as background, with 25 percentile as default.
#' @param chrom.level.background logical(1): TRUE or FALSE, default to TRUE,
#' use chromosome-level background to calculate z-score.
#' @param future.scheduling An integer(1), average number of futures ("chunks")
#' per worker for multithreading using the [future.apply::future_lapply()]
#' function. If 0.0, then a single future is used to process all elements of
#' X. If 1.0 or TRUE, then one future per worker is used. If 2.0, then each
#' worker will process two futures (if there are enough elements in X). If
#' Inf or FALSE, then one future per element of X is used.
#' @param ... Parameters could be passed to \link{butterFilter}.
#' @importFrom methods is
#' @importFrom future.apply future_lapply
#' @export
#' @return A \link[S4Vectors:SimpleList-class]{SimpleList} of
#' [SummarizedExperiment::RangedSummarizedExperiment-class] with smoothed
#' ratios.
#' @author Jianhong Ou, Haibo Liu and Julie Zhu
#' @examples
#' library(future.apply)
#' if (.Platform$OS.type == "windows")
#' {
#'     plan(sequential)
#' } else {
#'     plan(multisession)
#' }
#' data(single.count)
#' se <- single.count
#' dat <- log2se(se, nucleolusCols="N18.subsampled.srt.bam",
#'               genomeCols="G18.subsampled.srt.bam",
#'               transformation="log2CPMRatio")
#' dat1 <- smoothRatiosByChromosome(dat, N=100, chr = c("chr18", "chr19"))
#' dat2 <- smoothRatiosByChromosome(dat, N=100, chr = c("chr18", "chr19"),
#'                                  chrom.level.background = FALSE)
#'

smoothRatiosByChromosome <- function(se,
                                     chr = paste0("chr",
                                                  c(seq_len(21), "X", "Y")),
                                     ratioAssay = "ratio",
                                     backgroundCorrectedAssay = "bcRatio",
                                     smoothedRatioAssay = "smoothedRatio",
                                     zscoreAssay = "zscore",
                                     backgroundPercentage = 0.25,
                                     chrom.level.background = TRUE,
                                     future.scheduling = 1L,
                                     ...)
{
    stopifnot(is(se, "RangedSummarizedExperiment"))
    stopifnot(length(ratioAssay) == 1)
    stopifnot(ratioAssay %in% names(assays(se)))

    ## for computing z-score using global background
    if (!chrom.level.background)
    {
        all.ratios <- assays(se)$ratio
        bg.percentile <-
            quantile(all.ratios, probs = backgroundPercentage)
        bg.ratios <- all.ratios[all.ratios <= bg.percentile]

        pop.sd <-
            sd(bg.ratios) * sqrt(length(bg.ratios) - 1) /
            sqrt(length(bg.ratios))
        pop.mean <- mean(bg.ratios)
    }

    ## split se to get rse for each chromosome
    se <- split(se, as.character(seqnames(rowRanges(se))))
    se <- se[names(se) %in% chr]

    se <- future_lapply(se, function(.ele)
    {
        if (!is(assays(.ele)[[ratioAssay]], "matrix"))
        {
            assays(.ele)[[ratioAssay]] <- as.matrix(assays(.ele)[[ratioAssay]])
        }

        bkgcorrect <- apply(assays(.ele)[[ratioAssay]],
                            2, backgroundCorrection)
        rownames(bkgcorrect) <- rownames(assays(.ele)[[ratioAssay]])
        assays(.ele)[[backgroundCorrectedAssay]] <- bkgcorrect

        smooth <- apply(assays(.ele)[[backgroundCorrectedAssay]],
                        2, function(.e) butterFilter(.e))
        rownames(smooth) <- rownames(assays(.ele)[[ratioAssay]])
        assays(.ele)[[smoothedRatioAssay]] <- smooth

        if (chrom.level.background)
        {
            zscore <- apply(assays(.ele)[[smoothedRatioAssay]],
                      2, function(.e)
                          zscoreOverBck(.e, backgroundPercentage))
            rownames(zscore) <- rownames(assays(.ele)[[ratioAssay]])
            assays(.ele)[[zscoreAssay]] <- zscore

        } else {
            zscore <- apply(assays(.ele)[[smoothedRatioAssay]],
                      2, function(.e) {(.e - pop.mean) / pop.sd})
            rownames(zscore) <- rownames(assays(.ele)[[ratioAssay]])
            assays(.ele)[[zscoreAssay]] <- zscore
        }
        .ele
    }, future.scheduling = future.scheduling)
    se <- SimpleList(se)
    se
}
