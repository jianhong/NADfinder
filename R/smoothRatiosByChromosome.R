#' Backgound correction and signal smoothing per chromosome
#'
#' Split the ratios by chromosome and do background correction and signal smoothing.
#'
#' @param se An object of 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} 
#' with log2-transformed ratios, CPMRatios or OddRatios. Output of \link{log2se}
## Why not include these parameters in tileCount.R??
#' @param chr A character vector, used to filter out seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay The name of assay in se, which store the values (log2-transformed ratios,
#' CPMRatios or OddRatios) to be smoothed. 
#' @param backgroundCorrectedAssay,smoothedRatioAssay,zscoreAssay character(1).
#' Assays names for background corrected ratios, smoothed ratios and 
#' z-scores based on background corrected ratios.
#' @param backgroundPercentage numeric(1). Percentage of values for background, 
#' see \link{zscoreOverBck}. The percentage of values lower than this threshold 
#' will be treated as background, with 25 percentile as default.
#' @param chrom.level.background logical(1): TRUE or FALSE, default to TRUE, use chromosome-level 
#' background to calculate z-score
#' @param ... Parameters could be passed to \link{butterFilter}.
#' @importFrom methods is
#' @export
#' @return A \link[S4Vectors:SimpleList-class]{SimpleList} of 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' with smoothed ratios.
#' @author Jianhong Ou, Haibo Liu and Julie Zhu
#' @examples
#'
#' data(single.count)
#' se <- single.count
#' dat <- log2se(se, nucleolusCols="N18.subsampled.srt.bam", genomeCols="G18.subsampled.srt.bam", 
#' transformation="log2CPMRatio")
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
            sd(bg.ratios) * sqrt(length(bg.ratios) - 1) / sqrt(length(bg.ratios))
        pop.mean <- mean(bg.ratios)
    }
    
    ## split se to get rse for each chromosome
    se <- split(se, as.character(seqnames(rowRanges(se))))
    se <- se[names(se) %in% chr]
    
    se <- lapply(se, function(.ele) 
    {
        if (class(assays(.ele)[[ratioAssay]]) != "matrix") 
        {
            assays(.ele)[[ratioAssay]] <- as.matrix(assays(.ele)[[ratioAssay]])
        }
        assays(.ele)[[backgroundCorrectedAssay]] <-
            apply(assays(.ele)[[ratioAssay]],
                  2, backgroundCorrection)
        
        assays(.ele)[[smoothedRatioAssay]] <-
            apply(assays(.ele)[[backgroundCorrectedAssay]],
                  2, function(.e)
                      butterFilter(.e, ...))
        
        if (chrom.level.background)
        {
            assays(.ele)[[zscoreAssay]] <-
                apply(assays(.ele)[[smoothedRatioAssay]],
                      2, function(.e)
                          zscoreOverBck(.e, backgroundPercentage))
            
        } else
        {
            assays(.ele)[[zscoreAssay]] <-
                apply(assays(.ele)[[smoothedRatioAssay]],
                      2, function(.e) {
                          (.e - pop.mean) / pop.sd })
        }
        .ele
    })
    se <- SimpleList(se)
    se
}
