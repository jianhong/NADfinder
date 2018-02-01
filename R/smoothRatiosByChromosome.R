#' backgound correction and smooth by chromosome
#'
#' Split the ratios by chromosome and do background correction and smooth.
#'
#' @param se An object of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} 
#' with scores. Output of \link{log2se}
#' @param chr A character vector, used to filter out seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay The name of assay in se, which store the values 
#' to be smoothed. 
#' @param backgroundCorrectionAssay,smoothedRatioAssay,zscoreAssay character(1).
#' Assays names for background corrected ratios, smoothed ratios and 
#' z-scores based on background corrected ratios.
#' @param backgroundPercentage numeric(1). Percentage of values for background, 
#' see \link{zscoreOverBck}. How many percent lower values will be treated as
#' background, default lower 25 percent.
#' chrom.level.background TRUE or FALSE, default to TRUE, use chromosome-level background
#' to calculate z-score
#' @param chrom.level.background logical(1).
#' @param ... Parameters could be passed to \link{butterFilter}.
#' @importFrom methods is
#' @export
#' @return A \link[S4Vectors]{SimpleList} of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}
#' with smoothed ratios.
#' @author Jianhong Ou and Julie Zhu
#' @examples
#'
#' data(single.count)
#' se <- single.count
#' dat <- log2se(se, nucleosomeCols="nucleosome.bam", genomeCols="genome.bam", 
#' transformation="log2Ratio")
#' dat1 <- smoothRatiosByChromosome(dat, N=100)
#' dat2 <- smoothRatiosByChromosome(dat, N=100, chrom.level.background = FALSE)
#'
smoothRatiosByChromosome <- function(se,
                                     chr=paste0("chr", 
                                                c(seq_len(21), "X", "Y")),
                                     ratioAssay="ratio",
                                     backgroundCorrectionAssay="bcRatio",
                                     smoothedRatioAssay="smoothedRatio",
                                     zscoreAssay="zscore", 
                                     backgroundPercentage=0.25, chrom.level.background = TRUE,
                                     ...){
    stopifnot(is(se, "RangedSummarizedExperiment"))
    stopifnot(length(ratioAssay)==1)
    stopifnot(ratioAssay %in% names(assays(se)))
    if (!chrom.level.background)
    {
        all.ratios <- assays(se)$ratio
        bg.percentile <- quantile(all.ratios, probs = backgroundPercentage)
        bg.ratios <- all.ratios[all.ratios <= bg.percentile]
        pop.sd <- sd(bg.ratios) * sqrt(length(bg.ratios) - 1)/sqrt(length(bg.ratios))
        pop.mean <- mean(bg.ratios) 
    }
    se <- split(se, as.character(seqnames(rowRanges(se))))
    se <- se[names(se) %in% chr]
    se <- lapply(se, function(.ele){
        if(class(assays(.ele)[[ratioAssay]])!="matrix"){
            assays(.ele)[[ratioAssay]] <- as.matrix(assays(.ele)[[ratioAssay]])
        }
        assays(.ele)[[backgroundCorrectionAssay]] <- 
            apply(assays(.ele)[[ratioAssay]],
                  2, backgroundCorrection) 
        assays(.ele)[[smoothedRatioAssay]] <- 
            apply(assays(.ele)[[backgroundCorrectionAssay]],
                  2, function(.e) butterFilter(.e, ...))
        if (chrom.level.background)
        {
            assays(.ele)[[zscoreAssay]] <- 
                apply(assays(.ele)[[smoothedRatioAssay]],
                  2, function(.e) zscoreOverBck(.e, backgroundPercentage))
	   
        }
       else
       {
          assays(.ele)[[zscoreAssay]] <-
               apply(assays(.ele)[[smoothedRatioAssay]],
                 2, function(.e) { 
                 (.e - pop.mean)/pop.sd })
       }
        .ele
    })
    SimpleList(se) #do.call(rbind, se)
se
}
