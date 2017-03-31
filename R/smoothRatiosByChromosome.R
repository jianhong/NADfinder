#' smooth the ratios by chromosome
#'
#' Split the ratios by chromosome and do background correction and smooth.
#'
#' @param se An object of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} 
#' with scores. Output of \link{log2se}
#' @param chr A vector of character. Filter for seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay The name of assay in se, which store the values 
#' to be smoothed. 
#' @param backgroundCorrectionAssay,smoothedRatioAssay,zscoreAssay character(1).
#' Assays names for background correction ratios, smoothed ratios and 
#' z-score based on background correction ratios.
#' @param backgroundPercentage numeric(1). Percentage of values for background, 
#' see \link{zscoreOverBck}. How many percent lower values will be treated as
#' background.
#' @param ... Parameters could be passed to \link{butterFilter}.
#'
#' @export
#' @return A \link[S4Vectors]{SimpleList} of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}
#'  with smoothed ratios.
#' @examples
#'
#' data(single.count)
#' se <- single.count
#' dat <- log2se(se, nucleosomeCols="nucleosome.bam", genomeCols="genome.bam")
#' dat <- smoothRatiosByChromosome(dat, N=100)
#'
smoothRatiosByChromosome <- function(se,
                                     chr=paste0("chr", 
                                                c(seq_len(21), "X", "Y")),
                                     ratioAssay="ratio",
                                     backgroundCorrectionAssay="bcRatio",
                                     smoothedRatioAssay="smoothedRatio",
                                     zscoreAssay="zscore", 
                                     backgroundPercentage=0.25,
                                     ...){
    stopifnot(class(se)=="RangedSummarizedExperiment")
    stopifnot(length(ratioAssay)==1)
    stopifnot(ratioAssay %in% names(assays(se)))
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
        assays(.ele)[[zscoreAssay]] <- 
            apply(assays(.ele)[[smoothedRatioAssay]],
                  2, function(.e) zscoreOverBck(.e, backgroundPercentage))
        .ele
    })
    SimpleList(se) #do.call(rbind, se)
}
