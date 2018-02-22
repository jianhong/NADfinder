#' Trim peaks
#'
#' Filter the peaks by pvalue and trim the range of peaks 
#' for an NAD experiment without biological replicates.
#'
#' @param se An object of 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment}
#' with assays of raw counts, ratios, background corrected ratios,
#' smoothed ratios and z-scores. It should be an element of the output of 
#' \link{smoothRatiosByChromosome}
#' @param cutoffPvalue numeric(1). Cutoff p-value.
#' @param backgroundPercentage numeric(1). Cutoff value for the peaks height.
#' @param countFilter numeric(1) or integer(1). Cutoff value for mean of 
#' raw reads count in each window.
#' @param ratioAssay character(1). The name of assay in se, which store the 
#' values to be smoothed. 
#' @param backgroundCorrectionAssay,smoothedRatioAssay,zscoreAssay Assays names 
#' for background-corrected ratios, smoothed ratios and z-scores based on
#' background corrected ratios.
#'
#' @return An object of \link[GenomicRanges]{GRanges}.
#'
#' @export
#' @importFrom stats quantile
#' @examples
#'
#' data(single.count)
#' se <- single.count
#' ## Calculate ratios for peak calling. We use signal vs input.
#' dat <- log2se(se, nucleoleusCols="nucleoleus.bam", genomeCols="genome.bam",
#' transformation = "log2Ratio")
#' ## Smooth the ratios for each chromosome.
#' dat <- smoothRatiosByChromosome(dat, N=100)
#' peaks <- trimPeaks(dat[["chr18"]],
#'                 backgroundPercentage=.25,
#'                 cutoffPvalue=0.05, countFilter=1000)
#'
trimPeaks <- function(se, cutoffPvalue=0.05,
                      backgroundPercentage=.25,
                      countFilter=1000,
                      ratioAssay="ratio",
                      backgroundCorrectionAssay="bcRatio",
                      smoothedRatioAssay="smoothedRatio",
                      zscoreAssay="zscore"){
    backgroundPercentage2=.5
    stopifnot(backgroundPercentage>0 && backgroundPercentage<1)
    stopifnot(class(se)=="RangedSummarizedExperiment")
    asyNames <- c("nucleoleus", "genome", ratioAssay, 
                   backgroundCorrectionAssay,
                   smoothedRatioAssay, zscoreAssay)
    if(any(!asyNames %in% names(assays(se)))){
        stop("se must be a list contain assays of nucleoleus, genome,", 
             paste(ratioAssay, backgroundCorrectionAssay, 
                   smoothedRatioAssay, zscoreAssay,
             sep=", "))
    }
    chr <- rowRanges(se)
    if(!all(width(chr)==width(chr)[1])){
        stop("The width of rowRanges should be identical.")
    }
    windowSize <- width(chr)[1]
    sampleLen <- ncol(assays(se)[[zscoreAssay]])
    
    ## this function is for single sample. For multiple samples, use callPeaks
    if(sampleLen!=1){
        stop("This function is for single sample.", 
             "For multiple samples, please try callPeaks.")
    }
    ## prepare data
    mcols(chr) <- 
        do.call(cbind, assays(se)[asyNames])
    colnames(mcols(chr)) <- asyNames
    ## call peaks by zscore. 
    mcols(chr) <-
        cbind(mcols(chr),
              DataFrame(groupZscores(mcols(chr)[, zscoreAssay])[, -1]))
    r <- quantile(mcols(chr)[, smoothedRatioAssay],
                  probs = c(0, backgroundPercentage, 1))
    ## filter
    chr <-
        chr[chr$pvalue < cutoffPvalue & mcols(chr)[, smoothedRatioAssay]>r[2]]
    chr <-
        chr[rowMeans(as.data.frame(
            mcols(chr)[, seq_len(which(colnames(mcols(chr)) == 
                                           ratioAssay) - 1)])) > 
                countFilter]
    ## return if not peak
    if(length(chr)==0){
        mcols(chr) <- DataFrame(zscore=numeric(0), pvalue=numeric(0))
        colnames(mcols(chr)) <- c(zscoreAssay, "pvalue")
        return(chr)
    }
    chr <- split(chr, chr$group)
    ## trim peaks by cut from shoulder
    chr.rd <- lapply(chr, function(.e) {
        .idx <- which(.e$grp.zscore[1]==mcols(.e)[, zscoreAssay])[1]
        if(!is.na(.idx)){
            .leftmin <- min(mcols(.e)[seq_len(.idx), smoothedRatioAssay],
                            na.rm=TRUE)
            .rightmin <- min(mcols(.e)[.idx:length(.e), smoothedRatioAssay],
                             na.rm=TRUE)
            .z <- zscoreOverBck(mcols(.e)[, backgroundCorrectionAssay],
                                backgroundPercentage2)
            .e <- .e[.z>0 & mcols(.e)[, smoothedRatioAssay] >
                         max(.leftmin, .rightmin, na.rm=TRUE)]
        }
        .e <- .e[-length(.e)]
        ra <- range(.e)
        if(length(.e)>0){
            mcols(ra)[, zscoreAssay] <- .e$grp.zscore[1]
            ra$pvalue <- .e$pvalue[1]
        }else{
            mcols(ra) <- DataFrame(zscore=numeric(0), pvalue=numeric(0))
            colnames(mcols(ra)) <- c(zscoreAssay, "pvalue")
        }
        ra
    })
    chr.rd <- unlist(GRangesList(chr.rd))
    ## trim the peaks to avoid overlaps
    ## wh <- ceiling(windowSize/2)
    wh <- ceiling(windowSize/20)
    start(chr.rd) <- start(chr.rd) + wh
    end(chr.rd) <- end(chr.rd) - wh
    chr.rd <- sort(chr.rd)
    s <- start(chr.rd)
    e <- end(chr.rd)
    idx <- which(e[-length(e)]>s[-1])
    if(length(idx)>0){
        end(chr.rd)[idx] <- start(chr.rd)[idx+1] - wh
        start(chr.rd)[idx+1] <- start(chr.rd)[idx+1] + wh
    }
    chr.rd
}
