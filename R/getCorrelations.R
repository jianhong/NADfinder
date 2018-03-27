#' Get correlation coefficinets and p-values between biological replicates
#' @description Get correlations and p-values between biological replicates based on coverage signal for 
#'              peak regions. The signals will be filtered by the background cutoff value before
#'              calculated correlations. This function also output a correlation plots using the 
#'              \link[corrplot]{corrplot}.
#' @param se A \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
#' The output from \link{log2se}.
#' @param chr A vector of character. Filter for seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay character(1). 
#' Column name of ratio for correlation calculation.
#' @param window numeric(1) or integer(1). 
#' The window size for summary of the ratios.
#' @param cutoff numeric(1). All the coverage signals lower than cutoff value 
#' in a given window will be filtered out.
#' @param method character(1) indicating which correlation coefficient
#'  is to be computed. See \link[stats]{cor}.
#' @param file_name A file name for output correlation plots
#' @param ... Parameters not used.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom stats as.formula
#' @importFrom utils combn
#' @importFrom corrplot corrplot
#' @importFrom IRanges Views viewSums
#' @export
#' @return A list of matrixes of correlation coefficients and p-values.
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' data(triplicate.count)
#' se <- triplicate.count
#' se <- log2se(se, transformation = "log2CPMRatio",
#'              nucleolusCols = c("N18.subsampled.srt-2.bam",
#'              "N18.subsampled.srt-3.bam",
#'              "N18.subsampled.srt.bam"),
#'              genomeCols = c("G18.subsampled.srt-2.bam",
#'              "G18.subsampled.srt-3.bam",
#'              "G18.subsampled.srt.bam"))
#' getCorrelations(se, chr="chr18")
#'

getCorrelations <- function(se,
                            chr = paste0("chr", seq_len(19)),
                            ratioAssay = "ratio",
                            window = 10000L,
                            cutoff = 1,
                            method = c("spearman", "pearson", "kendall"),
                            file_name = "Correlation plots.pdf",
                            ...) 
{
    stopifnot(class(se) == "RangedSummarizedExperiment")
    stopifnot(length(ratioAssay) == 1)
    stopifnot(ratioAssay %in% names(assays(se)))
    method <- match.arg(method)
    gr <- rowRanges(se)
    se <- subset(se, seqnames(gr) %in% chr)
    gr <- keepSeqlevels(gr, chr, pruning.mode="coarse")
    seqlen <- seqlengths(gr)
    seqlen[is.na(seqlen)] <- sapply(names(seqlen)[is.na(seqlen)],
                                    function(.ele) {
                                        end(range(gr[seqnames(gr) %in% .ele]))})
    seqlen <- unlist(seqlen)
    tiles <- tileGenome(seqlen, tilewidth = window)
    tiles <- unlist(tiles)
    tiles <- split(tiles, seqnames(tiles))
    tiles <- lapply(tiles, ranges)
    ## resample the signals
    cn <- colnames(assays(se)[[ratioAssay]])
    if (length(cn) == 0) 
    {
        cn <- paste0("sample", seq_len(ncol(assays(se)[[ratioAssay]])))
        colnames(assays(se)[[ratioAssay]]) <- cn
    }
    resample <- lapply(cn, function(.ele) {
        .ele <- exportSignals(dat = se,
                              assayName = ratioAssay,
                              colName = .ele)
        .views <- Views(.ele, tiles[names(.ele)])
        viewSums(.views)})
    names(resample) <- cn
    ## swap resamples
    chr <- sort(unique(sapply(resample, names)))
    resample <- lapply(chr, function(.ele) {
        do.call(cbind, lapply(resample, function(.e)
            .e[[.ele]]))})
    names(resample) <- chr
    ## filter
    resample <- lapply(resample, function(.ele) {
        .ele[rowSums(.ele >= cutoff) == ncol(.ele), , drop = FALSE]})
    ## rbind
    resample <- do.call(rbind, resample)
    ## R-squared from linear model fitting
    ## Correlation coefficients can be calculaed using 3 different methods,
    ## But here normality is assumed. So I changed the code to call cor.test
    ## function to get both the correlation coefficients and the p-values.
    
    correlations <- function(df) {
        coln <- colnames(df)
        mat <- matrix(
            0,
            nrow = length(coln),
            ncol = length(coln),
            dimnames = list(coln, coln)
        )
        cb <- combn(coln, m = 2, simplify = FALSE)
        corNp <- lapply(cb, function(.ele) {
            corTestOut <-
                cor.test(as.formula(paste("~", .ele[1], "+", .ele[2])), data = df)
            c(rho = corTestOut$estimate, p = corTestOut$p.value)
        })
        
        fillMatrix <- function(cb, mat, j) {
            for (i in seq_along(cb)) {
                mat[cb[[i]][1], cb[[i]][2]] <-
                    mat[cb[[i]][2], cb[[i]][1]] <-
                    corNp[[i]][j]
            }
            mat
        }
        out <- lapply(1:2, fillMatrix, cb = cb, mat = mat)
        names(out) <- c("cor.coeff", "p.values")
        out
    }
    
    corOut <- correlations(as.data.frame(resample))
    
    ## plot correlation and color the squares if p <= 0.01
    
    col <-
        colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    
    ## correlation plots
    pdf(file_name)
    corrplot(
        corOut$cor.coeff,
        method = "color",
        col = col(200),
        type = "upper",
        order = "hclust",
        addgrid.col = "gray",
        addCoef.col = "black",
        # Add coefficient of correlation
        tl.col = "black",
        tl.srt = 45,
        #Text label color and rotation
        # Combine with significance
        p.mat = corOut$p.values,
        sig.level = 0.01,
        insig = "blank",
        # hide correlation coefficient on the principal diagonal
        diag = FALSE
    )
    dev.off()

    return(invisible(corOut))
}
