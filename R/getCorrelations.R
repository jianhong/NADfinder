#' get correlations for replicates
#' @description Get the correlations of replicates by the coverage of peaks.
#'              The signals will be filter by the background cutoff value and
#'              the correlations will be calculated.
#' @param se A \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
#' The output of \link{log2se}.
#' @param chr A vector of character. Filter for seqnames. It should be the
#' chromosome names to be kept.
#' @param ratioAssay character(1). 
#' Column name of ratio for correlation calculation.
#' @param window numeric(1) or integer(1). 
#' The window size for summary of the ratios.
#' @param cutoff numeric(1). All the coverages lower than cutoff value 
#' in a given window will be filtered out.
#' @param method A character string indicating which correlation coefficient
#'  is to be computed. See \link[stats]{cor}.
#' @param ... Parameters not used.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom stats as.formula
#' @importFrom utils combn
#' @importFrom IRanges Views viewSums
#' @export
#' @return A list of matrixes of correlation and coefficient.
#' @examples
#' data(triplicates.counts)
#' se <- triplicates.counts
#' gps <- c("26", "28", "29")
#' se <- log2se(se, 
#'              nucleosomeCols = paste0("N", gps, ".bam"),
#'              genomeCols = paste0("G", gps, ".bam"))
#' getCorrelations(se, chr="chr18")
#'

getCorrelations <- function(se, chr = paste0("chr", seq_len(21)),
                             ratioAssay = "ratio", window=10000, cutoff=1,
                             method=c("spearman", "pearson", "kendall"),
                             ...){
    stopifnot(class(se)=="RangedSummarizedExperiment")
    stopifnot(length(ratioAssay)==1)
    stopifnot(ratioAssay %in% names(assays(se)))
    method <- match.arg(method)
    gr <- rowRanges(se)
    se <- subset(se, seqnames(gr) %in% chr)
    seqlevels(gr) <- chr
    seqlen <- seqlengths(gr)
    seqlen[is.na(seqlen)] <- sapply(names(seqlen)[is.na(seqlen)],
                                    function(.ele){
                                        end(range(gr[seqnames(gr) %in% .ele]))
                                    })
    seqlen <- unlist(seqlen)
    tiles <- tileGenome(seqlen, tilewidth = window)
    tiles <- unlist(tiles)
    tiles <- split(tiles, seqnames(tiles))
    tiles <- lapply(tiles, ranges)
    ## resample the signals
    cn <- colnames(assays(se)[[ratioAssay]])
    if(length(cn)==0){
        cn <- paste0("sample", seq_len(ncol(assays(se)[[ratioAssay]])))
        colnames(assays(se)[[ratioAssay]]) <- cn
    }
    resample <- lapply(cn, function(.ele){
        .ele <- exportSignals(dat = se, 
                              assayName = ratioAssay, 
                              colName = .ele)
        .views <- Views(.ele, tiles[names(.ele)])
        viewSums(.views)
    })
    names(resample) <- cn
    ## swap resamples
    chr <- sort(unique(sapply(resample, names)))
    resample <- lapply(chr, function(.ele){
        do.call(cbind, lapply(resample, function(.e) .e[[.ele]]))
    })
    names(resample) <- chr
    ## filter
    resample <- lapply(resample, function(.ele){
        .ele[rowSums(.ele >= cutoff) == ncol(.ele), , drop=FALSE]
    })
    ## rbind
    resample <- do.call(rbind, resample)
    ## correlation
    cor <- cor(resample, method = method)
    ## coefficient
    coe <- function(df){
        coln <- colnames(df)
        out <- matrix(0, nrow=length(coln), ncol=length(coln),
                      dimnames = list(coln, coln))
        cb <- combn(coln, m=2, simplify = FALSE)
        r.squared <- lapply(cb, function(.ele){
            lm <- lm(as.formula(paste(.ele, collapse="~")), data=df)
            summary(lm)$r.squared
        })
        for(i in seq_along(cb)){
            out[cb[[i]][1], cb[[i]][2]] <-
                out[cb[[i]][2], cb[[i]][1]] <-
                r.squared[[i]]
        }
        for(i in seq_along(coln)){
            out[i, i] <- 1
        }
        out
    }
    coe <- coe(as.data.frame(resample))
    return(list(cor=cor, coe=coe))
}
