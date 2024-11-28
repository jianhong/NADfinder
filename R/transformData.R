#' transform counts to log2 cpm ratios, log2 ratios or log2 odds ratios
#'
#' @description calculate the log2 ratios, log2 cpm (count per million) ratios,
#' or log2 odds ratios for nucleolus vs genome. pseudo-count will be used to
#' avoid x/0 or log(0).
#'
#' @param A,B window-level counts for nucleolus and genome, extracted from
#' the assays of the output of the tileCounts function
#' @param seqnames.A,seqnames.B seqnames, extracted from the rowRanges of
#' the ouput of the tileCounts function
#' @param pseudo.count pseudo-count will be used to aviod x/0 or log0, defult
#' to 1.
#' @param transformation transformation type
#' @param chrom.level.lib indicating whether calculating CPM or odds using
#' sequence depth of the whole genome or the corresponding chromosome
#' @param lib.size.A,lib.size.B library size for A and B.
#' these two dataframes contain chromosome-level sequence depth for the
#' chromosomes, which can be extracted from the metadata of the output of
#' the tileCounts function
#' @export
#' @return a numeric vector of log2 ratios, log2 CPM ratios or log2 odds ratios.
#' @examples
#' transformData(seq_len(10), 10:1,
#'              seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2OddsRatio",
#' chrom.level.lib = FALSE,
#' lib.size.A = cbind(c("chr1", "chr2"), c(10000, 12000)),
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 12000)))
#' transformData(seq_len(10), 10:1,
#' seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2CPMRatio",
#' chrom.level.lib = FALSE,
#' lib.size.A = cbind(c("chr1", "chr2"), c(10000, 12000)),
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 12000)))
#' transformData(seq_len(10), 10:1,
#' seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2CPMRatio",
#' chrom.level.lib = TRUE,
#' lib.size.A = cbind(c("chr1", "chr2"), c(100, 12000)),
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 200)))
#' transformData(seq_len(10), 10:1,
#' seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2OddsRatio",
#' chrom.level.lib = TRUE,
#' lib.size.A = cbind(c("chr1", "chr2"), c(100, 12000)),
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 200)))
#' transformData(seq_len(10), 10:1, transformation = "log2Ratio")
#' @author Julie Zhu

transformData <-
    function(A, B,
             seqnames.A, seqnames.B,
             pseudo.count = 1L,
             transformation = c("log2OddsRatio",
                                "log2CPMRatio",
                                "log2Ratio"),
             chrom.level.lib = TRUE,
             lib.size.A, lib.size.B) {
    stopifnot(length(A) == length(B))
    transformation <- match.arg(transformation)
    stopifnot(inherits(A, c("numeric", "integer")))
    stopifnot(inherits(B, c("numeric", "integer")))
    stopifnot(all(A>=0))
    stopifnot(all(B>=0))
    if (transformation != "log2Ratio")
    {
        stopifnot(dim(lib.size.A)[2] >= 2 &&
                      all(as.numeric(lib.size.A[,2]) >=0))
        stopifnot(dim(lib.size.B)[2] >= 2 &&
                      all(as.numeric(lib.size.B[,2]) >=0))
        stopifnot(dim(lib.size.A)[1] == dim(lib.size.B)[1] &&
            length(seqnames.A) == length(seqnames.B) &&
                length(seqnames.A) == length(A))
    }
    if (transformation != "log2Ratio" && chrom.level.lib)
    {
        A <- cbind(A, as.character(seqnames.A),A)
        A[,3] <- lib.size.A[match(A[,2], lib.size.A[,1]), 2]
        B <- cbind(B, as.character(seqnames.B), B)
        B[,3] <- lib.size.B[match(B[,2], lib.size.B[,1]), 2]
    }
    if (transformation != "log2Ratio" && !chrom.level.lib)
    {
        genome.library.size.A <-  sum(as.numeric(lib.size.A[,2]))
        genome.library.size.B <-  sum(as.numeric(lib.size.B[,2]))
    }
    if (transformation == "log2OddsRatio")
    {
        if (!chrom.level.lib)
        {
             r <- log2((A + pseudo.count)/(genome.library.size.A - A)) -
                log2((B + pseudo.count)/(genome.library.size.B - B))
        } else {
             odds.A <-
                   (as.numeric(A[,1]) + pseudo.count ) /
                        (as.numeric(A[,3]) - as.numeric(A[,1]))

             odds.B <-
                   (as.numeric(B[,1]) + pseudo.count ) /
                        (as.numeric(B[,3]) - as.numeric(B[,1]))
	     r <- log2(odds.A / odds.B)
        }
    } else if (transformation == "log2CPMRatio") {
        if (!chrom.level.lib)
        {
            r <- log2((A + pseudo.count)/(genome.library.size.A +
                                              pseudo.count)) -
                log2((B + pseudo.count)/(genome.library.size.B +
                                             pseudo.count))
        } else {
             cpm.A <- (as.numeric(A[,1]) + pseudo.count ) /
                 (as.numeric(A[,3]) + pseudo.count) * 1e+6
             cpm.B <- (as.numeric(B[,1]) + pseudo.count ) /
                 (as.numeric(B[,3]) + pseudo.count) * 1e+6
             r <- log2(cpm.A / cpm.B)
        }
   } else if (transformation == "log2Ratio") {
        r <- log2(A + pseudo.count) - log2(B + pseudo.count)
   }
   r
}
