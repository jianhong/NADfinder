#' Transform counts to log2 cpm ratios, log2 ratios or log2 odds ratios
#'
#' @description calculate the log2 ratios, log2 cpm (count per million) ratios, or 
#' log2 odds ratios for nucleolus vs genome. Pseudo-count will be used to avoid 
#' zero division or log(0).
#'
#' @param A,B window-level counts for nucleolus and genome, extracted from the assays of 
#' the output of the tileCounts function
#' @param seqnames.A,seqnames.B seqnames, extracted from the rowRanges of the ouput of 
#' the tileCounts function
#' @param pseudo.count pseudo-count will be used to aviod x/0 or log0, defult to 1.
#' @param transformation transformation type
#' @param chrom.level.lib indicating whether calculating CPM or odds using 
#' sequence depth of the whole genome or the corresponding chromosome
#' @param lib.size.A,lib.size.B library size for A and B.  
#' these two dataframes contain chromosome-level sequence depth for the chromosomes,
#' which can be extracted from the metadata of the output of the tileCounts function
#' @export
#' @return a numeric vector of log2 ratios, log2 CPM ratios or log2 odds ratios.
#' @examples
#' transformData(seq_len(10), 10:1, seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2OddsRatio",
#' chrom.level.lib = FALSE, lib.size.A = cbind(c("chr1", "chr2"), c(10000, 12000)), 
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 12000)))
#' transformData(seq_len(10), 10:1, seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)), 
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2CPMRatio",
#' chrom.level.lib = FALSE, lib.size.A = cbind(c("chr1", "chr2"), c(10000, 12000)), 
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 12000)))
#' transformData(seq_len(10), 10:1, seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)), 
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2CPMRatio",
#' chrom.level.lib = TRUE, lib.size.A = cbind(c("chr1", "chr2"), c(100, 12000)), 
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 200)))
#' transformData(seq_len(10), 10:1, seqnames.A = Rle(c("chr1", "chr2" ) , c(5,5)),
#' Rle(c("chr1", "chr2" ) , c(5,5)), transformation = "log2OddsRatio",
#' chrom.level.lib = TRUE, lib.size.A = cbind(c("chr1", "chr2"), c(100, 12000)),
#' lib.size.B = cbind(c("chr1", "chr2"), c(10000, 200)))
#' transformData(seq_len(10), 10:1, transformation = "log2Ratio")
#' @author Julie Zhu and Haibo Liu

transformData <- function(A, B, seqnames.A, seqnames.B, pseudo.count = 1L, 
    transformation = c("log2OddsRatio", "log2CPMRatio", "log2Ratio"),
    chrom.level.lib = TRUE, lib.size.A, lib.size.B)
{ 
    stopifnot(length(A) == length(B))
    stopifnot(inherits(A, c("numeric", "integer")))
    stopifnot(inherits(B, c("numeric", "integer")))
    stopifnot(all(A >= 0), all(B >= 0))
    
    transformation <- match.arg(transformation)
    
    if (transformation != "log2Ratio")
    {
        stopifnot(dim(lib.size.A)[2] >= 1 && all(as.numeric(lib.size.A[, 1]) >= 0))
        stopifnot(dim(lib.size.B)[2] >= 1 && all(as.numeric(lib.size.B[, 1]) >= 0))
        stopifnot(dim(lib.size.A)[1] == dim(lib.size.B)[1] &&
            length(seqnames.A) == length(seqnames.B) && length(seqnames.A) == length(A))
    }
    
    ## if any element of A or B is zero, add psudocount to that element of A and B
    ## This will not bias the data
    A[A==0 | B==0] <- A[A==0 | B==0] + pseudo.count
    B[A==0 | B==0] <- B[A==0 | B==0] + pseudo.count
    
    if (transformation != "log2Ratio")
    {
        if(chrom.level.lib)
        {
            A <- data.frame(A, seqnames = as.character(seqnames.A), stringsAsFactors = FALSE)
            A <- merge(A, lib.size.A, by.x= "seqnames", by.y ="row.names", all.x = TRUE)
            colnames(A) <- c("seqnames", "count.A", "lib.size.A")
            
            B <- data.frame(B, seqnames = as.character(seqnames.B), stringsAsFactors = FALSE)
            B <- merge(B, lib.size.B, by.x = "seqnames", by.y = "row.names", all.x = TRUE)
            colnames(B) <- c("seqnames", "count.B", "lib.size.B")
        } else
        {
            genome.library.size.A <-  sum(as.numeric(lib.size.A[,1]))
            genome.library.size.B <-  sum(as.numeric(lib.size.B[,1]))
        }
        
        if (transformation == "log2OddsRatio")
        {
            if (!chrom.level.lib)
            {
                r <- log2(A/(genome.library.size.A - A)) - log2(B/(genome.library.size.B - B))
            } else 
            {

                odds.A <- as.numeric(A[, 2]) / (as.numeric(A[,3]) - as.numeric(A[, 2]))
                odds.B <- as.numeric(B[, 2]) / (as.numeric(B[,3]) - as.numeric(B[, 2]))
                r <- log2(odds.A / odds.B)
            }
        } else if (transformation == "log2CPMRatio")
        {
            if (!chrom.level.lib)
            {        
                r <- log2(A/genome.library.size.A) - log2(B/genome.library.size.B)
            } else 
            {
                cpm.A <- as.numeric(A[,2]) / as.numeric(A[,3])
                cpm.B <- as.numeric(B[,2]) / as.numeric(B[,3])
                r <- log2(cpm.A / cpm.B)
            }
        }
    } else 
    {
        r <- log2(A / B)
    }
    ## should I control the number of digits of r?
    r 
}
