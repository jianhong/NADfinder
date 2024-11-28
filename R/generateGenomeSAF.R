#' Generate sliding windows of a specified width and frequency
#'
#' @param bamfile A character(1), specifying a path to a BAM file, from which
#'   chromosome names and lengths are extracted for generating sliding windows.
#' @param windowSize An integer(1), specifying the width of window in bp.
#' @param step An integer(1), specifying the step in bp for generating sliding
#'   windows.
#' @param chrExcluded A character(n), specifying chromosomes to be excluded
#'   from generating sliding windows.
#'
#' @return A data frame containing sliding windows for each chromosome and
#'   individual chromosome in the SAF (simplified annotation format) format
#'   as described in the Rsubread package.
#' @export
#' @importFrom Rsamtools seqinfo BamFile
#' @importFrom GenomicRanges makeGRangesFromDataFrame slidingWindows
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="NADfinder"),
#'                       recursive=FALSE,
#'                       pattern="*bam$",
#'                       full=TRUE)
#' genomeSW_SAF <-
#'     generateGenomeSAF(bamfile = bamfiles[1],
#'                       windowSize = 50000L,
#'                       step = 5000L,
#'                       chrExcluded = c("chrM", "MT"))

generateGenomeSAF <-
    function(bamfile = NULL,
             windowSize = 50000L,
             step = 5000L,
             chrExcluded = c("chrM", "MT"))
    {
        stopifnot(file.exists(bamfile))
        stopifnot(windowSize %% 1 == 0,
                  step %% 1 ==0,
                  windowSize > 0,
                  step > 0,
                  step < windowSize)
        stopifnot(is.character(chrExcluded))

        genome_seqinfo <- seqinfo(BamFile(bamfile))
        genome_info <- as.data.frame(genome_seqinfo)
        genome_info <- genome_info[!rownames(genome_info) %in% chrExcluded, ]
        genome_df <- data.frame(seqnames = rownames(genome_info),
                                  start = 1,
                                  end = genome_info$seqlengths,
                                  strand = "*")
        genome_GR <- makeGRangesFromDataFrame(genome_df,
                                              seqinfo = genome_seqinfo,
                                              seqnames.field="seqnames",
                                              start.field="start",
                                              end.field= "end",
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)

        slidingWindows <- slidingWindows(genome_GR,
                                            width = windowSize,
                                            step= step)
        slidingWindows_saf <- lapply(slidingWindows, granges2saf)
        slidingWindows_saf <- do.call("rbind", slidingWindows_saf)
        slidingWindows_saf$GeneID <- seq_len(nrow(slidingWindows_saf))
        chr_saf <- granges2saf(genome_GR)
        chr_saf$GeneID<- paste0("genome_",chr_saf$GeneID)
        saf <- rbind(slidingWindows_saf, chr_saf)
        saf
    }

#' Convert genomic features from GRanges to Simplified Annotation Format (SAF)
#'
#' @param granges An object of [GenomicRanges::GRanges-class]
#' @export
#' @return A data frame with columns: "GeneID", "Chr", "Start", "End","Strand".
#'
#' @examples
#' library("GenomicRanges")
#' gr0 <- GRanges(Rle(
#'     c("chr2", "chr2", "chr1", "chr3"),
#'     c(1, 3, 2, 4)
#' ), IRanges(seq_len(10), width = 10:1))
#' saf <- granges2saf(gr0)
#'
granges2saf <- function(granges) {
    if (!is(granges, "GRanges")) {
        stop(deparse(substitute(granges)), " is not a GRanges object.")
    }
    saf <- as.data.frame(granges)[, c(seq_len(3), 5)]
    saf$GeneID <- rownames(saf)
    rownames(saf) <- NULL
    colnames(saf) <- c("Chr", "Start", "End", "Strand", "GeneID")
    saf_colnames <- c("GeneID", "Chr", "Start", "End", "Strand")
    saf <- saf[, saf_colnames]
    saf
}
