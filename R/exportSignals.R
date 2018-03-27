#' Output signals for visualization
#'
#' Output signals to bedgraph, bed, wig, etc, for track viewer
#'
#' @param dat An object of \link[GenomicRanges:GRanges-class]{GRanges}, 
#' or \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' with assays of raw counts, ratios, background correct ratios,
#' smoothed ratios and z-scores. It should be an element of output of 
#' \link{smoothRatiosByChromosome}
#' @param assayName character(1). Assay name for 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' @param colName character(1). Column name of metadata of dat or assay of dat 
#' for coverage weight, see
#' \link[GenomicRanges:coverage-methods]{coverage}, 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}.
#' @param con The connection to which data is saved. If this is a character
#' vector, it is assumed to be a filename and a corresponding file connection
#' is created and then closed after exporting the object.
#' If missing, a \link[IRanges:AtomicList-class]{SimpleRleList} will be returned.
#' @param format The format of the output. see \link[rtracklayer]{export}.
#' @param ... Parameters to be passed to \link[rtracklayer]{export}
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import rtracklayer
#' @export
#' @return If con is missing, a \link[IRanges:AtomicList-class]{SimpleRleList} will be returned.
#' Otherwise, nothing is returned.
#' @examples
#' gr <- GRanges("chr1", IRanges(seq_len(100), 201:300), reads=rep(1, 100))
#' myTrackLine <- new("TrackLine", name="my track",
#'                     description="description of my track",
#'                     color=col2rgb("red")[, 1],
#'                     visibility="full")
#' exportSignals(gr, colName="reads", 
#'               con="test.bedGraph", trackLine=myTrackLine)
#' data(triplicate.count)
#' exportSignals(triplicate.count, "counts", 
#'               "G18.subsampled.srt.bam", "test.bw", format="bigWig")
#'
exportSignals <- function(dat, assayName, colName,
                          con, format = "bedGraph", ...) 
{
    stopifnot(inherits(dat, c(
        "GRanges", "RangedSummarizedExperiment")))
    if (is(dat, "GRanges")) 
    {
        gr <- dat
    } else
    {
        ## RangedSummarizedExperiment
        stopifnot(!missing(assayName))
        stopifnot(assayName %in% names(assays(dat)))
        gr <- rowRanges(dat)
        mcols(gr) <- assays(dat)[[assayName]]
    }
    if (!colName %in% colnames(mcols(gr))) 
    {
        stop("colName is not valid.")
    }
    seqlevels(gr) <-
        seqlevels(gr)[seqlevels(gr) %in% unique(seqnames(gr))]
    dat1 <- coverage(gr, weight = mcols(gr)[, colName])
    dat2 <- coverage(gr)
    dat <- dat1 / dat2
    dat[is.na(dat)] <- 0
    if (missing(con)) 
    {
        return(dat)
    } else
    {
        return(export(
            object = dat,
            con = con,
            format = format,
            ...))
    }
}
