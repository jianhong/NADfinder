#' plot signals with ideograms
#'
#' Plot signals with ideograms for \link[GenomicRanges:GRangesList-class]{GRangesList}.
#'
#' @param ideo Output of \link[trackViewer]{loadIdeogram}.
#' @param grList A \link[GenomicRanges:GRangesList-class]{GRangesList} of data to plot.
#' @param mcolName Column name of metadata of GRangesList for plotting.
#' @param ... Parameters to pass to \link[trackViewer]{ideogramPlot}
#' @importFrom trackViewer ideogramPlot
#' @export
#' @return Invisible argument list for \link[trackViewer]{ideogramPlot}.
#' @examples
#'
#' library(trackViewer)
#' #ideo <- loadIdeogram("mm10")
#' ideo <- readRDS(system.file("extdata", "ideo.mm10.rds",
#'                              package = "NADfinder"))
#' gr1 <- gr2 <- ideo
#' mcols(gr1) <- DataFrame(score=runif(length(gr1)))
#' mcols(gr2) <- DataFrame(score=runif(length(gr2)))
#' grList <- GRangesList(gr1, gr2)
#' plotSig(ideo, grList, mcolName="score", layout=list("chr1"))
#'
#'

plotSig <- function(ideo, grList, mcolName, ...){
    stopifnot(is(grList, "GRangesList"))
    args <- list(...)
    ylabs <- lapply(seqlevels(ideo), function(.ele){
        c(.ele, names(grList))
    })
    names(ylabs) <- seqlevels(ideo)
    if(length(args$parameterList)){
        if(!length(args$parameterList$dataColumn))
            args$parameterList$dataColumn <- mcolName
        if(!length(args$parameterList$ylabs))
            args$parameterList$ylabs <- ylabs
    }else{
        args$parameterList <- list(dataColumn=mcolName,
                                   ylabs=ylabs)
    }
    args$ideo <- ideo
    args$dataList <- grList
    do.call(ideogramPlot, args = args)
    return(invisible(args))
}
