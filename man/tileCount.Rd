% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tileCount.R
\name{tileCount}
\alias{tileCount}
\title{Perform overlap queries between reads and genome by windows}
\usage{
tileCount(
  reads,
  genome,
  excludeChrs = c("chrM", "M", "Mt", "MT"),
  windowSize = 50000,
  step = 10000,
  mode = IntersectionNotStrict,
  dataOverSamples = FALSE,
  ...
)
}
\arguments{
\item{reads}{A \link[GenomicRanges:GRanges-class]{GRanges},
\link[GenomicRanges:GRangesList-class]{GRangesList} (should be one read per list element),
\link[GenomicAlignments:GAlignments-class]{GAlignments},
\link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList},
\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs} or
\link[Rsamtools:BamFile-class]{BamFileList} object that represents the data to be
counted by \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}. If reads are more than 1 bam files,
it should be a vector of character with full path, otherwise current working directory 
is the default directory.}

\item{genome}{A BSgenome object from/on which to get/set the sequence and metadata information.}

\item{excludeChrs}{A vector of string: chromosomes/scaffolds of no interest for NAD analysis.
see \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
default is countByOverlaps, alia of countOverlaps(features, reads, ignore.strand=ignore.strand)}

\item{windowSize}{numeric(1) or integer(1). Size of the windows.}

\item{step}{numeric(1) or integer(1). Step of generating silding windows.}

\item{mode}{One of the pre-defined count methods.}

\item{dataOverSamples}{logical(1). Data over several samples when use
\link[GenomicRanges:GRangesList-class]{GRangesList} as input.}

\item{...}{Additional arguments passed to
\code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.}
}
\value{
A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} object.
The assays slot holds the counts, rowRanges holds the annotation from the
sliding widows of genome.
metadata contains lib.size.chrom for holding chromosome-level sequence depth
}
\description{
tileCount extends \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps} by finding coverage for
each fixed window in the whole genome
}
\examples{
if (interactive())
{
    fls <- list.files(system.file("extdata", package="NADfinder"),
    recursive=FALSE, pattern="*bam$", full=TRUE)
    names(fls) <- basename(fls)
    if (!require(BSgenome.Mmusculus.UCSC.mm10))
    {
        if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
        BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
        library(BSgenome.Mmusculus.UCSC.mm10)
    }
    se <- tileCount(reads = fls, 
                    genome = Mmusculus,
                    excludeChrs = c("chrM", paste0("chr", c(1:17,19)), 
                                    "chrX", "chrY"), 
                    windowSize=50000, step=10000)
}


}
\author{
Jianhong Ou, Haibo Liu, Herve Pages and Julie Zhu
}
