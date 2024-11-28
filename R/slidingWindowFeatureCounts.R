#' Summarize reads for different genomic features
#'
#' Summarize reads in alignment files, SAM or BAM, to different genomic regions,
#' such as genic regions, intergenic regions, exonic regions, intronic regions,
#' rRNA genes, mitochrondrial genome, chloroplast genome (only for plants), and
#' gene-level exonic regions using [Rsubread::featureCounts()].
#'
#' @param bamfiles A named character(n), specifying paths to a set of BAM files.
#' @param saf A data frames containing sliding windows across the genome and
#'   individual chromosome-level annotation in the SAF format. It should be the
#'   output of the [generateGenomeSAF()] function.
#' @param isPairedEnd A logical(1), indicating whether the sequencing data is
#'   paired-end or not.
#' @param threads An integer(1), number of threads for
#'   [Rsubread::featureCounts()] calling.
#' @param checkFragLength A logical(1), indicating if the two ends from the
#'   same fragment are required to satisify the fragment length criteria
#'   before the fragment can be assigned to a feature or meta-feature. This
#'   parameter is only appliable when isPairedEnd is TRUE. The fragment
#'   length criteria are specified via minFragLength and maxFragLength.
#' @param requireBothEndsMapped A logical(1), indicating if both ends from
#'   the same fragment are required to be successfully aligned before the
#'   fragment can be assigned to a feature or meta-feature. This parameter
#'   is only appliable when isPairedEnd is TRUE.
#' @param minFragLength An integer(1), giving the minimum fragment length for
#'   paired-end reads. 50 by default.
#' @param maxFragLength An integer(1), giving the maximum fragment length for
#'   paired-end reads. 600 by default. minFragLength and maxFragLength are only
#'   applicable when isPairedEnd is TRUE. Note that when a fragment spans two
#'   or more exons, the observed fragment length might be much bigger than the
#'   nominal fragment length.
#' @param minMQS An integer(1), giving the minimum mapping quality score a read
#'   must satisfy in order to be counted. For paired-end reads, at least one end
#'   should satisfy this criteria. 0 by default.
#' @param verbose A logical(1) vector, indicating if verbose information for
#' debugging will be generated. This may include information such as unmatched
#' chromosomes/contigs between reads and annotation.
#'
#' @return A list of lists, each of the sublist contains an output of the
#'   [Rsubread::featureCounts()] function for each type of genomic features.
#'   For more details of each sublist, see the *value* section of the
#'   documentation for the [Rsubread::featureCounts()] function.
#'
#' @importFrom Rsubread featureCounts
#' @export
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="NADfinder"),
#'                        recursive=FALSE,
#'                        pattern="*bam$",
#'                        full=TRUE)
#' names(bamfiles) <- c("G18", "N18")
#' genomeSW_SAF <-
#'     generateGenomeSAF(bamfile = bamfiles[1],
#'                       windowSize = 50000L,
#'                       step = 5000L,
#'                       chrExcluded = c("chrM", "MT"))
#'
#' counts <- swFeatureCounts(bamfiles = bamfiles,
#'                           saf = genomeSW_SAF,
#'                           isPaired = TRUE)
#'

swFeatureCounts <- function(bamfiles = NULL,
                            saf = NULL,
                            isPairedEnd = TRUE,
                            threads = 1L,
                            checkFragLength = FALSE,
                            requireBothEndsMapped = FALSE,
                            minFragLength = 50L,
                            maxFragLength = 1000L,
                            minMQS = 0,
                            verbose = FALSE)
    {
       if (any(!file.exists(bamfiles)))
       {
           stop("Not all BAM files exist!")
       } else if (is.null(names(bamfiles))) {
           stop("BAM files must be a named vector of character!")
       }
       saf_colnames <- c("GeneID", "Chr", "Start", "End", "Strand")
       if (is.null(saf))
       {
           stop("saf is required!")
       }  else if (!is.data.frame(saf) ||
                   !all(saf_colnames == colnames(saf[[1]]))) {
           stop("saf must be a data frame!")
       }
       res <- featureCounts(
           files = bamfiles,
           # annotation
           annot.ext = saf,
           isGTFAnnotationFile = FALSE,
           # overlap between reads and features
           allowMultiOverlap = TRUE,
           minOverlap = 1,
           fracOverlap = 0,
           fracOverlapFeature = 0,
           largestOverlap = TRUE,
           nonOverlap = NULL,
           nonOverlapFeature = NULL,
           # multi-mapping reads
           countMultiMappingReads = FALSE,
           # fractional counting
           fraction = FALSE,
           # read filtering
           minMQS = minMQS,
           splitOnly = FALSE,
           nonSplitOnly = FALSE,
           primaryOnly = FALSE,
           ignoreDup = FALSE,
           # strandness
           strandSpecific = 0,
           genome = NULL,
           # parameters specific to paired end reads
           isPairedEnd = isPairedEnd,
           countReadPairs = TRUE,
           requireBothEndsMapped = requireBothEndsMapped,
           checkFragLength = checkFragLength,
           minFragLength = minFragLength,
           maxFragLength = maxFragLength,
           countChimericFragments = TRUE,
           autosort = TRUE,
           # number of CPU threads
           nthreads = threads,
           # read group
           byReadGroup = FALSE,
           # report assignment result for each read
           reportReads = NULL,
           reportReadsPath = NULL,
           # miscellaneous
           tmpDir = ".",
           verbose = verbose
       )
       colnames(res$counts) <- names(bamfiles)

       ## chromosome-level sequencing depth
       lib.size.chrom <- res$counts[grepl("genome_", rownames(res$counts)), ,
                                    drop = FALSE]
       rownames(lib.size.chrom) <-
           res$annotation$Chr[grepl("genome", res$annotation$GeneID)]

       slidingwindows_saf <-
           res$annotation[!grepl("genome", res$annotation$GeneID), ]
       rowRanges <-
           makeGRangesFromDataFrame(slidingwindows_saf,
                                    seqinfo = seqinfo(BamFile(bamfiles[1])),
                                    seqnames.field="Chr",
                                    start.field="Start",
                                    end.field= "End",
                                    strand.field="Strand",
                                    starts.in.df.are.0based=FALSE)
       counts <- res$counts[!grepl("genome_", rownames(res$counts)), ,
                            drop = FALSE]
       rse <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                   rowRanges = rowRanges)

       metadata(rse)$lib.size.chrom <- lib.size.chrom
       rse
    }
