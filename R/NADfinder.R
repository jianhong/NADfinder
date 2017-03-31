#' @title Call peaks for nucleolar-associated domains (NADs) sequencing data
#'
#' @description 
#' Call peaks for two purified nucleoli samples: target and control.
#' It will count the reads for tiles of the genome and NADs,
#' then convert it to ratios.
#' The ratios will be corrected and smoothed. The z-scores is calculated for
#' each counting windows over the background. The peaks will be detected based
#' on z-scores.
"_PACKAGE"
