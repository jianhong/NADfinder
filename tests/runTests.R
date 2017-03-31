require("NADfinder") || stop("unable to load Package:NADfinder")
require("SummarizedExperiment") || 
    stop("unable to load Package:SummarizedExperiment")
require("testthat") || stop("unable to load testthat")

test_check("NADfinder")
