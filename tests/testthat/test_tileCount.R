test_that("tileCount works not correct", {
    genes <- GRanges(
        seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
        ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600,
                           4000, 7500, 5000, 5400),
                         width=c(rep(500, 3), 600, 900, 500, 300, 900,
                                 300, 500, 500),
                         names=letters[1:11]))
    ## check tileCount
    reads <- GRangesList(A=GRanges(), B=GRanges())
    # 0 reads
    expect_error(tileCount(reads, genes, windowSize=50, step=10))
    # set seqlengths to genome
    seqlengths(genes) <- c("chr2L"=8000, "chr2R"=8000, "chr3L"=6000)
    tc <- tileCount(reads, genes, windowSize=50, step=10)
    # 100 reads, width=10
    reads$A <- GRanges("chr2L", IRanges(1:100, width=10))
    tc <- tileCount(reads$A, genes, windowSize=200, step=100)
    expect_equal(assays(tc)$counts[1:2, 1], c(100, 9))
    expect_error(tileCount(reads, genes, windowSize=200, step=100))
    tc2 <- tileCount(reads, genes, windowSize=200, step=100,
                     dataOverSamples=TRUE)
    expect_equal(assays(tc)$counts[, 1], assays(tc2)$counts[, 1])
    expect_true(all(assays(tc2)$counts[, "B"]==0))
    tc3 <- tileCount(reads, genes, windowSize=50, step=10, dataOverSamples=TRUE)
    expect_equal(max(assays(tc3)$counts[, "A"]), 59)
})
