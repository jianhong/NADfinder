test_that("smoothRatiosbyChromosome works not correct", {
    data(single.count)
    se <- single.count
    se <- log2se(se, 
                 nucleosomeCols = "nucleosome.bam", 
                 genomeCols = "genome.bam")
    se <- smoothRatiosByChromosome(se, N=100)
    chr18 <- assays(se[["chr18"]])
    len <- length(se[[1]])
    # original ratios is higher in 5ends and lower in 3ends
    r5 <- quantile(chr18$ratio[1:floor(len/2)])[2]
    r3 <- quantile(chr18$ratio[ceiling(len/2):len])[2]
    expect_true(r3 < 0)
    # background corrected ratios
    bcr5 <- quantile(chr18$bcRatio[1:floor(len/2)])[2]
    bcr3 <- quantile(chr18$bcRatio[ceiling(len/2):len])[2]
    expect_true(bcr3 > 0)
    ## background should be correct, how to check?
    # fit line
    data <- data.frame(x=(1:len)/1e4, y=chr18$ratio[, 1], z=chr18$bcRatio[, 1])
    lm1 <- lm(y ~ x, data=data)
    lm2 <- lm(z ~ x, data=data)
    expect_true(abs(lm2$coefficients[[2]]) < abs(lm1$coefficients[[2]]))
})
