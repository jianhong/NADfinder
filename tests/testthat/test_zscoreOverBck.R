test_that("zscoreOverBck works not correct", {
    x <- runif(10000)
    y <- zscoreOverBck(x, backgroundPercentage=1)
    expect_equal(mean(y), 0, tolerance=1e-10)
})
