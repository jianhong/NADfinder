test_that("log2ratio works not correct", {
    # log2ratio=log2(A) - log2(B)
    expect_equal(log2ratio(1, 1), 0)
    expect_equal(log2ratio(0, 1), -1)
    expect_equal(log2ratio(1:10, 10:1), log2(1:10)-log2(10:1))
    # B is 0
    expect_false(is.infinite(log2ratio(1, 0)))
})
