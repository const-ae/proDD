

context("mply_dbl")

library(proDD)

test_that("mply_dbl can handle simple inputs", {
    res <- mply_dbl(1:10, function(x) c(x*2, x^2), ncol=2)
    expect_equal(res, cbind((1:10)*2, (1:10)^2))

    res2 <- mply_dbl(1:10, identity, ncol=1)
    expect_equal(res2, matrix(1:10, ncol=1))

})

test_that("msply_dbl can handle simple inputs", {
    res <- msply_dbl(1:10, function(x) c(x*2, x^2))
    expect_equal(res, cbind((1:10)*2, (1:10)^2))


    res2 <- mply_dbl(1:10, identity)
    expect_equal(res2, matrix(1:10, ncol=1))
})

test_that("mply_dbl handles types correctly", {
    expect_equal(mply_dbl(TRUE, function(x) c(x, x*2), ncol=2), matrix(1:2, nrow=1))
    expect_error(mply_dbl("asd", function(x) c(x, x), ncol=2))
})
