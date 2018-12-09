

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
    expect_error(msply_dbl("asd", function(x) c(x, x)))
    expect_equal(msply_dbl(TRUE, function(x) c(x, x)), matrix(c(1,1), nrow=1))
})


test_that("mply_dbl handles matrix input", {

    mat <- matrix(1:10, nrow=5)

    expect_equal(mply_dbl(mat, sum, ncol=1), matrix(rowSums(mat), ncol=1))
    expect_equal(mply_dbl(mat, function(x) c(x[1]+x[2], x[1]*x[2]), ncol=2),
                 matrix(c(mat[,1] + mat[,2], mat[,1] * mat[,2]), ncol=2))

    expect_equal(msply_dbl(mat, function(x) c(x[1]+x[2], x[1]*x[2])),
                 matrix(c(mat[,1] + mat[,2], mat[,1] * mat[,2]), ncol=2))


})

test_that("mply_dbl can handle ncol=0", {

    input <- 1:10
    res <- mply_dbl(input, function(x) {numeric(0)}, ncol=0)
    expect_equal(dim(res), c(10, 0))
    res <- msply_dbl(input, function(x) {numeric(0)})
    expect_equal(dim(res), c(10, 0))

    mat <- matrix(1:10, nrow=5)


    expect_equal(dim(msply_dbl(input, function(x) numeric(0))), c(10, 0))
    expect_equal(dim(msply_dbl(mat, function(x) numeric(0))), c(5, 0))

    expect_equal(dim(mply_dbl(input, function(x) numeric(0), ncol=0)), c(10, 0))
    expect_equal(dim(mply_dbl(mat, function(x) numeric(0), ncol=0)), c(5, 0))

})


test_that("mply_dbl can handle empty inputs", {

    input <- numeric(0)

    expect_equal(dim(mply_dbl(input, function(x) x, ncol=3)), c(0, 3))
    expect_equal(dim(msply_dbl(input, function(x) x)), c(0, 0))

    mat <- matrix(numeric(0), nrow=0, ncol=5)

    expect_equal(dim(mply_dbl(mat, function(x) x, ncol=2)), c(0, 2))
    expect_equal(dim(msply_dbl(mat, function(x) x)), c(0, 0))

})

