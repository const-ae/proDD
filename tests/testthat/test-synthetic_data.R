context("test-synthetic_data")

test_that("generating new data works", {

    d1 <- generate_synthetic_data(n_rows=10, n_replicates=4)
    expect_equal(ncol(d1$X), 4)
    expect_equal(ncol(d1$t_X), 4)
    expect_equal(ncol(d1$mus), 1)
    expect_equal(length(d1$sigmas2), 10)
    expect_equal(length(d1$changed), 10)


    d2 <- generate_synthetic_data(n_rows=10, n_replicates=4, n_conditions = 1)
    expect_equal(ncol(d2$X), 4)

    d3 <- generate_synthetic_data(n_rows=10, n_replicates=4, n_conditions = 2)
    expect_equal(ncol(d3$X), 8)
    expect_equal(ncol(d3$mus), 2)

    d4 <- generate_synthetic_data(n_rows=10, n_replicates=c(2, 7))
    expect_equal(ncol(d4$X), 9)
    expect_equal(ncol(d4$mus), 2)

    d5 <- generate_synthetic_data(n_rows=10, n_replicates=c(2, 7),
                                  frac_changed = 0.5)
    expect_equal(sum(d5$changed), 5)

    # X contains the same values as t_X except for the NA's
    expect_equal(d1$X[! is.na(d1$X)],d1$t_X[! is.na(d1$X)])

    # If the sigmoid is far out, none of the values is missing
    d6 <- generate_synthetic_data(n_rows=100, n_replicates = 10, rho=-17)
    expect_true(all(! is.na(d6$X)))

    # Having some weird inputs
    expect_silent(generate_synthetic_data(n_rows=0, n_replicates = 10))
    expect_error(generate_synthetic_data(n_rows=10, n_replicates = c(1, 0, 3)))

    experimental_design <- c("c1","c2","c1","c2","c1")
    d7 <- generate_synthetic_data(n_rows=100, experimental_design = experimental_design)
    expect_equal(colnames(d7$X), c("c1-1","c2-1","c1-2","c2-2","c1-3"))
})



test_that("as_replicate works as expected", {

    expect_equal(as_replicate(c(1, 2, 1, 2)), c(1,1,2,2))
    expect_equal(as_replicate(c(1, 2, 2, 2, 2, 3, 1)), c(1,1,2,3,4,1,2))
    expect_equal(as_replicate(c(1, 2, 3, 4, 2, 1, 3, 3, 4)), c(1,1,1,1,2,2,2,3,2))

    x <- c("a", "b", "a", "b", "b", "d")
    expect_equal(as_replicate(x), c(1,1,2,2,3,1))
    x <- as.factor(c("c", "a", "b", "a", "b", "b", "B"))
    expect_equal(as_replicate(x), c(1, 1,1,2,2,3,1))


})
