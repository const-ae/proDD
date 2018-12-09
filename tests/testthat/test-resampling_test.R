context("test-resampling_test")

test_that("resampling test works", {
    set.seed(1)
    x <- rnorm(1000, mean=0, sd=1)
    y <- rnorm(1000, mean=2, sd=2)
    expect_equal(resampling_test(x, y, alternative="less"),   mean(outer(x, y, `>`)))
    expect_equal(resampling_test(x, y, nmax=30000, alternative="less"),   mean(outer(x, y, `>`)),
                 tolerance=0.1)


    expect_lt(resampling_test(x, y, alternative="less"), 0.5)

    expect_equal(resampling_test(x, y, alternative="less"),
                 1-resampling_test(x, y, alternative="greater"))

    expect_equal(2 * resampling_test(x, y, alternative="less"),
                 resampling_test(x, y, alternative="two.sided"))

    expect_equal(resampling_test(x, y, mu=3),
                 resampling_test(x, y + 3))

    expect_equal(resampling_test(x, 2, alternative="less"),
                 pt((mean(x) - 2) / sd(x), df=length(x)-1), tolerance=0.1)
})


test_that("test_diff works", {


    experimental_design <- rep(c("A", "B"), each=3)
    data <- generate_synthetic_data(n_rows=100, experimental_design)
    params <- fit_parameters(data$X, experimental_design)
    posteriors <- sample_protein_means(data$X, params, verbose=FALSE, control=list(adapt_delta=0.95),
                                           niter=1000, nchains=1)


    res <- test_diff(posteriors$A, posteriors$B)
    expect_equal(colnames(res), c("name", "pval", "adj_pval", "diff", "diff_quantile_0.025", "diff_quantile_0.975", "mean"))
    expect_equal(dim(res), c(100, 7))


})

