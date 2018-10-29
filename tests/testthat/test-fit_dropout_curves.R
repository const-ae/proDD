context("test-fit_dropout_curves")






test_that("fitting a global sigmoids work", {
    set.seed(1)
    rho <- rep(19, times=5)
    zeta <- rep(-0.5, times=5)
    data <- generate_synthetic_data(n_rows = 1000, n_replicates=5,
                                    rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_global_dropout_curves(data$X, data$mus, data$sigmas2,
                       experimental_design=rep(1, 5), prev_rho=rnorm(5, 18, 1),
                       prev_zeta=rep(-1, 5))
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.1)
    expect_equal(res$zeta, zeta, tolerance=0.1)

    # More data means better fit
    set.seed(1)
    rho <- rep(18.7, times=10)
    zeta <- rep(-0.5, times=10)
    data <- generate_synthetic_data(n_rows = 1e4, n_replicates=10,
                                    rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_global_dropout_curves(data$X, data$mus, data$sigmas2,
                       experimental_design=rep(1, 10), prev_rho=18,
                       prev_zeta=-1)
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.01)
    expect_equal(res$zeta, zeta, tolerance=0.01)
})


test_that("fitting a global sigmoids work", {
    set.seed(1)
    rho <- rep(19, times=5)
    zeta <- rep(-0.5, times=5)
    data <- generate_synthetic_data(n_rows = 1000, n_replicates=5,
                                    rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_sample_dropout_curves(data$X, data$mus, data$sigmas2,
                           experimental_design=rep(1, 5), prev_rho=rnorm(5, 18, 1),
                           prev_zeta=rep(-1, 5))
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.2)
    expect_equal(res$zeta, zeta, tolerance=0.2)

    # More data means better fit
    set.seed(1)
    rho <- rep(18.7, times=5)
    zeta <- rep(-0.5, times=5)
    data <- generate_synthetic_data(n_rows = 1e4, n_replicates=5,
                                    rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_sample_dropout_curves(data$X, data$mus, data$sigmas2,
                             experimental_design=rep(1, 5), prev_rho=rnorm(5, 18, 1),
                             prev_zeta=rep(-1, 5))
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.05)
    expect_equal(res$zeta, zeta, tolerance=0.05)
})




test_that("fitting a global scale and individual rhos work", {
    set.seed(1)
    rho <- rnorm(5, mean=19, sd=1)
    zeta <- rep(-0.5, times=5)
    data <- generate_synthetic_data(n_rows = 1000, n_replicates=5,
                                   rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_global_scale_dropout_curves(data$X, data$mus, data$sigmas2,
                        experimental_design=rep(1, 5), prev_rho=rnorm(5, 18, 1),
                        prev_zeta=-1)
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.1)
    expect_equal(res$zeta, zeta, tolerance=0.1)

    # More data means better fit
    data <- generate_synthetic_data(n_rows = 1e4, n_replicates=5,
                                    rho=rho, zeta=zeta)
    expect_silent({
        res <- fit_global_scale_dropout_curves(data$X, data$mus, data$sigmas2,
                        experimental_design=rep(1, 5), prev_rho=rnorm(5, 18, 1),
                        prev_zeta=-1, epsilon=1e-5)
    })
    expect_true(res$converged)
    expect_equal(res$rho, rho, tolerance=0.1)
    expect_equal(res$zeta, zeta, tolerance=0.01)
})
