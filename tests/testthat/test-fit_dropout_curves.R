context("test-fit_dropout_curves")

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
