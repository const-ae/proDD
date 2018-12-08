context("test-prodd_parameters")

test_that("printing looks beautiful", {
    n_rows <- 100
    n_cols <- 6
    n_conditions <- 2
    experimental_design <- rep(1:n_conditions, each=n_cols/n_conditions)
    hyper_params <- list(eta=0.2234234, nu=3.31213, mu0=18, sigma20=3.123, rho=rnorm(n_cols, mean=18), zeta=rep(-1, times=n_cols))
    feat_params <- list(mup=matrix(rnorm(n_conditions * n_rows, mean=18), nrow=n_rows, ncol=n_conditions),
                        sigma2p=rep(0.4, times=n_rows),
                        sigma2mup = matrix(0.1, nrow=n_rows, ncol=n_conditions))
    params <- list(hyper_params=hyper_params,
                   feature_params=feat_params,
                   experimental_design=experimental_design,
                   error=1e-5, converged=TRUE)
    class(params) <- "prodd_parameters"
    params

    expect_true(is_valid_prodd_parameters(params))

    trans_params <- transform_parameters(params, c(1, 3, 2, 2, NA, 1))
    expect_equal(trans_params$feature_params$mup[, "3"], params$feature_params$mup[, 1])
    expect_equal(trans_params$experimental_design, factor(c(1, 3, 2, 2, 1)))
    expect_equal(trans_params$error, params$error)
    expect_equal(trans_params$feature_params$mup[, "1"], trans_params$feature_params$mup[, "2"])

    nrep <- 2
    mu0 <- 20
    sigma20 <- 10
    nu <- 3
    eta <- 0.3
    rho <- rnorm(nrep*2, 18, sd=2)
    zeta <- rep(-0.1, nrep * 2)
    experimental_design <- c(rep(1, nrep), rep(2, nrep))
    data <- generate_synthetic_data(n_rows=100, n_replicates=nrep, n_conditions=2,
                                    mu0=mu0, sigma20=sigma20, rho=rho, zeta=zeta,
                                    frac_changed = 0.1)

    fit <- fit_parameters(data$X, experimental_design)
    fit
    expect_equal(class(format(fit)), "character")
})





