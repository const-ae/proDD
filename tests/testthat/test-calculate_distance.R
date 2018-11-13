context("calculate_distance")

test_that("distance calculations works", {

    data <- matrix(rnorm(1000), nrow=4, ncol=250)
    # Without missing data dist_approx behaves the same as dist
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA)$mean), c(dist(data)))
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA)$var), rep(0, 6))



    nrep <- 6
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

    fit <- fit_hyperparameters(data$X, experimental_design)


    mis <- find_approx_for_missing(data$X, fit, experimental_design=experimental_design)

    expect_silent({
        dist_approx(t(data$X), params=fit)
        dist_approx(t(data$X), params=fit, experimental_design = rep(1, times=ncol(data$X)))
        dist_approx(t(data$X), hyper_params = fit$hyper_params, experimental_design = rep(1, ncol(data$X)))
        dist_approx(t(data$X), hyper_params = fit$hyper_params, experimental_design = experimental_design)
        dist_approx(t(data$X), feature_params = fit$feature_params, rho=rho, zeta=zeta, experimental_design = experimental_design)
        dist_approx(data$X, hyper_params = fit$hyper_params, experimental_design = experimental_design)
        dist_approx(data$X, feature_params = fit$feature_params, rho=rho, zeta=zeta, experimental_design = experimental_design)
        dist_approx(data$X, mis$mu_mis, mis$var_mis)
        dist_approx(t(data$X), t(mis$mu_mis), t(mis$var_mis))
        dist_approx(t(data$X), 16, 0.1)
    })


})
