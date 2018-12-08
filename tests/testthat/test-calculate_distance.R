context("calculate_distance")

test_that("distance calculations works", {

    data <- matrix(rnorm(1000), nrow=4, ncol=250)
    # Without missing data dist_approx behaves the same as dist
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA, by_sample=FALSE)$mean), c(dist(data)))
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA, by_sample=FALSE)$var), rep(0, 6))



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

    fit <- fit_parameters(data$X, experimental_design, dropout_curve_calc = "global", max_iter = 2)


    mis <- find_approx_for_missing(data$X, fit, experimental_design=experimental_design)

    dmat <- dist_approx(data$X, params=fit)$mean
    expect_equal(nrow(as.matrix(dmat)), 12)

    dmat2 <- dist_approx(data$X, params=fit, by_sample = FALSE)$mean
    expect_equal(nrow(as.matrix(dmat2)), 100)

    dmat3 <- dist_approx(data$X, mu_mis=mis$mu_mis, var_mis=mis$var_mis)$mean
    expect_equal(dmat, dmat3)

    dmat4 <- dist_approx(data$X, mu_mis=mis$mu_mis, var_mis=mis$var_mis, by_sample=FALSE)$mean
    dmat5 <- dist_approx(t(data$X), mu_mis=t(mis$mu_mis), var_mis=t(mis$var_mis))$mean

    expect_equal(dmat4, dmat2)
    expect_equal(dmat4, dmat5)

    expect_silent(dist_approx(data$X, mu_mis=16, var_mis=0.1))

    expect_error(dist_approx(t(data$X), fit))


})
