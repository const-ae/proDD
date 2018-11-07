context("calculate_distance")

test_that("distance calculations works", {

    data <- matrix(rnorm(1000), nrow=4, ncol=250)
    # Without missing data dist_approx behaves the same as dist
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA)$mean), c(dist(data)))
    expect_equal(c(dist_approx(data, mu_mis=NA, var_mis=NA)$var), rep(0, 6))


    # Method also works if no mu_mis and var_mis are provided
    expect_equal(c(dist_approx(data)$mean), c(dist(data)))


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
    expect_silent({
        res <- dist_approx(t(data$X), hyperparameters = list(
            mu0, sigma20, mean(rho), mean(zeta)
        ))
        res2 <- dist_approx(t(data$X), hyperparameters = list(
            mu0, sigma20, rho, zeta
        ))
        res3 <- dist_approx(t(data$X), experimental_design = experimental_design)
    })


    dist(t(data$t_X))

})
