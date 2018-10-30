context("test-fit_feature_parameters")

test_that("fit_feature_variances works", {

    # If the mean is right of rho, and the only value is missing, then the variance
    # is the mode of the scaled inverse chisq.
    expect_equal(fit_feature_variances(matrix(NA, nrow=1, ncol=1), mup=matrix(0, nrow=1),
                          rho=rep(18, 1), zeta=rep(-1, 1), nu=10, eta=0.3,
                          experimental_design = rep(1, 1)),
                 0.3 * 10 / 12, tolerance=1e-5)

    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 1e3, n_replicates = N_rep,
                                    rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                    nu=10, eta=0.3)

    vars <- fit_feature_variances(data$X, data$mus, rho=rep(18, N_rep),
                                  zeta=rep(-1, N_rep), nu=1, eta=0.3,
                                  experimental_design = rep(1, N_rep))
    expect_true(! any(is.na(vars)))

    # There is no perfect, but at least decent correlation
    expect_gt(cor(data$sigmas2, vars), 0.5)
    # The moderation reduces the squared error
    simple_vars <- matrixStats::rowVars(data$X, na.rm=TRUE)
    expect_lt(sum((data$sigmas2[! is.na(simple_vars)] - vars[! is.na(simple_vars)])^2),
              sum((data$sigmas2[! is.na(simple_vars)] - simple_vars[! is.na(simple_vars)])^2))

})



test_that("fit_feature_means works", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 1e3, n_replicates = N_rep,
                                    rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)

    sigma2p <- fit_feature_variances(data$X, data$mus, rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                     nu=10, eta=0.3, experimental_design = rep(1, N_rep))

    sigma2mup <- fit_feature_mean_uncertainties(data$X, rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                                nu=10, eta=0.3, mu0=20, sigma20=10,
                                                experimental_design = rep(1, N_rep))

    mup <- fit_feature_means(data$X, sigma2p=sigma2p, sigma2mus=sigma2mup,
                             rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                             nu=10, eta=0.3, mu0=20, sigma20=10,
                             experimental_design = rep(1, N_rep))

    nobs <- rowSums(! is.na(data$X))
    expect_gt(cor(data$mus[nobs != 0], mup[nobs != 0]), 0.99)


})



test_that("unregularized variance estimates work", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 500, n_replicates = N_rep,
                                    rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)

    sigma2_unreg <- fit_unregularized_feature_variances(data$X,
                                        rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                        experimental_design=rep(1, N_rep))

    nobs <- rowSums(! is.na(data$X))
    expect_equal(which(is.na(sigma2_unreg)), which(nobs <= 1))

    sigma2p <- fit_feature_variances(data$X, data$mus, rho=rep(18, N_rep), zeta=rep(-1, N_rep),
                                     nu=10, eta=0.3, experimental_design = rep(1, N_rep))

    expect_gt(cor(sigma2p, sigma2_unreg, use = "complete.obs"), 0.8)
})
