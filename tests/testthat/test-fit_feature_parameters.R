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



test_that("unregularized mean estimates work", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 500, n_replicates = N_rep, n_conditions = 2,
                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)

    mu_unreg <- fit_unregularized_feature_means(data$X, data$mus, mu0=20,
                                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                                    experimental_design=rep(1:2, each=N_rep))


    expect_gt(diag(cor(data$mus, mu_unreg, use = "complete.obs"))[1], 0.99)
    expect_gt(diag(cor(data$mus, mu_unreg, use = "complete.obs"))[2], 0.99)
})


test_that("inferring the location prior works", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 1000, n_replicates = N_rep, n_conditions = 2,
                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    loc_prior <- fit_location_prior(data$X, data$mus, rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    experimental_design=rep(1:2, each=N_rep))

    expect_equal(loc_prior$mu0, 20, tolerance=0.1)
    expect_equal(loc_prior$sigma20, 10, tolerance=0.1)
})


test_that("inferring the variance prior works", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 100, n_replicates = N_rep, n_conditions = 2,
                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    var_prior <- fit_variance_prior(data$X, rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    experimental_design=rep(1:2, each=N_rep))

    expect_equal(var_prior$var_prior, 0.3, tolerance=0.1)
    expect_equal(var_prior$df_prior, 10, tolerance=1)
})




test_that("everything ties together", {
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 100, n_replicates = N_rep, n_conditions = 2,
                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    result <- fit_hyperparameters(data$X, experimental_design=rep(1:2, each=N_rep),
                                  dropout_curve_calc = "global", verbose=FALSE)

    expect_equal(result$hyper_params$eta, 0.3, tolerance=0.1)
    expect_equal(result$hyper_params$nu, 10, tolerance=1)
    expect_equal(result$hyper_params$mu0, 20, tolerance=1)
    expect_equal(result$hyper_params$sigma20, 10, tolerance=1)
    expect_equal(result$hyper_params$rho[1], 18, tolerance=1)
    expect_equal(result$hyper_params$zeta[1], -1, tolerance=0.2)

})


test_that("predicting features for new data works", {
    experimental_design <- rep(LETTERS[1:2], each=4)
    data <- generate_synthetic_data(n_rows = 200, experimental_design,
                                    rho=18, zeta=-1.2,
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    result <- fit_hyperparameters(data$X, experimental_design, frac_subsample = 1,
                                  dropout_curve_calc = "global", verbose=TRUE)


    new_res <- predict_feature_parameters(data$X, experimental_design, result$hyper_params,
                                          verbose=TRUE)


    expect_equal(result$feature_params$mup, new_res$mup, tolerance=0.01)
    expect_equal(result$feature_params$sigma2p, new_res$sigma2p, tolerance=0.01)
    expect_equal(result$feature_params$sigma2mup, new_res$sigma2mup, tolerance=0.01)

})








