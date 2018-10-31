context("calculate_posterior")

test_that("posterior calculation works", {
    set.seed(1)
    N_rep <- 10
    data <- generate_synthetic_data(n_rows = 100, n_replicates = N_rep, n_conditions = 2,
                                    frac_changed = 0,
                                    rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                    nu=10, eta=0.3, mu0=20, sigma20=10)

    posterior <- calculate_posterior(data$X, mu0=20, sigma20=10, nu=10, eta=0.3,
                                     rho=rep(18, N_rep*2), zeta=rep(-1, N_rep*2),
                                     experimental_design=rep(1:2, each=N_rep),
                                     niter=1000, verbose=FALSE)

    pval <- sapply(1:100, function(idx){
        resampling_test(posterior$`1`[idx, ], posterior$`2`[idx, ])
    })

    expect_gt(ks.test(pval, "punif")$p.value, 0.05)
})
