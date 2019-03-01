context("test-s4_objects")

test_that("Support for SummarizedExperiment works", {
    experimental_design <- rep(LETTERS[1:2], each=4)
    data <- generate_synthetic_data(n_rows = 35, experimental_design,
                                    rho=18, zeta=-1.2,
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    library(SummarizedExperiment)

    sample_info <- data.frame(
        condition = experimental_design,
        replicate = proDD:::as_replicate(experimental_design),
        row.names = colnames(data$X)
    )
    protein_info <- data.frame(
        ID = seq_len(35),
        changed = data$changed,
        variance = data$sigmas2,
        mu_a = data$mus[, "A"],
        mu_b = data$mus[, "B"],
        row.names = rownames(data$X)
    )
    se <- SummarizedExperiment(data$X, colData = sample_info, rowData = protein_info)
    norm_se <- median_normalization(se)
    expect_equal(assay(norm_se), median_normalization(data$X))
    set.seed(1)
    se_params <- fit_parameters(se, "condition")
    set.seed(1)
    std_params <- fit_parameters(data$X, experimental_design)
    expect_equal(se_params, std_params, tolerance = 1e-3)
    se_res <- sample_protein_means(norm_se, se_params, verbose = FALSE)
    se_dist <- dist_approx(se, se_params)
    expect_equal(se_dist, dist_approx(data$X, std_params), tolerance = 1e-3)
})

test_that("Support for MSnSet works", {
    experimental_design <- rep(LETTERS[1:2], each=4)
    data <- generate_synthetic_data(n_rows = 35, experimental_design,
                                    rho=18, zeta=-1.2,
                                    nu=10, eta=0.3, mu0=20, sigma20=10)
    library(MSnbase)

    sample_info_adf <- AnnotatedDataFrame(data.frame(
        condition = experimental_design,
        replicate = as_replicate(experimental_design),
        row.names = colnames(data$X)
    ))
    protein_info_adf <- AnnotatedDataFrame(data.frame(
        ID = seq_len(35),
        changed = data$changed,
        variance = data$sigmas2,
        mu_a = data$mus[, "A"],
        mu_b = data$mus[, "B"],
        row.names = rownames(data$X)
    ))

    ms <- MSnSet(exprs = data$X, pData = sample_info_adf, fData = protein_info_adf)
    norm_ms <- median_normalization(ms)
    expect_equal(exprs(norm_ms), median_normalization(data$X))
    set.seed(1)
    ms_params <- fit_parameters(ms, "condition")
    set.seed(1)
    std_params <- fit_parameters(data$X, experimental_design)
    expect_equal(ms_params, std_params, tolerance = 1e-3)
    ms_res <- sample_protein_means(norm_ms, ms_params, verbose = FALSE)
    ms_dist <- dist_approx(ms, ms_params)
    expect_equal(ms_dist, dist_approx(data$X, std_params), tolerance = 1e-3)
})



