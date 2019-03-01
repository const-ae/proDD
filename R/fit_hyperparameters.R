








#' Fit the probabilistic dropout parameters
#'
#' The method infers the position and scale of the dropout sigmoids, the
#' location prior of the means and the prior for the variance. In addition it
#' estimates some feature parameters (mean, uncertainty of mean and variance
#' for each protein and condition).
#'
#' @param X the numerical data where each column is one sample and each row
#'   is one protein. Missing values are coded as \code{NA}.
#' @param experimental_design a vector that assignes each sample to one condition.
#'   It has the same length as the number of columns in \code{X}. It can either be
#'   a factor, a character or a numeric vector. Each unique element is one condition.
#'   If \code{X} is a \code{SummarizedExperiment} or an \code{MSnSet} object,
#'   \code{experimental_design} can also be the name of a column in the
#'   \code{colData(X)} or \code{pData(X)}, respectively.
#' @param dropout_curve_calc string that specifies how the dropout curves are
#'   estimated. There are three different modes. "sample": number of curves=
#'   number of samples, "global": number of curves=1, "global_scale": estimate
#'   only a single scale of the sigmoid, but estimate the position per sample.
#'   Default: "sample".
#' @param frac_subsample number between 0 and 1. Often it is not necessary to
#'   consider each protein, but the computation can be significantly sped up
#'   by only considering a subset of the subsets. Default: 1.0.
#' @param n_subsample number between 1 and \code{nrow(X)}. Alternative way to
#'   specify how many proteins are considered for estimating the hyper-parameters.
#'   Default: \code{nrow(X) * frac_subsample}.
#' @param max_iter integer larger than 1. How many iterations are run at most
#'   trying to reach convergence. Default: 10.
#' @param epsilon number larger than 0. How big is the maximum remaining error
#'   for the algorithm to be considered converged. Default: 10^-3
#' @param verbose boolean. Specify how much extra information is printed
#'   while the algorithm is running. Default: \code{FALSE}.
#'
#' @return a list containing the infered parameters. The list is tagged
#'   with the class "prodd_parameters" for simpler handling in downstream
#'   methods
#'
#' @export
fit_parameters <- function(X, experimental_design,
                                dropout_curve_calc=c("sample", "global_scale", "global"),
                                frac_subsample=1.0, n_subsample=round(nrow(X) * frac_subsample),
                                max_iter=10, epsilon=1e-3, verbose=FALSE){

    dropout_curve_calc <- match.arg(dropout_curve_calc, c("sample", "global_scale", "global"))

    experimental_design_fct <- as.factor(experimental_design)
    experimental_design <- as.numeric(experimental_design_fct)

    N_cond <- length(unique(experimental_design))
    if(n_subsample >= nrow(X)){
        X_bckp <- X
        sel <- seq_len(nrow(X))
    }else if(n_subsample > 0 && n_subsample < nrow(X)){
        X_bckp <- X
        sel <- sample(seq_len(nrow(X_bckp)), n_subsample)
        X <- X_bckp[sel, ,drop=FALSE]
    }else{
        stop(paste0("Illegal argument n_subsample it must be larger than 0"))
    }


    #Initialize the parameters
    mup <- mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            mean(X[idx, which(experimental_design == cond)], na.rm=TRUE)
        }, FUN.VALUE=0.0)
    })

    frac_mis <- sum(is.na(X)) / prod(dim(X))
    mup[is.na(mup)] <- quantile(mup, 0.5 * frac_mis, na.rm=TRUE, names=FALSE)
    mu0 <- mean(mup, na.rm=TRUE)

    # Estimate variances
    sigma2p <-  apply(X, 1, function(row){
        res <- sapply(seq_len(N_cond), function(cond){
            x <- row[which(experimental_design == cond)]
            nobs <- sum(! is.na(x))
            if(nobs > 1){
                c((nobs-1) * var(x, na.rm=TRUE), nobs-1)
            }else{
                c(NA, NA)
            }
        })
        sum(res[1, ], na.rm=TRUE) / sum(res[2, ], na.rm=TRUE)
    })
    sigma20 <- mean(sigma2p, na.rm=TRUE)
    sigma2p[is.na(sigma2p)] <- sigma20
    sigma2_unreg <- sigma2p

    # Sigmoid parameters
    zeta <- rep(- sqrt(sigma20),  ncol(X))
    rho <- sapply(seq_len(ncol(X)), function(colidx){
        sel <- which(is.na(c(X[, colidx])))
        optimize(function(rho){
            sum(invprobit(X[, colidx], rho, zeta[colidx], log=TRUE, oneminus = TRUE), na.rm=TRUE) +
              sum(invprobit(mup[sel, experimental_design[colidx]], rho,
                       zeta[colidx] * sqrt(1 + sigma2p[sel] / zeta[colidx]^2), log=TRUE))
        }, lower=-100, upper=100, maximum=TRUE)$maximum
    })



    # Variance prior
    nu <- 3
    eta <- sigma20 + 0.1

    last_round_params <- list(eta,nu,mu0,sigma20,rho,zeta)
    converged <- FALSE
    iter <- 1

    while(! converged && iter < max_iter){
        if(verbose) message(paste0("Starting iter ", iter))

        # Estimate dropout sigmoid
        sigmoid_est <- if(dropout_curve_calc == "global") {
            fit_global_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }else if(dropout_curve_calc == "sample"){
            fit_sample_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }else if(dropout_curve_calc == "global_scale"){
            fit_global_scale_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }
        rho <- sigmoid_est$rho
        zeta <- sigmoid_est$zeta

        # Make a good estimate for each mean
        sigma2p <- fit_feature_variances(X, mup, rho, zeta, nu, eta, experimental_design)
        mu_vars <- fit_feature_mean_uncertainties(X, rho, zeta, nu, eta, mu0,
                                                    sigma20, experimental_design)
        mup <- fit_feature_means(X, sigma2p, mu_vars, rho, zeta, nu, eta, mu0,
                                 sigma20, experimental_design)

        # Fit the variance prior
        var_est <- fit_variance_prior(X, rho, zeta, experimental_design)
        nu <- var_est$df_prior
        eta <- var_est$var_prior

        # Fit the location prior
        loc_est <- fit_location_prior(X, mup, zeta, rho, experimental_design)
        mu0 <- loc_est$mu0
        sigma20 <- loc_est$sigma20


        if(verbose){
            message(paste0("eta0 estimate: ", round(eta, 3)))
            message(paste0("nu0 estimate: ", round(nu, 3)))
            message(paste0("mu0 estimate: ", round(mu0, 3)))
            message(paste0("sigma20 estimate: ", round(sigma20, 3)))
            message(paste0("rho estimate: ", paste0(round(rho, 3), collapse = ",")))
            message(paste0("zeta estimate: ", paste0(round(zeta, 3), collapse = ",")))
        }
        error <- sum(mapply(function(new, old){sum(new - old)/length(new)},
                            list(eta,nu,mu0,sigma20,rho,zeta), last_round_params )^2)
        if(verbose) message(paste0("Error: ", error))
        if(error < epsilon) {
            if(verbose) message("converged!")
            converged <- TRUE
        }
        last_round_params <- list(eta,nu,mu0,sigma20,rho,zeta)
        iter <- iter + 1
    }

    names(last_round_params) <- c("eta", "nu", "mu0", "sigma20", "rho", "zeta")
    names(last_round_params$rho) <- colnames(X)
    names(last_round_params$zeta) <- colnames(X)


    feature_params <- list(
        mup=matrix(NA, nrow=nrow(X_bckp), ncol=N_cond),
        sigma2p=rep(NA, times=nrow(X_bckp)),
        sigma2mup=matrix(NA, nrow=nrow(X_bckp), ncol=N_cond)
    )
    feature_params$mup[sel, ] <- mup
    feature_params$sigma2p[sel] <- sigma2p
    feature_params$sigma2mup[sel, ] <- mu_vars

    colnames(feature_params$mup) <- levels(experimental_design_fct)
    colnames(feature_params$sigma2mup) <- levels(experimental_design_fct)
    rownames(feature_params$mup) <- rownames(X_bckp)
    rownames(feature_params$sigma2mup) <-  rownames(X_bckp)

    antisel <- setdiff(seq_len(nrow(X_bckp)), sel)
    if(length(antisel) > 0){
        if(verbose) message(paste0("Predicting features which were not in the subsample: ",
                                   length(antisel), " elements."))
        add_feats <- predict_feature_parameters(X_bckp[antisel, ], experimental_design_fct,
                                                last_round_params, verbose=verbose)
        feature_params$mup[antisel, ] <- add_feats$mup
        feature_params$sigma2p[antisel] <- add_feats$sigma2p
        feature_params$sigma2mup[antisel, ] <- add_feats$sigma2mup
    }

    ret <- list(hyper_params = last_round_params,
         feature_params = feature_params,
         experimental_design=experimental_design_fct,
         error=error, converged=converged)

    class(ret) <- "prodd_parameters"
    ret
}



setGeneric("fit_parameters")

#' @describeIn fit_parameters S4 method of \code{fit_parameters} for
#'   \code{SummarizedExperiment}
setMethod("fit_parameters",
          c(X = "SummarizedExperiment"),
          function(X, experimental_design,
                   dropout_curve_calc=c("sample", "global_scale", "global"),
                   frac_subsample=1.0, n_subsample=round(nrow(X) * frac_subsample),
                   max_iter=10, epsilon=1e-3, verbose=FALSE){

              # First extract the experimental_design column from the colData
              if(length(experimental_design) == 1){
                  if(! experimental_design %in% colnames(colData(X))) {
                      stop(paste0("'experimental_design' must reference a ",
                                  "column in colData(X). Ie. one of: ",
                                  paste0(colnames(colData(X)), collapse=", ")))
                  }
                  experimental_design <- colData(X)[, experimental_design, drop=TRUE]

              }
              params <- fit_parameters(SummarizedExperiment::assay(X), experimental_design,
                                          dropout_curve_calc, frac_subsample, n_subsample,
                                          max_iter, epsilon, verbose)
              params
          })

#' @describeIn fit_parameters S4 method of \code{fit_parameters} for
#'   \code{MSnSet}
setMethod("fit_parameters",
          c(X = "MSnSet"),
          function(X, experimental_design,
                   dropout_curve_calc=c("sample", "global_scale", "global"),
                   frac_subsample=1.0, n_subsample=round(nrow(X) * frac_subsample),
                   max_iter=10, epsilon=1e-3, verbose=FALSE){

              # First extract the experimental_design column from the pData
              if(length(experimental_design) == 1){
                  if(! experimental_design %in% colnames(pData(X))) {
                      stop(paste0("'experimental_design' must reference a ",
                                  "column in pData(X). Ie. one of: ",
                                  paste0(colnames(pData(X)), collapse=", ")))
                  }
                  experimental_design <- pData(X)[, experimental_design, drop=TRUE]

              }
              params <- fit_parameters(Biobase::exprs(X), experimental_design,
                                          dropout_curve_calc, frac_subsample, n_subsample,
                                          max_iter, epsilon, verbose)
              params
          })




#' Fit a Gaussian location prior over all samples
#'
#' The function takes the regularized feature means. The global mu0 is the
#' mean of the regularized feature means. It then calculates the unregularized
#' feature means right of the global mean (so that missing values are less
#' problematic). The global variance estimate is the variance of those
#' unregularized values to the global mean.
#'
#' @return list with elements mu0 and sigma20
#' @keywords internal
fit_location_prior <- function(X, mup, zeta, rho, experimental_design){
    mu0 <- mean(mup, na.rm=TRUE)
    mu_unreg <- fit_unregularized_feature_means(X, mup, mu0, zeta, rho, experimental_design)
    sigma20 <- sum((mu_unreg[! is.na(mu_unreg)] - mu0)^2/sum(! is.na(mu_unreg)))
    list(mu0=mu0, sigma20=sigma20)
}


#' Fit a inverse Chi-squared variance prior
#'
#' Run maximum likelihood estimate on the density of the F distribution with
#' the unregularized variance estimates.
#' @keywords internal
fit_variance_prior <- function(X, rho, zeta, experimental_design){
    DF_eff <- calc_df_eff(X, experimental_design)
    sigma2_unreg <- fit_unregularized_feature_variances(X, rho, zeta, experimental_design)

    var_est <- optim(par=c(eta=1, nu=1), function(par){
        if(par[1] < 0 || par[2] < 0 ) return(Inf)
        - sum(sapply(seq_len(nrow(X))[DF_eff >= 1 & ! is.na(sigma2_unreg)], function(idx){
            df(sigma2_unreg[idx]/par[1], df1=DF_eff[idx], df2=par[2], log=TRUE) - log(par[1])
        }))
    })
    names(var_est$par) <- NULL
    list(var_prior=var_est$par[1], df_prior=var_est$par[2])
}
