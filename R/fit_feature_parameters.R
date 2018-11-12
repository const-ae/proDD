




#' Find feature parameters for elements of X with a fixed set of hyperparameters
#'
#' @param X data matrix
#' @param experimental_design a vector that specifies which samples belong
#'   to the same condition.
#' @param hyper_params a list with 5 elements (`eta`, `nu`, `mu0`, `sigma20`,
#'   `rho`, and `zeta`). Alternatively the `hyper_params` can be specified
#'   individually.
#' @param mu0 the global mean around which the row means are drawn
#' @param sigma20 the global variance specifying the spread of means around `mu0`.
#' @param nu degrees of freedom for the the global variance prior.
#' @param eta scale of the global variance prior.
#' @param rho vector specifying the intensity where the chance of a dropout is
#'   50/50. Length is either one or \code{ncol(X)}.
#' @param zeta vector specifying the scale of the dropout curve.
#'   Length is either one or \code{ncol(X)}.
#' @param mup Optional matrix that fixes the mean for each row and condition. Default `NULL`
#' @param sigma2p Optional vector that fixes the variance for each row. Default `NULL`
#' @param sigma2mup Optional matrix that fixes the uncertainty of the mean
#'  for each row and condition. Default `NULL`.
#' @param max_iter the maximum number of iterations. Default: 10
#' @param epsilon the error under which the result is considered converged.
#'   Default: 0.001
#' @param verbose boolean that indicates if verbose output is printed to the console.
#'
#' @return list with three elements
#'   \describe{
#'     \item{mup}{a matrix with size \code{nrow(X) * unique(experimental_design)}
#'       with the means for each feature}
#'     \item{sigma20}{a numeric vector with the variance for each feature}
#'     \item{sigma2mup}{a matrix with size \code{nrow(X) * unique(experimental_design)}
#'       with the uncertainty for each `mup`}
#'   }
#'
predict_feature_parameters <- function(X, experimental_design, hyper_params=NULL,
                                       mu0, sigma20, nu, eta, rho, zeta,
                                       mup=NULL, sigma2p=NULL, sigma2mup=NULL,
                                       max_iter=10, epsilon=1e-3, verbose=FALSE){

    if(! is.null(hyper_params)){
        mu0 <- hyper_params$mu0
        sigma20 <- hyper_params$sigma20
        nu <- hyper_params$nu
        eta <- hyper_params$eta
        rho <- hyper_params$rho
        zeta <- hyper_params$zeta
    }
    if(length(rho) == 1){
        rho <- rep(rho, times=ncol(X))
    }
    if(length(zeta) == 1){
        zeta <- rep(zeta, times=ncol(X))
    }

    fixed_mup <- ! is.null(mup)
    fixed_mu_vars <- ! is.null(sigma2mup)
    mu_vars <- sigma2mup
    fixed_sigma2p <- ! is.null(sigma2p)

    experimental_design_fct <- as.factor(experimental_design)
    experimental_design <- as.numeric(experimental_design_fct)

    N_cond <- length(unique(experimental_design))

    # Initialize the parameters
    if(! fixed_mup){
        mup <-  matrix(mu0, nrow=nrow(X), ncol=N_cond)
    }
    if(! fixed_sigma2p){
        sigma2p <-  rep(sigma20, times=nrow(X))
    }
    if(! fixed_mu_vars){
        mu_vars <- matrix(eta, nrow=nrow(X), ncol=N_cond)
    }


    last_round_params <- list(mup, sigma2p, mu_vars)
    converged <- FALSE
    iter <- 1
    while(! converged && iter < max_iter){
        if(verbose) message(paste0("Starting iter ", iter))
        # Make a good estimate for each feature
        if(! fixed_sigma2p){
            sigma2p <- fit_feature_variances(X, mup, rho, zeta, nu, eta, experimental_design)

        }
        if(! fixed_mu_vars){
            mu_vars <- fit_feature_mean_uncertainties(X, rho, zeta, nu, eta, mu0,
                                                      sigma20, experimental_design)
        }
        if(! fixed_mup){
            mup <- fit_feature_means(X, sigma2p, mu_vars, rho, zeta, nu, eta, mu0,
                                     sigma20, experimental_design)
        }


        # Calculate error
        error <- sum(mapply(function(new, old){sum(new - old)/max(length(new),1)},
                            list(mup, sigma2p, mu_vars), last_round_params )^2)
        if(verbose) message(paste0("Error: ", error))
        if(error < epsilon) {
            if(verbose) message("converged!")
            converged <- TRUE
        }
        last_round_params <- list(mup, sigma2p, mu_vars)
        iter <- iter + 1
    }


    colnames(last_round_params[[1]]) <- levels(experimental_design_fct)
    colnames(last_round_params[[3]]) <- levels(experimental_design_fct)
    rownames(last_round_params[[1]]) <- rownames(X)
    rownames(last_round_params[[3]]) <-  rownames(X)

    names(last_round_params) <- c("mup", "sigma2p", "sigma2mup")
    last_round_params
}







#' Find the most likely variance explaining the values for each row in X
#'
#'
fit_feature_variances <- function(X, mup, rho, zeta, nu, eta, experimental_design,
                                  upper=1e3){
    N_cond <- length(unique(experimental_design))

    vapply(seq_len(nrow(X)), function(idx){

        optimize(f=function(sigma2){
            reg <- extraDistr::dinvchisq(sigma2, nu, eta, log=TRUE)
            reg + sum(vapply(seq_len(N_cond), function(cond){
                # All the missing values for this condition
                sel <- which(experimental_design == cond & is.na(X[idx,]))
                x <- X[idx, which(experimental_design == cond)]
                zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
                if(sum(! is.na(x)) == 0){
                    norm <- 0
                }else{
                    norm <- sum(dnorm(x, mean=mup[idx, cond], sd=sqrt(sigma2), log=TRUE), na.rm=TRUE)
                }

                cumulative <- sum(invprobit(mup[idx, cond], rho[sel], zetastar[sel], log=TRUE))
                norm + cumulative
            }, FUN.VALUE=0.0))
        }, lower=0, upper=upper, maximum = TRUE)$maximum
    }, FUN.VALUE = 0.0)

}


#' Calculate the variance of the mean estimator
#'
fit_feature_mean_uncertainties <- function(X, rho, zeta, nu, eta, mu0, sigma20,
                                           experimental_design){
    N_cond <- length(unique(experimental_design))

    mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            sel <- which(experimental_design == cond & is.na(X[idx,]))
            nobs <- sum(! is.na(x))
            if(nobs == 0){
                sigma2obsreg <- if(nu > 2) nu * eta / (nu-2) else nu * eta / (nu + 2)

                sigma2mu <- sigma20
                mu <- mu0
            }else if(nobs == 1) {
                sigma2obsreg <- if(nu > 2) nu * eta / (nu-2) else nu * eta / (nu + 2)

                sigma2mu <- 1/(1/sigma20 + 1/(sigma2obsreg/nobs))
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu+ mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }else{
                sigma2obs <- var(x, na.rm=TRUE)
                sigma2obsreg <- ((nobs-1) * sigma2obs + nu * eta) / (nobs-1 + nu)

                sigma2mu <- 1/(1/sigma20 + 1/(sigma2obsreg/nobs))
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu + mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }
            zetastar <- zeta * sqrt(1 + sigma2obsreg/zeta^2)
            variance_probdropout(mu, sigma2=sigma2mu, rho=rho[sel], zeta=zetastar[sel])
        }, FUN.VALUE=0.0)
    })
}


#' Estimate the mean for each feature
#'
#'
fit_feature_means <- function(X, sigma2p, sigma2mus, rho, zeta, nu, eta,
                              mu0, sigma20, experimental_design){
    N_cond <- length(unique(experimental_design))

    mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            sel <- which(experimental_design == cond & is.na(X[idx,]))
            nobs <- sum(! is.na(x))
            if(nobs == 0){
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- mu0
            }else if(nobs == 1){
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu + mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }else{
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu+ mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }

            zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
            mean_probdropout(mu, sigma2=sigma2mu, rho=rho[sel], zeta=zetastar[sel])
        }, FUN.VALUE=0.0)
    })
}


#' Calculate the effective degrees of freedom
#'
#' The effective degrees of freedom are the sum of all non NA observation
#' minus the number of non-empty conditions.
#'
calc_df_eff <- function(X, experimental_design){
    N_cond <- length(unique(experimental_design))
    vapply(seq_len(nrow(X)), function(idx){
        df_eff <- vapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            nobs <- sum(! is.na(x))
            if(nobs <= 1) NA else nobs
        }, FUN.VALUE=0.0)
        sum(df_eff, na.rm=TRUE) - sum(!is.na(df_eff))
    }, FUN.VALUE = 0.0)
}


#' Find the variance without using the regularization
#'
#' This method is important for finding the variance prior hyperparameters.
#' The results are similar to the unregularized results, but I cannot reuse
#' those because the prior estimation needs to be unbiased from former rounds
#'
#' The methods returns \code{NA} for each feature with less than 2 observations.
fit_unregularized_feature_variances <- function(X, rho, zeta,
                                                experimental_design,
                                                upper=1000){
    DF_eff <- calc_df_eff(X, experimental_design)
    N_cond <- length(unique(experimental_design))
    vapply(seq_len(nrow(X)), function(idx){
        nobs <- sum(! is.na(X[idx, ]))
        if(nobs <= 1){
            # Not enough observations to do anything
            NA
        }else if(sum(is.na(X[idx, ])) == 0){
            # No missing values, avoid expensive calculations
            res <- vapply(seq_len(N_cond), function(cond){
                x <- X[idx, which(experimental_design == cond)]
                nobs <- sum(! is.na(x))
                if(nobs > 1){
                    c((nobs-1) * var(x, na.rm=TRUE), nobs-1)
                }else{
                    c(NA, NA)
                }
            }, FUN.VALUE=c(0.0, 0.0))
            sum(res[1, ], na.rm=TRUE) / sum(res[2, ], na.rm=TRUE)
        }else{
            mode <- optimize(f=function(.sigma2){
                zetastar <- zeta * sqrt(1 + .sigma2/ zeta^2)

                1/.sigma2 * prod(vapply(seq_len(N_cond), function(cond){
                    x <- X[idx, which(experimental_design == cond)]
                    nobs <- sum(! is.na(x))
                    if(nobs <= 1){
                        NA
                    }else{
                        obsvar <- var(x, na.rm=TRUE)
                        muobs <- mean(x, na.rm=TRUE)
                        v1 <- sqrt(.sigma2)^(-nobs) * exp(-(nobs-1) * obsvar / (2 * .sigma2))
                        # Numerically integrating out the mean, to get a better variance estimate
                        v2 <- tryCatch({integrate(function(mu){
                            exp(-1/(2 * .sigma2) * nobs * (muobs-mu-muobs)^2) * (
                                if(! any(is.na(x))) 1
                                else
                                    apply(msply_dbl(which(experimental_design == cond & is.na(X[idx, ])), function(hidx)
                                        invprobit(mu-muobs, rho[hidx], zetastar[hidx])), 2, prod)
                            )
                        }, lower=-Inf, upper=Inf)$value
                        }, error=function(err){warning("sigma2_adapted_approx:\n");warning(err); 0})
                        v1 * v2
                    }
                }, FUN.VALUE = 0.0), na.rm=TRUE)
            }, lower=0, upper=upper, maximum = TRUE)$maximum
            df_eff <- DF_eff[idx]
            # adapt_s2 <- mode * (df_eff+2)/(df_eff)   # Convert mode into tau^2 of the Inv.Chisq. with same mode
            adapt_s2 <- mode * (df_eff+3)/(df_eff) # Weird hack that seems to work
            adapt_s2
        }
    }, FUN.VALUE=0.0)
}




#' Calculate the unregularized means for all rows for which mup is larger than mu0
#'
#'
fit_unregularized_feature_means <- function(X, mup, mu0, zeta, rho, experimental_design){
    N_cond <- length(unique(experimental_design))
    mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        sapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            sel <- which(experimental_design == cond & is.na(X[idx,]))
            nobs <- sum(! is.na(x))
            if(nobs < 2 || mup[idx, cond] < mu0){
                NA
            }else{
                sigma2 <-  var(x, na.rm=TRUE)
                sigma2mu <- sigma2 / nobs
                mu <- mean(x, na.rm=TRUE)
                zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
                tryCatch(
                    mean_probdropout(mu, sigma2=sigma2mu, rho=rho[sel],
                                     zeta=zetastar[sel]),
                    error=function(err) NA
                )
            }
        })
    })
}



#' Fit a Gaussian location prior over all samples
#'
#' The function takes the regularized feature means. The global mu0 is the
#' mean of the regularized feature means. It then calculates the unregularized
#' feature means right of the global mean (so that missing values are less
#' problematic). The global variance estimate is the variance of those
#' unregularized values to the global mean.
#'
#' @return list with elements mu0 and sigma20
#'
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
#'
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
