




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
                                                DF_eff=calc_df_eff(X, experimental_design),
                                                upper=1000){

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

