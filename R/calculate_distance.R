




#' Distance matrix computation respecting nonrandom missing values
#'
#' The method calculates the distances between each row of the matrix X and
#' returns a distance matrix. Because of the missing values no exact distance can
#' be calculated instead realistic values for the missing values are considered
#' and the mean with the corresponding variance is calculated for each distance.
#'
#' If no appropriate replacement values are available you can call it with the
#' `params` argument, which must be an object of type `prodd_parameters`. Usually
#' the object that is returned by \code{fit_hyperparameters()} is used.
#'
#' Alternatively, if more flexibility is necessary you can also call the method using
#' `feature_params` which provide the mean and variance for each feature and
#' condition. Combined with `rho` and `zeta` they can be used to determine
#' the replacement values. If also no `feature_params` are available, you
#' can also just call the function with `hyper_params` and all of the above
#' values are calculated. In both cases it is necessary to provide the
#' `experimental_design` parameter so there is a clear mapping between
#' sample and condition. In case mu_mis and var_mis are inferred, the
#' \code{length(experimental_design)} is used to determine if `X` needs to be
#' transposed before inferring the parameters.
#'
#' @param X the data matrix
#' @param mu_mis mean of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X.
#' @param var_mis variance of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X.
#' @param params an object of class `prodd_parameters` which is returned by the
#'   \code{fit_hyperparameters()} function.
#' @param feature_params a list with three elements (mup, sigma2mup and sigma2p)
#' @param rho the dropout curve positions. Necessary if you call the function with
#'   `feature_params`.
#' @param zeta the dropout curve width Necessary if you call the function with
#'   `feature_params`.
#' @param hyper_params a list with four elements (mu0, sigma20, rho and zeta),
#'   that are used to calculate a good replacement value.
#' @param experimental_design assignment to which condition each sample belongs.
#' @seealso stats::dist
#' @export
dist_approx <- function(X, mu_mis=NULL, var_mis=NULL,
                        params=NULL,
                        feature_params=NULL, rho=NULL, zeta=NULL,
                        hyper_params=NULL,
                        experimental_design=NULL,
                        ...){



    if((is.null(mu_mis) || is.null(var_mis))){
        if(! is.null(params)){
            if(! is.prodd_parameters(x)){
                stop("params must be an object of class prodd_parameters, which",
                     " is for example returned by fit_hyperparameters()")
            }
            if(is.null(experimental_design)){
                experimental_design <- params$experimental_design
            }
            # Check if X needs to be rotated
            if(length(experimental_design) == ncol(X)){
                mis <- find_approx_for_missing(X, params, experimental_design=experimental_design)
                mu_mis <- mis$mu_mis
                var_mis <- mis$var_mis
            }else if(length(experimental_design) == nrow(X)){
                mis <- find_approx_for_missing(t(X), params, experimental_design=experimental_design)
                mu_mis <- t(mis$mu_mis)
                var_mis <- t(mis$var_mis)
            }else{
                stop("Length of experimental_design doesn't match ncol(X) nor nrow(X)")
            }
        }else{
            stopifnot(! is.null(experimental_design))
            if(is.null(rho)){
                rho <- hyper_params$rho
            }
            if(is.null(zeta)){
                zeta <- hyper_params$zeta
            }
            if(! is.null(feature_params)){

                # Check if X needs to be rotated
                if(length(experimental_design) == ncol(X)){
                    mis <- find_approx_for_missing(X, mup=feature_params$mup, sigma2mup=feature_params$sigma2mup,
                                           sigma2p=feature_params$sigma2p, rho=rho, zeta=zeta,
                                           experimental_design=experimental_design)
                    mu_mis <- mis$mu_mis
                    var_mis <- mis$var_mis
                }else if(length(experimental_design) == nrow(X)){
                    mis <- find_approx_for_missing(t(X), mup=feature_params$mup, sigma2mup=feature_params$sigma2mup,
                                           sigma2p=feature_params$sigma2p, rho=rho, zeta=zeta,
                                           experimental_design=experimental_design)
                    mu_mis <- t(mis$mu_mis)
                    var_mis <- t(mis$var_mis)
                }else{
                    stop("Length of experimental_design doesn't match ncol(X) nor nrow(X)")
                }



            }else if(!is.null(hyper_params)){
                if(is.null(names(hyper_params))){
                    # Assume they are correctly ordered
                    names(hyper_params) <- c("eta", "nu", "mu0", "sigma20", "rho", "zeta")
                }
                # Check if X needs to be rotated
                if(length(experimental_design) == ncol(X)){
                    feature_params <- predict_feature_parameters(X, experimental_design, hyper_params)
                    mis <- find_approx_for_missing(X, mup=feature_params$mup, sigma2mup=feature_params$sigma2mup,
                                           sigma2p=feature_params$sigma2p, rho=rho, zeta=zeta,
                                           experimental_design=experimental_design)
                    mu_mis <- mis$mu_mis
                    var_mis <- mis$var_mis
                }else if(length(experimental_design) == nrow(X)){
                    feature_params <- predict_feature_parameters(t(X), experimental_design, hyper_params)
                    mis <- find_approx_for_missing(t(X),mup=feature_params$mup, sigma2mup=feature_params$sigma2mup,
                                           sigma2p=feature_params$sigma2p, rho=rho, zeta=zeta,
                                           experimental_design=experimental_design)
                    mu_mis <- t(mis$mu_mis)
                    var_mis <- t(mis$var_mis)
                }else{
                    stop("Length of experimental_design doesn't match ncol(X) nor nrow(X)")
                }


            }else {
                stop(paste0("Either need (mu_mis and var_mis) or (feature_params and rho",
                            " and zeta) or hyper_params"))
            }
        }
    }

    if(length(mu_mis) == 1){
        # One replacement for every missing value
        mu_mis <- matrix(mu_mis, nrow=nrow(X), ncol=ncol(X))
    }else if(length(mu_mis) == nrow(X)){
        # One replacement for every sample
        mu_mis <- matrix(rep(mu_mis, times=ncol(X)), nrow=nrow(X), ncol=ncol(X))
    }else if(! is.matrix(mu_mis) || nrow(mu_mis) != nrow(X) || ncol(mu_mis) != ncol(X)){
        stop("Dimensions of mu_mis and X don't match")
    }

    if(length(var_mis) == 1){
        # One replacement for every missing value
        var_mis <- matrix(var_mis, nrow=nrow(X), ncol=ncol(X))
    }else if(length(var_mis) == nrow(X)){
        # One replacement for every sample
        var_mis <- matrix(rep(var_mis, times=ncol(X)), nrow=nrow(X), ncol=ncol(X))
    }else if(! is.matrix(var_mis) || nrow(var_mis) != nrow(X) || ncol(var_mis) != ncol(X)){
        stop("Dimensions of var_mis and X don't match")
    }

    mu_res <- matrix(NA, nrow=nrow(X), ncol=nrow(X),
                     dimnames = list(rownames(X), rownames(X)))
    var_res <- matrix(NA, nrow=nrow(X), ncol=nrow(X),
                      dimnames = list(rownames(X), rownames(X)))


    for(idx in seq_len(nrow(X)-1)){
        for(idx2 in idx:nrow(X)){
            x <- X[idx, ]
            mu1 <- x
            mu1[is.na(x)] <- mu_mis[idx, is.na(x)]
            sigma1 <- rep(0, length(x))
            sigma1[is.na(x)] <- var_mis[idx, is.na(x)]

            y <- X[idx2, ]
            mu2 <- y
            mu2[is.na(y)] <- mu_mis[idx2, is.na(y)]
            sigma2 <- rep(0, length(y))
            sigma2[is.na(y)] <- var_mis[idx2, is.na(y)]


            ds <- distance_sq(mu1, sigma1, mu2, sigma2)
            mu_res[idx2, idx] <- ds$mean
            var_res[idx2, idx] <- ds$var
        }
    }


    # The sqrt can be linearly approximated at the mean to update the variance
    list(mean=as.dist(sqrt(mu_res)), var=as.dist(var_res * 1/(4 * mu_res)))
}



#' Square distance between two Gaussian distributions
#'
#' The function takes the mean and the diagonal of the covariance matrix
#' as vector and calculates the mean and variance of their distance distribution.
#' The formulas are based on [1] page 53.
#'
#' @return a list with elements `mean` and `var`
#'
#' 1. Mathai, A. & Provost, S. Quadratic Forms in Random Variables. (1992).
#'
distance_sq <- function(mu1, sigma1, mu2, sigma2){
    mu <- mu2 - mu1
    sigma <- sigma1 + sigma2

    analyt_mean <- sum(sigma) + sum(mu^2)
    analyt_var <- 2 *  sum(sigma^2) + 4 * sum(mu^2 * sigma)
    list(mean=analyt_mean, var=analyt_var)
}


#' Find a Gaussian approximation for each missing value in X
#'
#' @return a list with with two elements `mu_mis` and `var_mis`
#'   both of which are matrices with the same dimensions as X
#'   and values wherever X has a missing value.
#'
find_approx_for_missing <- function(X, params=NULL, experimental_design, mup, sigma2mup, sigma2p,
                                    rho, zeta){
    if(! is.null(params)){
        if(! is.prodd_parameters(x)){
            stop("params must be an object of class prodd_parameters, which",
                 " is for example returned by fit_hyperparameters()")
        }
        if(is.null(experimental_design)){
            experimental_design <- params$experimental_design
        }
        mup <- params$feature_params$mup
        sigma2mup <- params$feature_params$sigma2mup
        sigma2p <- params$feature_params$sigma2p
        rho <- params$hyper_params$rho
        zeta <- params$hyper_params$zeta
    }

    mu_mis <- proDD:::mply_dbl(seq_len(nrow(X)), ncol=ncol(X), function(idx){
        vapply(seq_len(ncol(X)), function(sample){
            if(is.na(X[idx, sample])){
                mean_probdropout(mup[idx, experimental_design[sample]],
                                 sigma2mup[idx, experimental_design[sample]] + sigma2p[idx],
                         rho=rho[sample], zeta=zeta[sample])
            }else{
                NA
            }
        }, FUN.VALUE=0.0)
    })

    var_mis <- proDD:::mply_dbl(seq_len(nrow(X)), ncol=ncol(X), function(idx){
        vapply(seq_len(ncol(X)), function(sample){
            if(is.na(X[idx, sample])){
                variance_probdropout(mup[idx, experimental_design[sample]],
                         sigma2mup[idx, experimental_design[sample]] + sigma2p[idx],
                         rho=rho[sample], zeta=zeta[sample])
            }else{
                NA
            }
        }, FUN.VALUE=0.0)
    })

    list(mu_mis=mu_mis, var_mis=var_mis)

}

