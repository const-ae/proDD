




#' Distance matrix computation respecting nonrandom missing values
#'
#' The method calculates the distances between each row of the matrix X and
#' returns a distance matrix. Because of the missing values no exact distance can
#' be calculated instead realistic values for the missing values are considered
#' and the mean with the corresponding variance is calculated for each distance.
#'
#' If no appropriate replacement values are available and the method is called
#' without any parameter in addition to X, it calls `fit_hyperparameters` and
#' infers them.
#'
#' @param X the transposed (!) data matrix
#' @param mu_mis mean of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X.
#' @param var_mis variance of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X.
#' @param experimental_design if no good replacement values are available
#'   you can provide the experimental_design vector
#' @param hyperparameters a list with four elements (mu0, sigma20, rho and zeta),
#'   that are used to calculate a good replacement value.
#' @seealso stats::dist
#' @export
dist_approx <- function(X, mu_mis=NULL, var_mis=NULL,
                        experimental_design=NULL,
                        hyperparameters=NULL,
                        ...){

    if(is.null(hyperparameters) &&  (is.null(mu_mis) || is.null(var_mis))){
        # If no parameters are provided, just infer them
        if(is.null(experimental_design)){
            experimental_design <- rep(1, nrow(X))
        }
        res <- fit_hyperparameters(t(X), experimental_design, calculate_distance=TRUE, ...)
        return(res$distance)
    }else if(!is.null(hyperparameters) && (is.null(mu_mis) || is.null(var_mis))){
        if(is.null(names(hyperparameters))){
            # Assume they are correctly ordered
            names(hyperparameters) <- c("mu0", "sigma20", "rho", "zeta")
        }
        # Just make a good guess based on the dropout curves
        if(is.null(mu_mis)){
            mu_mis <- mapply(function(.mu0, .sigma20, .rho, .zeta)
                        mean_probdropout(mu=.mu0, sigma2=.sigma20,rho=.rho, zeta=.zeta),
                .mu0=hyperparameters$mu0, .sigma20=hyperparameters$sigma20,
                .rho=hyperparameters$rho, .zeta=hyperparameters$zeta)
        }
        if(is.null(var_mis)){
            var_mis <- mapply(function(.mu0, .sigma20, .rho, .zeta)
                    variance_probdropout(mu=.mu0, sigma2=.sigma20,rho=.rho, zeta=.zeta),
                .mu0=hyperparameters$mu0, .sigma20=hyperparameters$sigma20,
                .rho=hyperparameters$rho, .zeta=hyperparameters$zeta)
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
find_approx_for_missing <- function(X, mup, mu_vars, sigma2p,
                                    rho, zeta, experimental_design){

    mu_mis <- proDD:::mply_dbl(seq_len(nrow(X)), ncol=ncol(X), function(idx){
        vapply(seq_len(ncol(X)), function(sample){
            if(is.na(X[idx, sample])){
                mean_probdropout(mup[idx, experimental_design[sample]],
                         mu_vars[idx, experimental_design[sample]] + sigma2p[idx],
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
                         mu_vars[idx, experimental_design[sample]] + sigma2p[idx],
                         rho=rho[sample], zeta=zeta[sample])
            }else{
                NA
            }
        }, FUN.VALUE=0.0)
    })

    list(mu_mis=mu_mis, var_mis=var_mis)

}

