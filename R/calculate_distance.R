




#' Distance matrix computation respecting non-random missing values
#'
#' The method calculates the distances between each samples or proteins of the matrix X and
#' returns a distance matrix and the corresponding uncertainty. Because of the missing
#' values no exact distance can be calculated instead realistic values for the missing
#' values are considered and the mean with the corresponding variance is
#' calculated for each distance.
#'
#' Usually the method is called with the data matrix `X `and
#' the object that is returned by \code{fit_parameters()} or the result
#' obtained by first calling \code{transform_parameters()} to remove the group
#' information. The `params` object must be  of type `prodd_parameters`.
#'
#' If particular information are available, where the missing values would have
#' been (ie. mean and variance for each missing value), they can instead of the
#' params object be provided in form of two matrices (with the same dimensions as X)
#' or individual values.
#'
#'
#' Unlike the \code{stats::dist} function which always calculates the distance
#' between the rows of the matrix and one transposes X to find the distances
#' between columns, this method uses the \code{by_sample} parameter. In `X` (and
#' correspondingly for `mu_mis` and `var_mis`)
#' the columns always correspond to the samples and the rows to the proteins. By default
#' the distances are calculated between the samples, to calculate the distances
#' between the proteins set \code{by_sample=FALSE}.
#'
#' @param X the numerical data where each column is one sample and each row
#'   is one protein. Missing values are coded as \code{NA}.
#' @param params an object of class `prodd_parameters` which for example is returned
#'    by the \code{fit_parameters()} function.
#' @param by_sample boolean. Indicate if the distances between samples (columns) or
#'   proteins (rows) is calculated. Default: TRUE.
#' @param blind boolean. If one provides the params argument infered by the
#'   \code{find_parameters()} function, the feature parameters contain information
#'   about the condition of each sample. This can be undesirable if one wants to
#'   infer unsupervised sample similarities for quality control. This is the most
#'   common use case for the \code{dist_approx()} function. Thus the function by default
#'   removes condition information to give unbiased distance estimates by internally
#'   transforming the params object using
#'   \code{transform_parameters(params, rep(1, length(params$experimental_design)))}.
#'   Default: TRUE.
#' @param mu_mis mean of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X. Can be provided instead of the `params` parameters.
#' @param var_mis variance of the replacement values. Can be a single number, a
#'   vector with one number for each sample or a matrix with the same dimensions
#'   as X. Can be provided instead of the `params` parameters.
#' @seealso \code{\link[stats]{dist}}
#'
#' @return a list with two elements:
#'   \describe{
#'     \item{mean}{a distance matrix with the mean of the distance estimate}
#'     \item{var}{the corresponding uncertainty to each distance estimate}
#'   }
#'
#' @export
dist_approx <- function(X, params=NULL, by_sample=TRUE, blind=TRUE,
                        mu_mis=NULL, var_mis=NULL){


    if((is.null(mu_mis) || is.null(var_mis))){
        if(is.null(params) || ! is.prodd_parameters(params)){
            stop("params must be an object of class prodd_parameters, which",
                 " is for example returned by fit_parameters()")
        }else{
            if(blind){
                params <- transform_parameters(params, rep(1, length(params$experimental_design)))
            }

            # Find approximation for the missing values
            mis <- find_approx_for_missing(X, params)
            mu_mis <- mis$mu_mis
            var_mis <- mis$var_mis
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

    if(by_sample){
        X <- t(X)
        mu_mis <- t(mu_mis)
        var_mis <- t(var_mis)
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
#' @keywords internal
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
#' @keywords internal
find_approx_for_missing <- function(X, params=NULL, experimental_design, mup, sigma2mup, sigma2p,
                                    rho, zeta){
    if(! is.null(params)){
        if(! is.prodd_parameters(params)){
            stop("params must be an object of class prodd_parameters, which",
                 " is for example returned by fit_parameters()")
        }
        if(missing(experimental_design) || is.null(experimental_design)){
            experimental_design <- params$experimental_design
        }
        if(ncol(X) != length(experimental_design)){
            if(nrow(X) == length(experimental_design)){
                stop("Number of columns in X doesn't match length of ",
                     "experimental_design, but it looks like you provided ",
                     "the transpose of X. Please instead use the `by_sample` ",
                     "parameter to control if distances are calculated between ",
                     "samples or proteins.")
            }else{
                stop("Number of columns in X doesn't match length of ",
                     "experimental_design")
            }
        }
        stopifnot(ncol(X) == length(experimental_design))
        mup <- params$feature_params$mup
        sigma2mup <- params$feature_params$sigma2mup
        sigma2p <- params$feature_params$sigma2p
        rho <- params$hyper_params$rho
        zeta <- params$hyper_params$zeta
    }

    mu_mis <- mply_dbl(seq_len(nrow(X)), ncol=ncol(X), function(idx){
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

    var_mis <- mply_dbl(seq_len(nrow(X)), ncol=ncol(X), function(idx){
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

    colnames(mu_mis) <- colnames(X)
    colnames(var_mis) <- colnames(X)
    rownames(mu_mis) <- rownames(X)
    rownames(var_mis) <- rownames(X)

    list(mu_mis=mu_mis, var_mis=var_mis)

}

