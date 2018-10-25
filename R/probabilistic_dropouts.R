
#' Density of the probabilistic dropout distribution
#'
#' Calculate the density under a censoring mechanism that probabilistically
#' causes dropouts described with \code{rho} and \code{zeta}. The data is drawn
#' from a normal with \code{mu} and \code{sigma2 * nobs}, but each value drops out with
#' probability according to a sigmoid with its center at \code{rho} and width
#' \code{zeta}:
#' \deqn{
#'   p(x | \mu, \sigma^2, \rho, \zeta) \propto f_{\text{Normal}}(x;\mu, \sigma^2)
#'        \prod_i f_{\Phi}(x; \rho_i, \zeta_i)
#' }{
#'   p(x | \mu, \sigma2, \rho, \zeta) ~ Normal(x; \mu, \sigma2) Prod_i \Phi(x; \rho_i, \zeta_i)
#' }
#' The distribution is related to the extended skewed normal distribution
#' and formally belongs to the class of closed skewed normals.
#'
#' @param x vector of input positions.
#' @param mu double. The mean of the observed values.
#' @param sigma2 double. The variance of the mu estimate.
#' @param rho vector. The positions of the inflection points of the dropout
#'    sigmoids for each sample. Can also be a single number that is repeated
#'    \code{nmis} times. Defaults to an empty vector.
#' @param zeta vector.  The flattness of the dropout sigmoids for each sample.
#'    Can also be a single number that is repeated \code{nmis} times.
#'    Defaults to an empty vector.
#' @param nmis integer The number of missing values.
#'    Defaults to \code{length(rho)}.
#' @param log boolean if the log of the density is returned.
#'
#' @examples
#'   xg <- seq(-5, 5, length.out=101)
#'   plot(xg, dprobdropout(xg, mu=0, sigma2=3, rho=0, zeta=1))
#'
#' @references
#'   1. Azzalini, A. & Capitanio, A. The Skew-Normal and Related Families.
#'     (Cambridge University Press, 2013). doi:10.1017/CBO97811392488911.
#'
#'   2. González-Farías, G., Domínguez-Molina, A. & Gupta, A. K.
#'     Additive properties of skew normal random vectors.
#'     J. Stat. Plan. Inference 126, 521–534 (2004).
#'     doi:10.1016/j.jspi.2003.09.008
dprobdropout <- function(x, mu, sigma2,
                         rho=numeric(0),
                         zeta=numeric(0),
                         nmis=length(rho),
                         log=FALSE){
    if(length(mu) != 1 || length(sigma2) != 1){
        stop("This function is only vectorized for x")
    }
    if(length(rho) != length(zeta)){
        stop("Assuming one sigmoid per sample means",
             "that rho and zeta should be equal length")
    }
    if(abs(sum(sign(zeta))) != length(zeta)){
        stop("The sign of all zetas has to be the same!")
    }

    if(length(rho) != 1 && nmis != length(rho)){
        stop("nmis must be the length of rho/zeta or rho/zeta must",
             " be length one.")
    }else if(length(rho) == 1 && nmis > 1){
        rho <- rep(rho, times=nmis)
        zeta <- rep(zeta, times=nmis)
    }

    if(nmis == 0){
        # Easy way out
        return(dnorm(x, mean=mu, sd=sqrt(sigma2), log=log))
    }

    norm <- dnorm(x, mean=mu, sd=sqrt(sigma2), log=log)
    cumulative <- mply_dbl(x, function(.x) {
        invprobit(.x, rho, zeta, log=log)
    }, ncol=nmis)


    if(all(zeta >= 0)){
        normalizer <- mvtnorm::pmvnorm(lower=rep(-Inf, times=nmis), upper=rep(0, times=nmis), mean=rho - mu,
                                       sigma=zeta^2 * diag(nmis) + matrix(sigma2, nrow=nmis, ncol=nmis))
    }else{
        normalizer <- mvtnorm::pmvnorm(lower=rep(0, times=nmis), upper=rep(Inf, times=nmis), mean=rho - mu,
                                       sigma=zeta^2 * diag(nmis) + matrix(sigma2, nrow=nmis, ncol=nmis))
    }

    if(log){
        as.numeric(norm + apply(cumulative, 1, sum) - log(normalizer))
    }else{
        as.numeric(norm * apply(cumulative, 1, prod) / normalizer)
    }
}
