
#' Density of the probabilistic dropout distribution
#'
#' Calculate the density under a censoring mechanism that probabilistically
#' causes dropouts described with \code{rho} and \code{zeta}. The data is drawn
#' from a normal with \code{mu} and \code{sigma2 * nobs}, but each value drops out with
#' probability according to a inverse probit with its center at \code{rho} and width
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
#'    curves for each sample. Can also be a single number that is repeated
#'    \code{nmis} times. Defaults to an empty vector.
#' @param zeta vector.  The scale of the dropout curves for each sample.
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
#' @export
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
    # Avoid returning nonesense
    if(normalizer == 0){
        mode <- mode_probdropout(mu, sigma2, rho, zeta)
        normalizer <- integrate(function(x){
            norm <- dnorm(x + mode, mean=mu, sd=sqrt(sigma2))
            cumulative <- mply_dbl(x + mode, function(.x) {
                invprobit(.x, rho, zeta)
            }, ncol=nmis)
            norm * apply(cumulative, 1, prod)
        }, lower=-Inf, upper=Inf)$value
    }
    if(log){
        as.numeric(norm + apply(cumulative, 1, sum) - log(normalizer))
    }else{
        as.numeric(norm * apply(cumulative, 1, prod) / normalizer)
    }
}



#' Find the mode of the probabilistic dropout distribution function
#'
#' @inheritParams dprobdropout
#'
#' @section Warning:
#' This function is not vectorized
#'
mode_probdropout <- function(mu, sigma2, rho, zeta){

    if(zeta[1] < 0){
        # This is a dirty hack
        lower <- mu - (abs(mean(rho)-mu)+1) * 50
        upper <- mu
    }else{
        lower <- mu
        upper <- mu + (abs(mean(rho)-mu)+1) * 50
    }

    mode <- uniroot(f=function(x){
        (x-mu)/sigma2 - sign(zeta[1]) *
            sum(exp(dnorm(x, mean=rho, sd=abs(zeta), log=TRUE) -
                        invprobit(x, rho, zeta, log=TRUE)))
    }, lower=lower, upper=upper)$root
    mode
}








#' Find the mean of  of the probabilistic dropout distribution function
#'
#' @inheritParams dprobdropout
#' @param approx boolean. If TRUE match the probabilistic dropout
#'   distribution to a skewed normal distribution and use its mean as an
#'   approximation. If FALSE calculate the first moment using numerical
#'   integration
#'
#' @section Warning:
#' This function is not vectorized
mean_probdropout <- function(mu, sigma2, rho, zeta, log=FALSE, approx=TRUE){

    nmis <- length(rho)
    if(nmis == 0){
        # Easy way out
        return(mu)
    }

    if(approx){
        sn_params <- match_probdropout_skewed_normal(mu, sigma2, rho, zeta)
        # Use the analytical formula to calculate the mean of the skewed normal
        # which approximates the mean of the probabilistic dropout distribution
        ret <- sn_params[1] + sqrt(2/pi) * sn_params[2] *
            sn_params[3] / sqrt(1 + sn_params[3]^2)
        names(ret) <- NULL
        ret
    }else{
        # Calculate the first moment with numerical integration
        integrate(function(x){
            x * dprobdropout(x, mu=0, sigma2=sigma2, rho=rho-mu, zeta=zeta)
        }, lower=-Inf, upper=Inf)$value + mu
    }
}





#' Find the variance of  of the probabilistic dropout distribution function
#'
#'
#' @inheritParams dprobdropout
#' @param approx boolean. If TRUE match the probabilistic dropout
#'   distribution to a skewed normal distribution and use its variance as an
#'   approximation. If FALSE calculate the second moment using numerical
#'   integration
#'
#' @section Warning:
#' This function is not vectorized
variance_probdropout <- function(mu, sigma2, rho, zeta, log=FALSE, approx=TRUE){

    nmis <- length(rho)
    if(nmis == 0){
        # Easy way out
        return(sigma2)
    }

    if(approx){
        sn_params <- match_probdropout_skewed_normal(mu, sigma2, rho, zeta)
        # Use the analytical formula to calculate the variance of the skewed normal
        # which approximates the variance of the probabilistic dropout distribution
        ret <- sn_params[2]^2 *
            (1 - 2 * (sn_params[3] / sqrt(1 + sn_params[3]^2))^2 / pi)
        names(ret) <- NULL
        ret
    }else{
        mom2 <- integrate(function(x){
            x^2 * dprobdropout(x, mu=0, sigma2=sigma2, rho=rho-mu, zeta=zeta)
        }, lower=-Inf, upper=Inf)$value

        mom2 - mean_probdropout(mu=0, sigma2=sigma2, rho=rho-mu, zeta=zeta, approx=FALSE)^2
    }

}


#' Match a skewed normal distribution to the probabilistic dropout density
#'
#' Find the mode of probdropout density and take 5 points around the mode
#' with a distance according to the Laplace approximation variance and calculate
#' the parameters of the skewed normal that most closely matches those
#' 5 points.
#'
#' @inheritParams dprobdropout
match_probdropout_skewed_normal <- function(mu, sigma2, rho, zeta){
    # Laplace Approximation
    mode <- mode_probdropout(mu, sigma2, rho, zeta)
    variance <- -1 / (-1/sigma2 +
                          sign(zeta[1]) * (sum((dnorm(mode, mean=rho, sd=abs(zeta)) /
                                                 invprobit(mode, rho, zeta))^2)) -
                          sign(zeta[1]) * (sum((mode-rho)/abs(zeta)^2 *
                                                dnorm(mode, mean=rho, sd=abs(zeta)) /
                                                invprobit(mode, rho, zeta))) )
    # 5 Point matching
    points <-  c(0,1,-1,0.5,-0.5)
    f <- dprobdropout(mode + points * sqrt(variance), mu, sigma2, rho=rho, zeta=zeta)
    res <- optim(par=c(location=mode, omega=1, alpha=-1), function(par){
        location <- par[1]
        omega <- par[2]
        alpha <- par[3]
        if(omega < 0) return(Inf)
        approx_f <- sn::dsn(mode+ points * sqrt(variance), xi=location,
                            omega=omega, alpha=alpha, tau=0)
        sum((f - approx_f)^2)
    })
    res$par
}

