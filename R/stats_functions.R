

#' Inverse probit function
#'
#' Calculate the values of the sigmoidal function that is defined by the
#' cumulative normal distribution function (\code{\link{pnorm}}). This
#' method provides a convient wrapper for the \code{pnorm} that automatically
#' handles negative zeta and is more consistent in its naming.
#'
#' @param x numeric vector
#' @param rho numeric vector of length 1 or the same length as x. Specifies
#'  the inflection point of the inverse probit curve.
#' @param zeta numeric vector of length 1 or the same length as x. Specifies
#'  how flat the curve is at the inflection point of the inverse probit curve.
#' @param log boolean if the log of the result is returned
#' @param oneminus boolean if one minus the result is returned
#'
#' @return a numeric vector of \code{length(x)}.
#' @examples
#'  xg <- seq(-5, 5, length.out=101)
#'  plot(xg, invprobit(xg, rho=-2, zeta=-0.3))
#' @export
invprobit <- function(x, rho, zeta, log=FALSE, oneminus=FALSE){
    stopifnot(length(x) == length(rho) || length(rho) == 1 || length(x) == 1)
    stopifnot(length(x) == length(zeta) || length(zeta) == 1 || length(x) == 1)

    sign_sum <- sum(sign(zeta), na.rm=TRUE)
    if(length(zeta) > 1 && abs(sign_sum) < sum(!is.na(zeta), na.rm=TRUE))
        stop("All zeta must have the same sign")
    if(sign_sum < 0){
        pnorm(x, mean=rho, sd=-zeta, lower.tail = oneminus, log.p=log)
    }else{
        pnorm(x, mean=rho, sd=zeta, lower.tail = ! oneminus, log.p=log)
    }
}
