




#' Find the inverse probit parameters that best explain the pattern of missing
#' values.
#'
#' There are three different modes for fitting the inverse probit dropout
#' curves:
#' \describe{
#' \item{fit_global_dropout_curves}{One curve for all samples}
#' \item{fit_sample_dropout_curves}{One curve for each sample}
#' \item{fit_global_scale_dropout_curves}{One curve for each sample,
#'   but all have a common scale}
#' }
#'
fit_global_dropout_curves <- function(){

}

#' @rdname fit_global_dropout_curves
fit_sample_dropout_curves <- function(){

}

#' @rdname fit_global_dropout_curves
fit_global_scale_dropout_curves <- function(X, mup, sigma2p, experimental_design,
                             prev_zeta, prev_rho,
                             max_iter=10, epsilon=1e-7){

    zeta <- mean(prev_zeta)
    rho <- prev_rho

    last_round_params  <- list(zeta, rho)
    iter <- 1
    converged <- FALSE
    while(! converged && iter < max_iter){
        zeta <- nlm(function(zeta){
            if(zeta > 0) return(1e100) # Aka. Inf, but Inf triggers a warning
            -sum(vapply(seq_len(ncol(X)), function(colidx){
                sum(invprobit(X[, colidx], rho[colidx], zeta, log=TRUE, oneminus = TRUE), na.rm=TRUE) +
                    sum(invprobit(mup[which(is.na(c(X[, colidx]))), experimental_design[colidx]],
                         rho[colidx], zeta * sqrt(1 + sigma2p[which(is.na(c(X[, colidx])))] / zeta^2), log=TRUE))
            }, FUN.VALUE = 0.0))

        }, p=zeta, stepmax=1)$estimate

        rho <- sapply(seq_len(ncol(X)), function(colidx){
            nlm(function(rho){
                -(sum(invprobit(X[, colidx], rho, zeta, log=TRUE, oneminus = TRUE), na.rm=TRUE) +
                      sum(invprobit(mup[which(is.na(c(X[, colidx]))), experimental_design[colidx]],
                         rho, zeta * sqrt(1 + sigma2p[which(is.na(c(X[, colidx])))] / zeta^2), log=TRUE)))
            }, p=rho[colidx])$estimate
        })

        iter <- iter + 1
        error <- sum(mapply(function(new, old){
            sum(new - old) / length(new)
        }, list(zeta, rho), last_round_params)^2)
        last_round_params <- list(zeta, rho)
        converged <- error < epsilon
    }

    zeta <- rep(zeta, ncol(X))
    list(zeta=zeta, rho=rho, converged=converged)
}


