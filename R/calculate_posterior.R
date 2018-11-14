

#' Calculate the posterior distribution of the mean for each condition
#'
#'
#' The method uses MCMC sampling to find a good approximation to the posterior
#' for each mean. The method is called using the parameters inferred from
#' \code{fit_hyperparameters()}. The parameter \code{niter} specifies how many
#' samples are drawn for each mean. The number of samples that are available
#' after burnin are \eqn{samples = niter / 2 * nchains}
#'
#' @param X the data matrix
#' @param params an object of class `prodd_parameters`
#' @param niter the number of iteration for each posterior and chain. In
#'   the end you will have \code{niter * nchains / 2} samples. Default: 1000
#' @param nchains the number of chains to run in parallel. Default: 4
#' @param batch_size Often it is faster to run the inference for multiple
#'   features at the same time, but there is a limit after which too many features
#'   slow down the inference. Set the number of features that are considered in
#'   a single run. Default: 1000
#' @param verbose boolean indicating how much output the function generates. D
#'   Default: `TRUE`
#' @export
calculate_posterior <- function(X, params=NULL,
                    mu0, sigma20, nu, eta, rho, zeta, experimental_design,
                    niter=1000, nchains=4, ncores=nchains, batch_size=1000, verbose=TRUE){

    if(! is.null(params)){
        if(! is.prodd_parameters(params)){
            stop("params must be an object of class prodd_parameters, which",
                 " is for example returned by fit_hyperparameters()")
        }
        if(missing(experimental_design) || is.null(experimental_design)){
            experimental_design <- params$experimental_design
        }
        mu0 <- params$hyper_params$mu0
        sigma20 <- params$hyper_params$sigma2
        nu <- params$hyper_params$nu
        eta<- params$hyper_params$eta
        rho <- params$hyper_params$rho
        zeta <- params$hyper_params$zeta

    }

    options(mc.cores = ncores)

    stopifnot(is.matrix(X))
    if(length(rho) != ncol(X)){
        stop("Need one rho for each sample")
    }
    if(length(zeta) != ncol(X)){
        stop("Need one zeta for each sample")
    }
    result <- lapply(unique(experimental_design), function(cond){
        cond_posts <- lapply(seq_len(ceiling(nrow(X)/batch_size)), function(batch_idx){
            .X <- X[((batch_idx-1) * batch_size+1):min(nrow(X), batch_idx*batch_size),
                    which(experimental_design == cond), drop=FALSE]

            # Hack to work around the lack of NA support of stan
            .X[is.na(.X)] <- Inf


            sam_f <- function(){
                rstan::sampling(stanmodels$batch_skewed_posterior, data=list(
                  nsamples=ncol(.X),
                  nrows=nrow(.X),
                  totalmissing=sum(is.infinite(c(.X))),

                  X=.X,

                  mu0=mu0,
                  sigma20=sigma20,
                  zeta=zeta[which(experimental_design == cond)],
                  rho=rho[which(experimental_design == cond)],
                  nu=nu,
                  eta=eta
                ), iter=niter, open_progress=FALSE, verbose=FALSE,
                refresh=if(verbose)  max(niter/10, 1) else 0)
            }

             if(verbose){
                 sam <- sam_f()
            }else{
                sink(tempfile())
                sam <- sam_f()
                sink()
            }

            list(mu=rstan::extract(sam, "mu")$mu)
        })
        t(do.call(cbind, lapply(cond_posts, function(e) e$mu)))
    })
    names(result) <- unique(experimental_design)
    result
}
