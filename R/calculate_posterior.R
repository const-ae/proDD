

#' Calculate the posterior distribution of the mean for each condition
#'
#' The method uses MCMC sampling to find a good approximation to the posterior
#' for each mean. The method is called using the parameters inferred from
#' \code{fit_hyperparameters()}. The parameter \code{niter} specifies how many
#' samples are drawn for each mean. The number of samples that are available
#' after burnin are \eqn{samples = niter / 2 * nchains}
#'
#' @export
calculate_posterior <- function(X, mu0, sigma20, nu, eta, rho, zeta, experimental_design,
                                niter=1000, nchains=4, ncores=nchains, batch_size=1000, verbose=TRUE){
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
