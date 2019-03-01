












#' Estimate the mean of each protein and condtion using MCMC sampling.
#'
#'
#' The method uses MCMC sampling implemented by Stan to find a good approximation
#' to the posterior for each mean. The method is called using the parameters inferred
#' from \code{fit_parameters()}. The parameter \code{niter} specifies how many
#' samples are drawn for each mean. The number of samples that are available
#' after burnin are \eqn{samples = niter / 2 * nchains}.
#'
#' @param X the numerical data where each column is one sample and each row
#'   is one protein. Missing values are coded as \code{NA}.
#' @param params an object of class `prodd_parameters`
#' @param niter the number of iteration for each posterior and chain. In
#'   the end you will have \code{niter * nchains / 2} samples. Default: 1000
#' @param nchains the number of chains to run in parallel. Default: 4
#' @param ncores the number of cores that are used in parallel. Deafault: nchains.
#' @param batch_size Often it is faster to run the inference for multiple
#'   proteins at the same time, but there is a limit after which too many proteins
#'   slow down the inference. Set the number of proteins that are considered in
#'   a single run. Default: 1000
#' @param verbose boolean indicating how much output the function generates. D
#'   Default: `TRUE`
#' @param ... additional parameters for passed to \code{rstan::sampling}.
#'
#'
#' @return   a list with one matrix per condtion. Each matrix has one row
#'   per protein and one column per MCMC sample.
#'
#' @export
sample_protein_means <- function(X, params,
                                niter=1000, nchains=4, ncores=nchains,
                                batch_size=1000, verbose=TRUE, ...){
    calculate_posterior(X, params,
                        niter=niter, nchains=nchains, ncores=ncores,
                        batch_size = batch_size, verbose=verbose, ...)
}



setGeneric("sample_protein_means")

#' @describeIn sample_protein_means S4 method of \code{sample_protein_means} for
#'   \code{SummarizedExperiment}
setMethod("sample_protein_means",
          c(X = "SummarizedExperiment"),
          function(X, params,
                   niter=1000, nchains=4, ncores=nchains,
                   batch_size=1000, verbose=TRUE, ...){

              sample_protein_means(SummarizedExperiment::assay(X), params, niter,
                                   nchains, ncores, batch_size, verbose)
          })

#' @describeIn sample_protein_means S4 method of \code{sample_protein_means} for
#'   \code{MSnSet}
setMethod("sample_protein_means",
          c(X = "MSnSet"),
          function(X, params,
                   niter=1000, nchains=4, ncores=nchains,
                   batch_size=1000, verbose=TRUE, ...){
              sample_protein_means(Biobase::exprs(X), params, niter,
                                   nchains, ncores, batch_size, verbose)
          })




#' Internal function that does the actual work for \code{sample_protein_means}
#'
#' @inheritParams sample_protein_means
#' @param mu0,sigma20,nu,eta,rho,zeta,experimental_design a simple way to
#'   to overwrite the hyper-parameters in `params`.
#'
#'
#' @return   a list with one matrix per condtion. Each matrix has one row
#'   per protein and one column per MCMC sample.
#'
#' @seealso \code{\link{sample_protein_means}}
calculate_posterior <- function(X, params=NULL,
                    mu0, sigma20, nu, eta, rho, zeta, experimental_design,
                    niter=1000, nchains=4, ncores=nchains, batch_size=1000, verbose=TRUE,
                    ...){

    if(! is.null(params)){
        if(! is.prodd_parameters(params)){
            stop("params must be an object of class prodd_parameters, which",
                 " is for example returned by fit_parameters()")
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
    dots <- list(...)
    if(is.null(dots$iter)){
        dots$iter <- niter
    }
    if(is.null(dots$open_progress)){
        dots$open_progress <- FALSE
    }
    if(is.null(dots$verbose)){
        dots$verbose <- FALSE
    }
    if(is.null(dots$refresh)){
        dots$refresh <- if(verbose)  max(niter/10, 1) else 0
    }

    experimental_design_fct <- as.factor(experimental_design)
    options(mc.cores = ncores)

    stopifnot(is.matrix(X))
    if(length(rho) != ncol(X)){
        stop("Need one rho for each sample")
    }
    if(length(zeta) != ncol(X)){
        stop("Need one zeta for each sample")
    }
    batch_results <- lapply(seq_len(ceiling(nrow(X)/batch_size)), function(batch_idx){
        .X <- X[((batch_idx-1) * batch_size+1):min(nrow(X), batch_idx*batch_size), , drop=FALSE]

        # Hack to work around the lack of NA support of stan
        .X[is.na(.X)] <- Inf

        args <- c(list(object=stanmodels$batch_skewed_posterior), list(data=list(
            nsamples=ncol(.X),
            nrows=nrow(.X),
            ncond=length(levels(experimental_design_fct)),
            totalmissing=sum(is.infinite(c(.X))),

            X=.X,
            experimental_design=as.numeric(experimental_design_fct),

            mu0=mu0,
            sigma20=sigma20,
            zeta=zeta,
            rho=rho,
            nu=nu,
            eta=eta
        )), dots)

        sam_f <- function(){
            do.call(rstan::sampling, args)
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
    result <- lapply(seq_along(levels(experimental_design_fct)), function(cond){
        mat <- t(do.call(cbind, lapply(batch_results, function(e) e$mu[,,cond])))
        rownames(mat) <- rownames(X)
        mat
    })
    names(result) <- levels(experimental_design_fct)
    result
}
