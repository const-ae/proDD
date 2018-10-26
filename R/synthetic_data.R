


#' Generate a data set according to the probabilistic dropout model
#'
#' @param n_rows integer. The number of rows in the new dataset
#' @param n_replicates integer or vector. The number of replicates in each
#'   condition
#' @param n_conditions
#'
#' @export
generate_synthetic_data <- function(n_rows, n_replicates,
            n_conditions=length(n_replicates),
            frac_changed=0.1,
            n_changed=round(n_rows * min(1, frac_changed)),
            experimental_design=c(unlist(mapply(rep, seq_len(n_conditions), n_replicates))),
            mu0=20, sigma20=10, nu=3, eta=0.3,
            rho=rep(18, times=if(length(n_replicates) == 1) n_replicates*n_conditions else sum(n_replicates)),
            zeta=rep(-1,times=if(length(n_replicates) == 1) n_replicates*n_conditions else sum(n_replicates))){

    stopifnot(all(n_replicates > 0))
    stopifnot(length(n_rows) == 1)
    stopifnot(length(n_changed) == 1)
    stopifnot(length(mu0) == 1)
    stopifnot(length(sigma20) == 1)
    stopifnot(length(nu) == 1)
    stopifnot(length(eta) == 1)

    if(length(n_replicates) == 1){
        n_replicates <- rep(n_replicates, n_conditions)
    }else if(length(n_replicates) != n_conditions){
        stop("The number of elements in n_replicates must match n_conditons")
    }
    stopifnot(length(experimental_design) == sum(n_replicates))
    stopifnot(length(unique(experimental_design)) == n_conditions)

    if(length(rho) == 1){
        rho <- rep(rho, times = sum(n_replicates))
    }else if(length(rho) != sum(n_replicates)){
        stop("length(rho)=", length(rho)," must either be one or match the",
             "  number of samples (sum(n_replicates)=", sum(n_replicates), ")")
    }
    if(length(zeta) == 1){
        zeta <- rep(zeta, times = sum(n_replicates))
    }else if(length(zeta) != sum(n_replicates)){
        stop("length(zeta)=", length(zeta)," must either be one or match the",
             " number of samples (sum(n_replicates)=", sum(n_replicates), ")")
    }


    sigmas2 <- pmin(1/(rchisq(n_rows, df=nu)) * nu * eta, 10)
    mus <- matrix(rep(rnorm(n_rows-n_changed, mu0, sd=sqrt(sigma20)), times=n_conditions),
                  ncol=n_conditions)
    mus <- rbind(mus, t(mply_dbl(seq_len(n_conditions), ncol=n_changed, function(cond){
        rnorm(n_changed, mu0, sd=sqrt(sigma20))
    })))

    t_X <- do.call(cbind, lapply(seq_len(n_conditions), function(cond){
        mply_dbl(seq_len(n_rows), ncol=sum(experimental_design == cond), function(idx){
            rnorm(sum(experimental_design == cond), mus[idx, cond], sd=sqrt(sigmas2[idx]))
        })
    }))

    X <- t(vapply(seq_len(n_rows), function(idx){
        x <- t_X[idx, ]
        dropout_prob <-  vapply(seq_along(x), function(j){
            invprobit(x[j], rho[j], zeta[j])
        }, FUN.VALUE = 0.0)
        x[runif(length(x)) < dropout_prob] <- NA
        x
    }, FUN.VALUE = rep(0, length(experimental_design))))

    colnames(X) <- c(LETTERS, outer(LETTERS, LETTERS, paste0))[experimental_design]
    colnames(t_X) <- c(LETTERS, outer(LETTERS, LETTERS, paste0))[experimental_design]
    colnames(mus) <- c(LETTERS, outer(LETTERS, LETTERS, paste0))[seq_len(n_conditions)]

    return(list(X=X, t_X=t_X, mus=mus, sigmas2=sigmas2,
                changed=c(rep(FALSE, n_rows-n_changed), rep(TRUE, n_changed))))
}




