


#' Generate a data set according to the probabilistic dropout model
#'
#' Specify the number of rows in the dataset, the number of conditions and
#' replicates, how many proteins have a different mean and a few additional
#' hyperparameters and get a synthetic dataset the is similar to data from a
#' real label-free mass spectrometry experiment.
#'
#' @param n_rows integer. The number of rows in the new dataset
#' @param experimental_design a vector that specifies which samples belong
#'   to the same condition. Default: `NULL` in which case `n_replicates` must
#'   be specified
#' @param n_replicates integer or vector. The number of replicates in each
#'   condition.
#' @param n_conditions The number of conditions. Setting `n_replicates=3` and
#'   `n_conditions=2` is equal to specifying `experimental_design=c(1,1,1,2,2,2)`.
#' @param frac_changed the fraction of rows for which different means are drawn
#'   for each conditon.
#' @param n_changed alternative way to specify for how many rows have different
#'   means in each condition.
#' @param mu0 the global mean around which the row means are drawn. Default `20`
#' @param sigma20 the global variance specifying the spread of means around `mu0`.
#'   Default `10`.
#' @param nu degrees of freedom for the the global variance prior. Default `3`.
#' @param eta scale of the global variance prior. Default `0.3`.
#' @param rho vector specifying the intensity where the chance of a dropout is
#'   50/50. Either length one or same length as `n_replicates * n_conditons` or
#'   `length(experimental_design)` respectively. Default `18`.
#' @param zeta vector specifying the scale of the dropout curve.
#'   Either length one or same length as `n_replicates * n_conditons` or
#'   `length(experimental_design)` respectively. Default `18`.
#'
#'
#' @return a list with 5 elements
#'   \describe{
#'     \item{X}{the data matrix with missing values}
#'     \item{t_X}{the true data matrix, before data dropped out}
#'     \item{mus}{matrix of size `n_rows * n_conditions`. The true means
#'      for each condition}
#'     \item{sigmas2}{a vector of size `n_rows`. The true variance for each row.}
#'     \item{changed}{a boolean vector of size `n_rows`, with the label if a row
#'     has different means for each condition}
#'   }
#'
#' @examples
#'  data <- generate_synthetic_data(n_rows=10,
#'                 n_replicates=3, n_conditions=2)
#'
#'  data2 <- generate_synthetic_data(n_rows=10,
#'                 experimental_design=c(1,1,1,2,2,2))
#'
#'  data3 <- generate_synthetic_data(n_rows=10,
#'                 rep(letters[1:3], each=4))
#'
#'
#' @export
generate_synthetic_data <- function(
            n_rows,
            experimental_design=NULL,
            n_replicates=as.numeric(table(experimental_design)),
            n_conditions=length(n_replicates),
            frac_changed=0.1,
            n_changed=round(n_rows * min(1, frac_changed)),
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

    if(is.null(experimental_design)){
        experimental_design <- c(unlist(mapply(rep, seq_len(n_conditions), n_replicates)))
        conditons_names <- c(LETTERS, outer(LETTERS, LETTERS, paste0))[experimental_design]
    }else{
        stopifnot(c(table(experimental_design)) == n_replicates)
        experimental_design_fct <- as.factor(experimental_design)
        conditons_names <- levels(experimental_design_fct)
        experimental_design <- as.numeric(experimental_design_fct)
    }

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

    t_X <- do.call(cbind, lapply(experimental_design, function(cond){
        rnorm(n_rows, mus[, cond], sd=sqrt(sigmas2))
    }))

    X <- t(vapply(seq_len(n_rows), function(idx){
        x <- t_X[idx, ]
        dropout_prob <-  vapply(seq_along(x), function(j){
            invprobit(x[j], rho[j], zeta[j])
        }, FUN.VALUE = 0.0)
        x[runif(length(x)) < dropout_prob] <- NA
        x
    }, FUN.VALUE = rep(0, length(experimental_design))))

    colnames(X) <- conditons_names[experimental_design]
    colnames(t_X) <- conditons_names[experimental_design]
    colnames(mus) <- conditons_names[seq_len(n_conditions)]

    return(list(X=X, t_X=t_X, mus=mus, sigmas2=sigmas2,
                changed=c(rep(FALSE, n_rows-n_changed), rep(TRUE, n_changed))))
}




