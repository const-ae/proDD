#' Test the difference between x and y using resampling.
#'
#' The method is an alternative to \code{mean(outer(x, y, `<`))}
#' which doesn't allocate memory and is thus much faster.
#'
#' @param x first numeric vector
#' @param y second numeric vector
#' @param alternative how are x and y compared. Same principles as
#'   \code{t.test}
#' @param mu test if there is at least distance mu between x and y.
#'   Equivalent to calling \code{resampling_test(x, y + mu)}
#' @param nmax the maximum number of comparisons. Can be used if the full
#'   precision is not needed. Defaults to making all possible comparisons
#'
#' @return the p-value, ie. the fraction of elements in x that are smaller than y.
#' @export
resampling_test <- function(x, y,
                            alternative=c("two.sided", "less", "greater"),
                            mu = 0,  nmax = length(x) * length(y)) {
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    if(alternative == "less"){
        resampling_test_impl(x, y + mu, nmax)
    }else if(alternative == "greater"){
        resampling_test_impl(y+ mu, x, nmax)
    }else{
        res <- resampling_test_impl(y + mu, x, nmax)
        min(res, 1-res) * 2
    }

}





#' Test difference between two posterior distribution
#'
#' The function calculates the probability that x is less, greater or just
#' different from y. In addition the function calculates the False sign or
#' small rate (practically equivalent to the false discovery rate) and
#' estimates the mean differrence between x and y.
#'
#' @param x a numeric matrix of with one row for each protein and one column
#'   for each MCMC sample. Usually obtained the result of a call to
#'   \code{sample_posterior_means()}.
#' @param y another numeric matrix of with one row for each protein and one column
#'   for each MCMC sample. Usually obtained the result of a call to
#'   \code{sample_posterior_means()}. It must have the same number of rows
#'   as \code{x}.
#' @param sort boolean indicating if the result is ordered. Default: FALSE. If
#'   TRUE the result is ordered by adj_pval. Alternatively you can provide the
#'   column name that is used for ordering the resulting data.frame.
#' @inheritParams resampling_test
#' @param pval_adjust_method the method used to adjusted the p-value. See
#'   \code{\link{p.adjust}}.
#' @param quantiles a vector of numbers between 0 and 1 that indicate which
#'   quantiles of the difference are calculated in addition to the mean. This
#'   is helpful to get an impression of the uncertainty of the difference.
#'   By default the 95\% credibility interval is calculated ie.
#'   \code{quantiles=c(0.025, 0.975)}.
#'
#' @return a data.frame with one row per protein and the following columns
#'   \describe{
#'     \item{name}{the names of each protein extracted from the rownames of x
#'       and y}
#'     \item{pval}{the probability that x is less/greater/different (depending on)
#'       the value of \code{alternative} than y. As the package is following a Bayesian
#'       paradigm this not strictly a p-value, but it for practical reasons
#'       can be considered as one.}
#'     \item{adj_pval}{the multiple testing adjusted p-values using the
#'       \code{pval_adjust_method}. Again in the Bayesian paradigm it can also
#'       be called the False Sign or Small rate.}
#'     \item{diff}{an estimate of the mean difference between x and y.}
#'     \item{`quantile_0.025`}{an estiamte of the upper bound of the difference
#'       estimate. The actual number of quantile columns is determined by
#'       \code{quantiles} argument}
#'     \item{`quantile_0.975`}{an estiamte of the upper bound of the difference
#'       estimate}
#'   }
#'
#' @seealso \code{\link{resampling_test}}
#'
#' @export
test_diff <- function(x, y, sort=FALSE, alternative=c("two.sided", "less", "greater"),
                      mu = 0,  nmax = ncol(x) * ncol(y),
                      pval_adjust_method="BH", quantiles=c(0.025, 0.975)) {

    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    stopifnot(nrow(x) == nrow(y))
    stopifnot(is.null(rownames(x)) || is.null(rownames(y)) ||
                  all(rownames(x) == rownames(y)))

    pvals <- vapply(seq_len(nrow(x)), function(idx){
        resampling_test(x[idx, , drop=FALSE], y[idx, , drop=FALSE],
                        alternative=alternative, mu=mu, nmax=nmax)
    }, FUN.VALUE = 0.0)
    adj_pvals <- p.adjust(pvals, method = pval_adjust_method)
    mean_diff <- rowMeans(x - y)
    quantiles_diff <- mply_dbl(x-y, quantile, probs=quantiles, ncol=length(quantiles))
    means <- rowMeans((x+y)/2)

    names <- if(! is.null(rownames(x))){
        rownames(x)
    }else if(! is.null(rownames(y))){
        rownames(y)
    }else{
        as.character(seq_len(nrow(x)))
    }
    result <- data.frame(names, pvals, adj_pvals, mean_diff, quantiles_diff, means,
                         stringsAsFactors = FALSE)
    colnames <- c("name", "pval", "adj_pval", "diff",
                  (if(length(quantiles) >= 1) paste0("diff_quantile_", quantiles) else NULL),
                  "mean")
    colnames(result) <- colnames

    if(sort == FALSE){
        result
    }else if(sort == TRUE){
        result[order(result$pval), ]
    }else{
        stopifnot(sort %in% colnames(result))
        result[order(result[, sort]), ]

    }
}


