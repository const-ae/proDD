#' Test the difference between x and y using resampling.
#'
#' The method is an alternative to \code{mean(outer(x, y, `<`))}
#' which doesn't allocate memory and is thus much faster.
#'
#' @param x first numeric vector
#' @param y second numeric vector
#' @param nmax the maximum number of comparisons. Can be used if the full
#'   precision is not needed. Defaults to making all possible comparisons
#' @param alternative how are x and y compared. Same principles as
#'   \code{t.test}
#' @param mu test if there is at least distance mu between x and y.
#'   Equivalent to calling \code{resampling_test(x, y + mu)}
#' @export
resampling_test <- function(x, y, nmax = length(x) * length(y),
                            alternative=c("two.sided", "less", "greater"),
                            mu = 0) {
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
