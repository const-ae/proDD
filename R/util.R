

#' @useDynLib proDD, .registration = TRUE
#' @import rstantools
#' @import Rcpp
#' @import methods
#' @importFrom Rcpp sourceCpp loadModule
#' @importFrom rstan sampling
NULL

#' @import stats
NULL


#' apply function that always returns a numeric matrix
#'
#' The function is modelled after `vapply`, but always returns a matrix
#' with one row for each iteration. You need to provide the number
#' of elements each function call produces beforehand (i.e. the number of
#' resulting columns). For a more flexible version where you don't need to
#' provide the number of columns see \code{\link{msply_dbl}}
#' @describeIn mply_dbl apply function that always returns a numeric matrix
#'
#' @param x a vector that will be passed to `vapply` or a matrix that will be
#'   passed to apply with \code{MARGIN=1}.
#' @param FUN the function that returns a vector of length ncol
#' @param ncol the length of the vector returned by `FUN`.
#' @param ... additional arguments to FUN
#'
#' @return a matrix of size \code{length(x) x ncol}
mply_dbl <- function(x, FUN, ncol=1, ...){
    if(is.vector(x)){
        res <- vapply(x, FUN, FUN.VALUE=rep(0.0, times=ncol), ...)
    }else{
        res <- apply(x, 1, FUN, ...) * 1.0
        if((ncol == 1 && ! is.vector(res)) || (ncol > 1 && nrow(res) != ncol)){
            stop(paste0("values must be length ", ncol,
                       ", but result is length ", nrow(res)))
        }
    }

    if(ncol == 1){
        as.matrix(res, nrow=length(res), ncol=1)
    }else{
        t(res)
    }
}

#' @describeIn mply_dbl flexible version that automatically infers the number
#'   of columns
msply_dbl <- function(x, FUN, ...){
    if(is.vector(x)){
        res <- sapply(x, FUN, ...) * 1.0
    }else{
        res <- apply(x, 1, FUN, ...) * 1.0
    }

    if(is.matrix(res)){
        t(res)
    }else{
        as.matrix(res, nrow=length(res))
    }
}

