

#' apply function that always returns a numeric matrix
#'
#' The function is modelled after `vapply`, but always returns a matrix
#' with one row for each iteration. You need to provide the number
#' of elements each function call produces beforehand (i.e. the number of
#' resulting columns). For a more flexible version where you don't need to
#' provide the number of columns see \code{\link{msply_dbl}}
#' @describeIn mply_dbl apply function that always returns a numeric matrix
#'
#' @param x a vector that will be passed to `vapply`
#' @param FUN the function that returns a vector of length ncol
#' @param ncol the length of the vector returned by `FUN`.
#'
#' @return a matrix of size \code{length(x) x ncol}
mply_dbl <- function(x, FUN, ncol=1, ...){
    res <- vapply(x, FUN, FUN.VALUE=rep(0.0, times=ncol), ...)
    if(ncol == 1){
        as.matrix(res, nrow=length(res), ncol=1)
    }else{
        t(res)
    }
}

#' @describeIn mply_dbl flexible version that automatically infers the number
#'   of columns
msply_dbl <- function(x, FUN, ...){
    res <- sapply(x, FUN, ...)
    if(is.matrix(res)){
        t(res)
    }else{
        as.matrix(res, nrow=length(res))
    }
}

