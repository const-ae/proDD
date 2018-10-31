




#' Column wise median normalization of the data matrix
#'
#' The method calculates for each sample the median change (i.e. the difference
#' between the observed value and the row average) and substracts it from
#' each row. Missing values are ignored in the procedure. The method is based
#' on the assumption that a majority of the rows did not change.
#'
#'
#' @param X a matrix of proteins x samples
#' @return the normalized matrix
#' @export
median_normalization <- function(X){
    stopifnot(is.matrix(X))
    Xnorm <- X
    for(idx in 1:ncol(X)){
        Xnorm[, idx, drop=FALSE] <- X[, idx, drop=FALSE] -
            median(X[, idx, drop=FALSE] - rowMeans(X, na.rm=TRUE), na.rm=TRUE )
    }
    Xnorm
}



