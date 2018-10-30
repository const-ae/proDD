




#' Find the most likely variance explaining the values for each row in X
#'
#'
fit_feature_variances <- function(X, mup, rho, zeta, nu, eta, experimental_design,
                                  upper=1e3){
    N_cond <- length(unique(experimental_design))

    vapply(seq_len(nrow(X)), function(idx){

        optimize(f=function(sigma2){
            reg <- extraDistr::dinvchisq(sigma2, nu, eta, log=TRUE)
            reg + sum(vapply(seq_len(N_cond), function(cond){
                # All the missing values for this condition
                sel <- which(experimental_design == cond & is.na(X[idx,]))
                x <- X[idx, which(experimental_design == cond)]
                zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
                if(sum(! is.na(x)) == 0){
                    norm <- 0
                }else{
                    norm <- sum(dnorm(x, mean=mup[idx, cond], sd=sqrt(sigma2), log=TRUE), na.rm=TRUE)
                }

                cumulative <- sum(invprobit(mup[idx, cond], rho[sel], zetastar[sel], log=TRUE))
                norm + cumulative
            }, FUN.VALUE=0.0))
        }, lower=0, upper=upper, maximum = TRUE)$maximum
    }, FUN.VALUE = 0.0)

}









