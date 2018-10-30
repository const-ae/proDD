




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


#' Calculate the variance of the mean estimator
#'
fit_feature_mean_uncertainties <- function(X, rho, zeta, nu, eta, mu0, sigma20,
                                           experimental_design){
    N_cond <- length(unique(experimental_design))

    mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            sel <- which(experimental_design == cond & is.na(X[idx,]))
            nobs <- sum(! is.na(x))
            if(nobs == 0){
                sigma2obsreg <- if(nu > 2) nu * eta / (nu-2) else nu * eta / (nu + 2)

                sigma2mu <- sigma20
                mu <- mu0
            }else if(nobs == 1) {
                sigma2obsreg <- if(nu > 2) nu * eta / (nu-2) else nu * eta / (nu + 2)

                sigma2mu <- 1/(1/sigma20 + 1/(sigma2obsreg/nobs))
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu+ mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }else{
                sigma2obs <- var(x, na.rm=TRUE)
                sigma2obsreg <- ((nobs-1) * sigma2obs + nu * eta) / (nobs-1 + nu)

                sigma2mu <- 1/(1/sigma20 + 1/(sigma2obsreg/nobs))
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu + mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }
            zetastar <- zeta * sqrt(1 + sigma2obsreg/zeta^2)
            variance_probdropout(mu, sigma2=sigma2mu, rho=rho[sel], zeta=zetastar[sel])
        }, FUN.VALUE=0.0)
    })
}


#' Estimate the mean for each feature
#'
#'
fit_feature_means <- function(X, sigma2p, sigma2mus, rho, zeta, nu, eta,
                              mu0, sigma20, experimental_design){
    N_cond <- length(unique(experimental_design))

    mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            x <- X[idx, which(experimental_design == cond)]
            sel <- which(experimental_design == cond & is.na(X[idx,]))
            nobs <- sum(! is.na(x))
            if(nobs == 0){
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- mu0
            }else if(nobs == 1){
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu + mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }else{
                sigma2 <- sigma2p[idx]
                sigma2mu <- sigma2mus[idx, cond]
                mu <- (mean(x, na.rm=TRUE) * 1/sigma2mu+ mu0 * 1/sigma20) / (1/sigma2mu + 1/sigma20)
            }

            zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
            mean_probdropout(mu, sigma2=sigma2mu, rho=rho[sel], zeta=zetastar[sel])
        }, FUN.VALUE=0.0)
    })
}






