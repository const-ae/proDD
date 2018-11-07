










#' Iteratively update the hyperparameter until convergence
#'
#'
fit_hyperparameters <- function(X, experimental_design,
                                dropout_curve_calc=c("global", "global_scale", "sample"),
                                frac_subsample=1.0, n_subsample=round(nrow(X) * frac_subsample),
                                max_iter=10, epsilon=1e-3, verbose=FALSE){

    dropout_curve_calc <- match.arg(dropout_curve_calc, c("global", "global_scale", "sample"))

    experimental_design_fct <- as.factor(experimental_design)
    experimental_design <- as.numeric(experimental_design_fct)

    N_cond <- length(unique(experimental_design))
    if(n_subsample >= nrow(X)){
        # do nothing ...
    }else if(n_subsample > 0 && n_subsample < nrow(X)){
        X_bckp <- X
        sel <- sample(seq_len(nrow(X_bckp)), n_subsample)
        X <- X_bckp[sel, ,drop=FALSE]
    }else{
        stop(paste0("Illegal argument n_subsample it must be larger than 0"))
    }


    #Initialize the parameters
    mup <- mply_dbl(seq_len(nrow(X)), ncol=N_cond, function(idx){
        vapply(seq_len(N_cond), function(cond){
            mean(X[idx, which(experimental_design == cond)], na.rm=TRUE)
        }, FUN.VALUE=0.0)
    })
    mup[is.na(mup)] <- min(mup, na.rm=TRUE)
    mu0 <- mean(mup, na.rm=TRUE)

    # Estimate variances
    sigma2p <-  apply(X, 1, function(row){
        res <- sapply(seq_len(N_cond), function(cond){
            x <- row[which(experimental_design == cond)]
            nobs <- sum(! is.na(x))
            if(nobs > 1){
                c((nobs-1) * var(x, na.rm=TRUE), nobs-1)
            }else{
                c(NA, NA)
            }
        })
        sum(res[1, ], na.rm=TRUE) / sum(res[2, ], na.rm=TRUE)
    })
    sigma20 <- mean(sigma2p, na.rm=TRUE)
    sigma2p[is.na(sigma2p)] <- sigma20
    sigma2_unreg <- sigma2p

    # Sigmoid parameters
    zeta <- rep(- sqrt(sigma20),  ncol(X))
    rho <- sapply(seq_len(ncol(X)), function(colidx){
        sel <- which(is.na(c(X[, colidx])))
        optimize(function(rho){
            sum(invprobit(X[, colidx], rho, zeta[colidx], log=TRUE, oneminus = TRUE), na.rm=TRUE) +
              sum(invprobit(mup[sel, experimental_design[colidx]], rho,
                       zeta[colidx] * sqrt(1 + sigma2p[sel] / zeta[colidx]^2), log=TRUE))
        }, lower=-100, upper=100, maximum=TRUE)$maximum
    })



    # Variance prior
    nu <- 3
    eta <- sigma20 + 0.1

    last_round_params <- list(eta,nu,mu0,sigma20,rho,zeta)
    converged <- FALSE
    iter <- 1

    while(! converged && iter < max_iter){
        if(verbose) message(paste0("Starting iter ", iter))

        # Estimate dropout sigmoid
        sigmoid_est <- if(dropout_curve_calc == "global") {
            fit_global_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }else if(dropout_curve_calc == "sample"){
            fit_sample_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }else if(dropout_curve_calc == "global_scale"){
            fit_global_scale_dropout_curves(X, mup, sigma2p, experimental_design,
                                      prev_zeta = zeta, prev_rho=rho)
        }
        rho <- sigmoid_est$rho
        zeta <- sigmoid_est$zeta

        # Make a good estimate for each mean
        sigma2p <- fit_feature_variances(X, mup, rho, zeta, nu, eta, experimental_design)
        mu_vars <- fit_feature_mean_uncertainties(X, rho, zeta, nu, eta, mu0,
                                                    sigma20, experimental_design)
        mup <- fit_feature_means(X, sigma2p, mu_vars, rho, zeta, nu, eta, mu0,
                                 sigma20, experimental_design)

        # Fit the variance prior
        var_est <- fit_variance_prior(X, rho, zeta, experimental_design)
        nu <- var_est$df_prior
        eta <- var_est$var_prior

        # Fit the location prior
        loc_est <- fit_location_prior(X, mup, zeta, rho, experimental_design)
        mu0 <- loc_est$mu0
        sigma20 <- loc_est$sigma20


        if(verbose){
            message(paste0("eta0 estimate: ", round(eta, 3)))
            message(paste0("nu0 estimate: ", round(nu, 3)))
            message(paste0("mu0 estimate: ", round(mu0, 3)))
            message(paste0("sigma20 estimate: ", round(sigma20, 3)))
            message(paste0("rho estimate: ", paste0(round(rho, 3), collapse = ",")))
            message(paste0("zeta estimate: ", paste0(round(zeta, 3), collapse = ",")))
        }
        error <- sum(mapply(function(new, old){sum(new - old)/length(new)},
                            list(eta,nu,mu0,sigma20,rho,zeta), last_round_params )^2)
        if(verbose) message(paste0("Error: ", error))
        if(error < epsilon) {
            if(verbose) message("converged!")
            converged <- TRUE
        }
        last_round_params <- list(eta,nu,mu0,sigma20,rho,zeta)
        iter <- iter + 1
    }

    names(last_round_params) <- c("eta", "nu", "mu0", "sigma20", "rho", "zeta")
    list(params = last_round_params,
         error=error, converged=converged)

}



