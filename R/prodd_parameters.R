

# This file provides a few useful helper function for handling the prodd_parameter function

is_valid_prodd_parameters <- function(x, ...){
    base_elements_present <- ! is.null(x$converged) &&  ! is.null(x$error) && ! is.null(x$hyper_params) &&
        ! is.null(x$feature_params) && ! is.null(x$experimental_design)
    hyper_elements_present <- ! is.null(x$hyper_params$mu0) && ! is.null(x$hyper_params$sigma20) &&
        ! is.null(x$hyper_params$rho) && ! is.null(x$hyper_params$zeta) &&
        ! is.null(x$hyper_params$nu) && ! is.null(x$hyper_params$eta)
    feature_elements_present <- ! is.null(x$feature_params$mup) && ! is.null(x$feature_params$sigma2mup) &&
        ! is.null(x$feature_params$sigma2p)
    feature_elements_right_size <- is.matrix(x$feature_params$mup) && is.matrix(x$feature_params$sigma2mup) &&
        nrow(x$feature_params$mup) == nrow(x$feature_params$sigma2mup) &&
        ncol(x$feature_params$mup) == ncol(x$feature_params$sigma2mup) &&
        length(x$feature_params$sigma2p) == nrow(x$feature_params$mup) &&
        length(unique(x$experimental_design)) == ncol(x$feature_params$mup)
    base_elements_present && hyper_elements_present && feature_elements_present && feature_elements_right_size
}

#' @rdname print.prodd_parameters
#' @export
format.prodd_parameters <- function(x, ...){
    stopifnot(is_valid_prodd_parameters(x))

    header <- "\tParameters of the probabilistic dropout model\n"

    exp_txt <- paste0("There were ", length(x$experimental_design), " samples",
                      " in ", length(unique(x$experimental_design)), " conditons.",
                      " In total there were ", nrow(x$feature_params$mup), " rows.")

    if(x$converged){
        converged_txt <- "The model has successfully converged."
    }else{
        converged_txt <- "Attention: the model has not converged."
    }

    error_txt <- paste0("The error of the last iteration was ", formatC(x$error, digits=4, format="g"))

    hyper_para_txt <- paste0("\nThe inferred parameters are:\n",
     paste0(vapply(seq_along(x$hyper_params), function(idx){
        paste0(names(x$hyper_params)[idx], ":",
               paste0(rep(" ", times=9-nchar(names(x$hyper_params)[idx])), collapse=""),
               paste0(formatC(x$hyper_params[[idx]], digits=3, width=1, format="g"), collapse=", "))
    }, FUN.VALUE = ""), collapse = "\n"))


    paste0(c(header, exp_txt, converged_txt, error_txt, hyper_para_txt), collapse = "\n")
}

#' Collection of parameters generated when fitting the probabilistic dropout model
#'
#' @param x an object of calss  prodd_parameters
#' @param ... additional arguments to be passed to or from the method
#'
#' @export
print.prodd_parameters <- function(x, ...){
    cat(format(x, ...), "\n")
    invisible(x)
}


is.prodd_parameters <- function(x) inherits(x, "prodd_parameters")


#' Update the experimental_design of a parameter object
#'
#' Helper method to transform the parameters from one experimental_design
#' to another. It avoids expensive re-calculation of the parameters and
#' instead uses a averaging of the parameters to produce the new values.
#'
#' The method only works for changing the condition assignment of samples
#' or removing samples. The method cannot add new samples and infer the
#' parameters for them.
#'
#' The method is for example useful if you want to calculate the distance
#' with unbiased feature_parameters that don't already include the information
#' to which condition a sample belongs.
#'
#' @param params an object of class `prodd_parameters` which is returned by the
#'   \code{fit_parameters()} function.
#' @param new_experimental_design a vector that assignes each sample to one condition.
#'   It has to have the same length as the old experimental_design. It can either be
#'   a factor, a character or a numeric vector. Each unique element is one condition.
#'   To remove samples code them as \code{NA}.
#' @return the updated list of class prodd_parameters
#'
#' @export
transform_parameters <- function(params, new_experimental_design){

    stopifnot(is.prodd_parameters(params))
    stopifnot(is_valid_prodd_parameters(params))
    stopifnot(length(new_experimental_design) == length(params$experimental_design))

    new_experimental_design_fct <- as.factor(new_experimental_design)
    new_experimental_design <- as.numeric(new_experimental_design_fct)
    N_cond <- length(levels(new_experimental_design_fct)) # To make sure to ignore NA's

    old_experimental_design <- params$experimental_design

    helper_mup <- params$feature_params$mup[ ,old_experimental_design, drop=FALSE]
    new_mup <- t(mply_dbl(seq_len(N_cond), ncol=nrow(helper_mup), function(condition){
        rowMeans(helper_mup[, which(new_experimental_design == condition),drop=FALSE])
    }))

    helper_sigma2mup <- params$feature_params$sigma2mup[ ,old_experimental_design, drop=FALSE]
    new_sigma2mup <- t(mply_dbl(seq_len(N_cond), ncol=nrow(helper_mup), function(condition){
        rowMeans(helper_sigma2mup[, which(new_experimental_design == condition),drop=FALSE])
    }))

    colnames(new_mup) <- levels(new_experimental_design_fct)
    colnames(new_sigma2mup) <- levels(new_experimental_design_fct)
    rownames(new_mup) <- rownames(params$feature_params$mup)
    rownames(new_sigma2mup) <- rownames(params$feature_params$sigma2mup)


    hyper_params <- params$hyper_params
    hyper_params$rho <- hyper_params$rho[! is.na(new_experimental_design)]
    hyper_params$zeta <- hyper_params$zeta[! is.na(new_experimental_design)]

    ret <- list(hyper_params = hyper_params,
                feature_params = list(
                    mup=new_mup,
                    sigma2p=params$feature_params$sigma2p,
                    sigma2mup=new_sigma2mup
                ),
                experimental_design=new_experimental_design_fct[! is.na(new_experimental_design)],
                error=params$error, converged=params$converged)

    class(ret) <- "prodd_parameters"
    ret
}



