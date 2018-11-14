

# This file provides a few useful helper function for handling the prodd_parameter function

is_valid_prodd_parameters <- function(x, ...){
    ! is.null(x$converged) &&  ! is.null(x$error) && ! is.null(x$hyper_params) &&
        ! is.null(x$feature_params) && ! is.null(x$experimental_design)
}

#' @rdname print.prodd_parameters
#' @export
format.prodd_parameters <- function(x, ...){
    stopifnot(is_valid_prodd_parameters(x))

    header <- "\tParameters of the probabilistic dropout model\n"

    exp_txt <- paste0("There were ", length(experimental_design), " samples",
                      " in ", length(unique(experimental_design)), " conditons.",
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
