// Adapted from the output of the rstantools package


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>



// [[Rcpp::init]]
void rstan_additional_init(DllInfo *dll){
  R_useDynamicSymbols(dll, TRUE); // necessary for .onLoad() to work
}



