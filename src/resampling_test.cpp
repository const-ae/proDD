#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double resampling_test_impl(NumericVector x, NumericVector y, R_xlen_t nmax) {

    nmax = std::min(nmax, x.size() * y.size());

    long pval = 0;
    long counter = 0;
    for(int i = 0; i < x.size(); i++){
        for(int j = 0; j < y.size(); j++){
            pval += x[i] > y[j];
            if(++counter >= nmax) goto leave_loop;
            if (counter % 100000 == 0)
                Rcpp::checkUserInterrupt();
        }
    }
    leave_loop:

        if(pval == nmax){
            return 1.0 - 1.0 / nmax;
        }else if(pval == 0){
            return 1.0 / nmax;
        }else{
            return pval * 1.0 / nmax;
        }
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
set.seed(1)
x <- rnorm(100, mean=0)
y <- rnorm(100, mean=2)
resampling_test(x, y, alternative= "less")
mean(outer(x, y, `>`))
t.test(x, y, alternative = "less")$p.value
*/
