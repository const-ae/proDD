data {
  int nsamples;
  int nrows;
  int ncond;
  int totalmissing;
  matrix[nrows, nsamples] X;
  int experimental_design[nsamples];

  real zeta[nsamples];
  real rho[nsamples];
  real mu0;
  real sigma20;
  real eta;
  real nu;
}
parameters {
  matrix[nrows, ncond] mu;
  real<lower=0> sigma2[nrows];
}
transformed parameters{
  real zetastar[totalmissing];
  {
    int counter = 1;
    for(i in 1:nrows){
      for(j in 1:nsamples){
        if(is_inf(X[i, j])){
          zetastar[counter] = zeta[j] * sqrt(1 + sigma2[i]/zeta[j]^2);
          counter = counter + 1;
        }
      }
    }
  }
}
model {
  for(c in 1:ncond){
    mu[,c] ~ normal(mu0, sqrt(sigma20));
  }

  sigma2 ~ scaled_inv_chi_square(nu, sqrt(eta));
  {
    int counter = 1;
    for(i in 1:nrows){
      for(j in 1:nsamples){
        if(is_inf(X[i, j])){
          target += normal_lccdf(mu[i, experimental_design[j]] | rho[j], fabs(zetastar[counter]));
          counter += 1;
        }else{
          X[i, j] ~ normal(mu[i, experimental_design[j]], sqrt(sigma2[i]));
        }
      }
    }
  }

}
