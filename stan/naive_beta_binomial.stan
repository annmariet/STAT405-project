// Model 1: independent Beta(1,1) prior per (gene, cancer)
data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  array[N] int<lower=1> n;
}
parameters {
  vector<lower=0, upper=1>[N] theta;
}
model {
  theta ~ beta(1, 1);
  y ~ binomial(n, theta);
}
generated quantities {
  array[N] int y_rep;
  vector[N] log_lik;
  for (i in 1:N) {
    y_rep[i]   = binomial_rng(n[i], theta[i]);
    log_lik[i] = binomial_lpmf(y[i] | n[i], theta[i]);
  }
}
