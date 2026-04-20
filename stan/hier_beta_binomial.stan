// Model 2: hierarchical Beta-Binomial with per-gene partial pooling
// Parameterization: (mu_g, s_g) -> (alpha_g, beta_g) via alpha = mu*s, beta = (1-mu)*s
data {
  int<lower=1> N;
  int<lower=1> G;
  array[N] int<lower=1, upper=G> gene;
  array[N] int<lower=0> y;
  array[N] int<lower=1> n;
}
parameters {
  vector<lower=0, upper=1>[G] mu;
  vector<lower=0>[G]          s;
  vector<lower=0, upper=1>[N] theta;
}
transformed parameters {
  vector<lower=0>[G] alpha_ = mu .* s;
  vector<lower=0>[G] beta_  = (1 - mu) .* s;
}
model {
  mu ~ beta(1, 1);
  s  ~ exponential(0.1);  // mean 10
  for (i in 1:N)
    theta[i] ~ beta(alpha_[gene[i]], beta_[gene[i]]);
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
