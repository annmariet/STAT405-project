// Hierarchical Beta-Binomial with data-driven prior on s_g for sensitivity analysis.
// prior_type: 1 = Exponential(prior_scale), 2 = Gamma(2, prior_scale), 3 = Half-Cauchy(0, prior_scale)
data {
  int<lower=1> N;
  int<lower=1> G;
  array[N] int<lower=1, upper=G> gene;
  array[N] int<lower=0> y;
  array[N] int<lower=1> n;
  int<lower=1, upper=3> prior_type;
  real<lower=0>         prior_scale;
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
  if (prior_type == 1)      s ~ exponential(prior_scale);
  else if (prior_type == 2) s ~ gamma(2, prior_scale);
  else                      s ~ cauchy(0, prior_scale);

  for (i in 1:N)
    theta[i] ~ beta(alpha_[gene[i]], beta_[gene[i]]);
  y ~ binomial(n, theta);
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = binomial_lpmf(y[i] | n[i], theta[i]);
}
