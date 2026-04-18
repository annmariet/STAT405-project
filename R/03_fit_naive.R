# fit Model 1 (naive independent Beta-Binomial) via HMC

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(ggplot2)
  library(bayesplot)
})
source("R/utils.R")

dir.create("results", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

dat <- load_counts()
mod <- cmdstan_model("stan/naive_beta_binomial.stan")

fit <- mod$sample(
  data            = dat$stan_data,
  seed            = 405,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 250
)

fit$save_object("results/naive_fit.rds")

# diagnostics
diag <- diag_summary(fit, vars = "theta")

# posterior predictive check: density overlay of y_rep vs y
y_rep <- fit$draws("y_rep", format = "matrix")
y_obs <- dat$stan_data$y
n_obs <- dat$stan_data$n

# work on rates so genes with different n are comparable
rate_obs <- y_obs / n_obs
rate_rep <- sweep(y_rep, 2, n_obs, FUN = "/")
idx <- sample(nrow(rate_rep), 100)

p <- ppc_dens_overlay(rate_obs, rate_rep[idx, ]) +
  labs(title = "Naive model: posterior predictive check (mutation rate)", x = "rate") +
  theme_minimal()

ggsave(
  "figs/ppc_naive.png",
  p,
  width = 7,
  height = 4,
  dpi = 150
)
