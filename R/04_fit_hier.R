# fit Model 2 (hierarchical Beta-Binomial) via HMC, produce shrinkage + PPC + LOO

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(bayesplot)
  library(loo)
})
source("R/utils.R")

dir.create("results", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

dat <- load_counts()
mod <- cmdstan_model("stan/hier_beta_binomial.stan")

fit <- mod$sample(
  data            = dat$stan_data,
  seed            = 405,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.95,
  max_treedepth   = 12,
  refresh         = 250
)
fit$save_object("results/hier_fit.rds")

diag <- diag_summary(fit, vars = c("mu", "s", "theta"))

# hyperparameter summary
hyper <- fit$summary(variables = c("mu", "s"),
                     mean,
                     median,
                     ~ quantile(.x, c(0.05, 0.95), na.rm = TRUE))
print(hyper, n = 20)

# shrinkage plot
theta_mean <- fit$summary("theta", mean)$mean
df <- dat$df %>%
  mutate(theta_post = theta_mean, raw = y / n)

p_shrink <- ggplot(df, aes(raw, theta_post, color = cancer)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ gene, scales = "free") +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Partial pooling: hierarchical posterior mean vs raw rate",
    subtitle = "points below identity line have been shrunk toward gene-level mean",
    x = "raw rate y/n",
    y = expression(hat(theta)[gc])
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

ggsave(
  "figs/shrinkage.png",
  p_shrink,
  width = 9,
  height = 6,
  dpi = 150
)

# PPC
y_rep   <- fit$draws("y_rep", format = "matrix")
rate_obs <- dat$stan_data$y / dat$stan_data$n
rate_rep <- sweep(y_rep, 2, dat$stan_data$n, FUN = "/")
idx <- sample(nrow(rate_rep), 100)

p_ppc <- ppc_dens_overlay(rate_obs, rate_rep[idx, ]) +
  labs(title = "Hierarchical model: posterior predictive check", x = "rate") +
  theme_minimal()
ggsave(
  "figs/ppc_hier.png",
  p_ppc,
  width = 7,
  height = 4,
  dpi = 150
)

# LOO comparison
fit_naive <- readRDS("results/naive_fit.rds")

loo_hier  <- fit$loo(variables = "log_lik")
loo_naive <- fit_naive$loo(variables = "log_lik")

cmp <- loo_compare(list(hier = loo_hier, naive = loo_naive))
print(cmp)
saveRDS(list(
  hier = loo_hier,
  naive = loo_naive,
  compare = cmp
),
"results/loo_compare.rds")
