# second posterior computation method: Stan ADVI (mean-field variational inference)
# fits the SAME hierarchical model as HMC; compares marginal posteriors for key parameters.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
source("R/utils.R")

dir.create("results", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

dat <- load_counts()
mod <- cmdstan_model("stan/hier_beta_binomial.stan")

# mean-field: factorized Gaussian approximation on unconstrained space
fit_mf <- mod$variational(
  data      = dat$stan_data,
  seed      = 405,
  algorithm = "meanfield",
  output_samples = 4000
)
fit_mf$save_object("results/advi_meanfield_fit.rds")

# full-rank: multivariate Gaussian (captures parameter correlations)
fit_fr <- mod$variational(
  data      = dat$stan_data,
  seed      = 405,
  algorithm = "fullrank",
  output_samples = 4000
)
fit_fr$save_object("results/advi_fullrank_fit.rds")

# load HMC reference
fit_hmc <- readRDS("results/hier_fit.rds")

# compare marginal densities
# TP53 (gene_id = 1) is most interesting: high mean + spread
# thyroid (cancer = "thca") has TP53 y=2, n=500 -> low rate, informative for shrinkage
df <- dat$df
thca_tp53_idx <- which(df$gene == "TP53" & df$cancer == "thca")

pars <- c("mu[1]", "s[1]", paste0("theta[", thca_tp53_idx, "]"))

extract_long <- function(fit, method) {
  as_draws_df(fit$draws(pars)) %>%
    as_tibble() %>%
    select(all_of(pars)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(method = method)
}

long <- bind_rows(
  extract_long(fit_hmc, "HMC"),
  extract_long(fit_mf, "ADVI (mean-field)"),
  extract_long(fit_fr, "ADVI (full-rank)")
)

param_labels <- c("mu[1]"  = "mu[TP53]  (gene-level mean)", "s[1]"   = "s[TP53]  (gene-level concentration)")
param_labels[paste0("theta[", thca_tp53_idx, "]")] <- "theta[TP53, thyroid]"

long <- long %>%
  mutate(parameter = factor(parameter, levels = pars, labels = param_labels[pars]))

p <- ggplot(long, aes(value, color = method, fill = method)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  facet_wrap( ~ parameter, scales = "free", ncol = 3) +
  scale_color_manual(
    values = c(
      "HMC" = "#1f77b4",
      "ADVI (mean-field)" = "#d62728",
      "ADVI (full-rank)"  = "#2ca02c"
    )
  ) +
  scale_fill_manual(
    values  = c(
      "HMC" = "#1f77b4",
      "ADVI (mean-field)" = "#d62728",
      "ADVI (full-rank)"  = "#2ca02c"
    )
  ) +
  labs(
    title = "Posterior marginals: HMC vs ADVI (mean-field, full-rank)",
    subtitle = "mean-field tracks HMC closely; full-rank fails to converge despite higher flexibility",
    x = "value",
    y = "density"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(
  "figs/hmc_vs_advi.png",
  p,
  width = 9,
  height = 4.5,
  dpi = 150
)

# numerical summary: SD and 90% CI width for each method
summ <- long %>%
  group_by(parameter, method) %>%
  summarise(
    mean   = mean(value),
    sd     = sd(value),
    q05    = quantile(value, 0.05),
    q95    = quantile(value, 0.95),
    ci90_w = q95 - q05,
    .groups = "drop"
  )

print(summ, n = 30)

saveRDS(summ, "results/hmc_vs_advi_summary.rds")
