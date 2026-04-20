# prior sensitivity: refit hierarchical model with 3 alternative priors on s_g
# primary: Exp(0.1)  mean 10, weakly informative
# tighter: Gamma(2, 0.1)  mean 20, pulls toward stronger pooling
# diffuse: Half-Cauchy(25) heavy tails, lets data dominate

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
mod <- cmdstan_model("stan/hier_flex_prior.stan")

configs <- tibble(
  label       = c(
    "Primary Exp(0.1)",
    "Tighter Gamma(2, 0.1)",
    "Diffuse Half-Cauchy(25)"
  ),
  prior_type  = c(1, 2, 3),
  prior_scale = c(0.1, 0.1, 25)
)

fits <- lapply(seq_len(nrow(configs)), function(i) {
  message("fitting: ", configs$label[i])
  mod$sample(
    data = c(
      dat$stan_data,
      list(
        prior_type  = configs$prior_type[i],
        prior_scale = configs$prior_scale[i]
      )
    ),
    seed            = 405,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    adapt_delta     = 0.95,
    refresh         = 0
  )
})
names(fits) <- configs$label

# summarize mu and s for each config
summ <- bind_rows(lapply(seq_along(fits), function(i) {
  fit <- fits[[i]]
  s <- fit$summary(c("mu", "s"), mean, median, ~ quantile(.x, c(0.05, 0.95), na.rm = TRUE))
  s %>% mutate(config = configs$label[i])
}))

# extract gene index and param name
summ <- summ %>%
  mutate(
    param = sub("\\[.*", "", variable),
    gene_id = as.integer(sub(".*\\[(\\d+)\\].*", "\\1", variable)),
    gene    = GENES[gene_id],
    config  = factor(config, levels = configs$label)
  )

saveRDS(summ, "results/prior_sensitivity.rds")

# posterior mean + 90% CI for mu_g under each config
summ_mu <- filter(summ, param == "mu")

p <- ggplot(summ_mu, aes(mean, gene, color = config)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `5%`, xmax = `95%`),
                 position = position_dodge(width = 0.5),
                 height = 0) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  labs(
    title = "Prior sensitivity: posterior of gene-level mean rate mu_g",
    subtitle = "points = posterior mean; bars = 90% CI under each s_g prior",
    x = expression(mu[g]),
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(
  "figs/prior_sensitivity.png",
  p,
  width = 8,
  height = 5,
  dpi = 150
)

# print summary table
mu_wide <- summ_mu %>%
  select(gene, config, mean) %>%
  pivot_wider(names_from = config, values_from = mean) %>%
  arrange(desc(`Primary Exp(0.1)`))
print(mu_wide)

s_wide <- filter(summ, param == "s") %>%
  select(gene, config, mean) %>%
  pivot_wider(names_from = config, values_from = mean)
print(s_wide)
