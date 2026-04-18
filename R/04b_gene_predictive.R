# gene-level posterior predictive: for each gene g, draw theta* ~ Beta(alpha_g, beta_g)
# integrating over the posterior of (alpha_g, beta_g). This shows the distribution of
# mutation rates we'd expect for a NEW cancer type.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
})
source("R/utils.R")

fit <- readRDS("results/hier_fit.rds")
dat <- load_counts()

draws  <- fit$draws(c("alpha_", "beta_"), format = "draws_matrix")
alpha_ <- as.matrix(draws[, grepl("^alpha_", colnames(draws))])
beta_  <- as.matrix(draws[, grepl("^beta_", colnames(draws))])

n_draws <- nrow(alpha_)
G       <- ncol(alpha_)
n_per   <- 50

gene_pp <- do.call(rbind, lapply(1:G, function(g) {
  idx <- sample(n_draws, n_per)
  theta_new <- unlist(lapply(idx, function(i)
    rbeta(20, alpha_[i, g], beta_[i, g])))
  data.frame(gene = GENES[g], theta = theta_new)
}))

observed <- dat$df %>% transmute(gene, rate = y / n)

p <- ggplot(gene_pp, aes(theta)) +
  geom_density(fill = "#4c72b0",
               alpha = 0.35,
               color = "#4c72b0") +
  geom_rug(
    data = observed,
    aes(rate),
    color = "firebrick",
    alpha = 0.7,
    length = unit(0.05, "npc")
  ) +
  facet_wrap( ~ factor(gene, levels = GENES), scales = "free_y") +
  labs(
    title = "Gene-level posterior predictive distribution of mutation rate",
    subtitle = "blue: Beta(alpha_g, beta_g) integrated over posterior | red rug: observed rates",
    x = expression(theta),
    y = "density"
  ) +
  theme_minimal(base_size = 10)

ggsave(
  "figs/gene_predictive.png",
  p,
  width = 9,
  height = 6,
  dpi = 150
)
