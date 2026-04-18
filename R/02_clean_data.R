# clean raw cBioPortal pull -> processed mutation_counts.rds with integer indices for Stan

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/utils.R")

raw <- read_csv("data/raw/cbioportal_pull_2026-04-14.csv", show_col_types = FALSE)

d <- raw %>%
  filter(n >= 150) %>%
  mutate(
    rate     = round(y / n, 4),
    gene_id   = as.integer(factor(gene,   levels = GENES)),
    cancer_id = as.integer(factor(cancer, levels = unique(cancer)))
  ) %>%
  arrange(gene_id, cancer_id) %>%
  select(gene, cancer, gene_id, cancer_id, y, n, rate, study)

stopifnot(
  nrow(d) > 100,
  all(d$y <= d$n),
  all(d$y >= 0),
  all(d$gene_id %in% 1:length(GENES))
)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(d, "data/processed/mutation_counts.rds")

