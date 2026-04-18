suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

# 10 oncogenes / tumor suppressors (from proposal)
GENES <- c("TP53", "KRAS", "PIK3CA", "BRAF", "PTEN",
           "APC", "ARID1A", "EGFR", "IDH1", "CDKN2A")

# 20 TCGA PanCancer Atlas studies with N >= 150 sequenced samples
STUDY_IDS <- c(
  "brca_tcga_pan_can_atlas_2018",  # breast
  "luad_tcga_pan_can_atlas_2018",  # lung adenocarcinoma
  "lusc_tcga_pan_can_atlas_2018",  # lung squamous
  "skcm_tcga_pan_can_atlas_2018",  # skin melanoma
  "coadread_tcga_pan_can_atlas_2018",  # colorectal
  "stad_tcga_pan_can_atlas_2018",  # stomach
  "ucec_tcga_pan_can_atlas_2018",  # uterine endometrial
  "ov_tcga_pan_can_atlas_2018",    # ovarian
  "blca_tcga_pan_can_atlas_2018",  # bladder
  "hnsc_tcga_pan_can_atlas_2018",  # head & neck
  "kirc_tcga_pan_can_atlas_2018",  # kidney clear cell
  "kirp_tcga_pan_can_atlas_2018",  # kidney papillary
  "lihc_tcga_pan_can_atlas_2018",  # liver
  "prad_tcga_pan_can_atlas_2018",  # prostate
  "thca_tcga_pan_can_atlas_2018",  # thyroid
  "lgg_tcga_pan_can_atlas_2018",   # low-grade glioma
  "gbm_tcga_pan_can_atlas_2018",   # glioblastoma
  "paad_tcga_pan_can_atlas_2018",  # pancreatic
  "cesc_tcga_pan_can_atlas_2018",  # cervical
  "sarc_tcga_pan_can_atlas_2018"   # sarcoma
)

study_to_cancer <- function(sid) gsub("_tcga_pan_can_atlas_2018", "", sid)

# pull (gene, sample) mutation counts for one study, with retry
pull_cbio <- function(sid, genes = GENES, sleep = 0.5, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      samples <- cbioportalR::available_samples(study_id = sid)
      n_total <- nrow(samples)

      muts <- cbioportalR::get_mutations_by_sample(
        sample_id = samples$sampleId,
        study_id  = sid,
        genes     = genes
      )

      counts <- muts %>%
        distinct(hugoGeneSymbol, sampleId) %>%
        count(hugoGeneSymbol, name = "y")

      tibble(gene = genes) %>%
        left_join(counts, by = c("gene" = "hugoGeneSymbol")) %>%
        mutate(y      = replace_na(y, 0),
               n      = n_total,
               study  = sid,
               cancer = study_to_cancer(sid))
    }, error = function(e) {
      message(sprintf("[%s] attempt %d failed: %s", sid, attempt, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      Sys.sleep(sleep)
      return(res)
    }
    Sys.sleep(2 ^ attempt)
  }
  stop(sprintf("pull_cbio failed for %s after %d attempts", sid, max_retries))
}

# load processed counts and build the Stan data list
load_counts <- function(rds = "data/processed/mutation_counts.rds") {
  d <- readRDS(rds)
  list(
    df = d,
    stan_data = list(
      N    = nrow(d),
      G    = max(d$gene_id),
      gene = d$gene_id,
      y    = d$y,
      n    = d$n
    )
  )
}

# extract MCMC diagnostics from a cmdstanr fit into a tidy tibble
diag_summary <- function(fit, vars = NULL) {
  s <- fit$summary(variables = vars)
  d <- fit$diagnostic_summary()
  list(
    summary = s,
    n_divergent     = sum(d$num_divergent),
    n_max_treedepth = sum(d$num_max_treedepth),
    ebfmi           = d$ebfmi
  )
}
