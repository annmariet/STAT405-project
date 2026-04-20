# pull mutation counts for 10 genes x 20 TCGA pan-cancer studies via cBioPortal API
# run once; output cached to data/raw/. report knit does not call this script.

suppressPackageStartupMessages({
  library(cbioportalR)
  library(dplyr)
  library(jsonlite)
})
source("R/utils.R")

out_csv      <- "data/raw/cbioportal_pull_2026-04-14.csv"
out_manifest <- "data/raw/pull_manifest.json"

set_cbioportal_db("public")

message(sprintf("pulling %d genes x %d studies", length(GENES), length(STUDY_IDS)))
raw <- bind_rows(lapply(STUDY_IDS, function(sid) {
  message("  ", sid)
  pull_cbio(sid, genes = GENES)
}))

write.csv(raw, out_csv, row.names = FALSE)

manifest <- list(
  pulled_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  cbioportalR_ver  = as.character(packageVersion("cbioportalR")),
  R_version        = R.version.string,
  genes            = GENES,
  study_ids        = STUDY_IDS,
  n_observations   = nrow(raw),
  n_genes          = length(GENES),
  n_studies        = length(STUDY_IDS)
)
write_json(manifest,
           out_manifest,
           pretty = TRUE,
           auto_unbox = TRUE)

