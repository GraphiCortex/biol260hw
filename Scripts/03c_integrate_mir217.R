#!/usr/bin/env Rscript

# ================================================================
# BIOL 260 — CDC14B project
# Stage 3C
#
# Query hsa-miR-217 using the same top-20% TargetScan criterion
# used for the other five candidates, then append those predictions
# to the existing Stage-3 combined dataset.
#
# Usage:
# Rscript 03c_integrate_mir217.R \
#   stage3_multimir/targetscan_combined_raw.csv \
#   stage3_complete
# ================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript 03c_integrate_mir217.R ",
      "stage3_multimir/targetscan_combined_raw.csv ",
      "stage3_complete\n"
    )
  )
}

old_file <- args[1]
outdir <- args[2]

if (!file.exists(old_file)) {
  stop("Existing TargetScan file not found: ", old_file)
}

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------
# Load multiMiR
# ------------------------------------------------

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  stop(
    "multiMiR is not installed. It should already be installed ",
    "from Stage 3."
  )
}

suppressPackageStartupMessages(
  library(multiMiR)
)


message("==============================================")
message("BIOL 260 — Stage 3C")
message("Adding miR-217 with matched TargetScan cutoff")
message("==============================================")


# ------------------------------------------------
# Read existing Stage-3 data
# ------------------------------------------------

old <- read.csv(
  old_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message(
  "[ok] existing TargetScan rows: ",
  nrow(old)
)


# ------------------------------------------------
# Query OLD TargetScan alias hsa-miR-217
# using same 20% cutoff as all other candidates
# ------------------------------------------------

message("[query] hsa-miR-217, top 20% TargetScan predictions")

result <- get_multimir(

  org = "hsa",

  mirna = "hsa-miR-217",

  table = "targetscan",

  predicted.cutoff = 20,

  predicted.cutoff.type = "p",

  predicted.site = "all",

  summary = FALSE,

  add.link = FALSE,

  use.tibble = FALSE
)


mir217 <- as.data.frame(
  result@data
)


if (nrow(mir217) == 0) {
  stop(
    "No hsa-miR-217 records returned at the 20% cutoff."
  )
}


message(
  "[ok] miR-217 TargetScan rows: ",
  nrow(mir217)
)


# ------------------------------------------------
# Add metadata from our circRNA sequence analysis
#
# Stage-2 values:
# rank = 1
# seed = ACUGCAU
# weighted seed score = 7
# ------------------------------------------------

mir217$sequence_rank <- 1

mir217$circRNA_seed_2_8 <- "ACUGCAU"

mir217$circRNA_weighted_seed_score <- 7

mir217$queried_miRNA <- "hsa-miR-217-5p"


metadata_cols <- c(
  "sequence_rank",
  "queried_miRNA",
  "circRNA_seed_2_8",
  "circRNA_weighted_seed_score"
)


mir217 <- mir217[
  ,
  c(
    metadata_cols,
    setdiff(
      colnames(mir217),
      metadata_cols
    )
  ),
  drop = FALSE
]


write.csv(

  mir217,

  file.path(
    outdir,
    "01_hsa-miR-217-5p_TargetScan.csv"
  ),

  row.names = FALSE
)


# ------------------------------------------------
# Align columns before combining
# ------------------------------------------------

all_columns <- union(
  colnames(old),
  colnames(mir217)
)


for (nm in setdiff(all_columns, colnames(old))) {
  old[, nm] <- NA
}


for (nm in setdiff(all_columns, colnames(mir217))) {
  mir217[, nm] <- NA
}


old <- old[
  ,
  all_columns,
  drop = FALSE
]

mir217 <- mir217[
  ,
  all_columns,
  drop = FALSE
]


combined <- rbind(
  mir217,
  old
)


combined <- combined[
  order(
    combined$sequence_rank
  ),
  ,
  drop = FALSE
]


write.csv(

  combined,

  file.path(
    outdir,
    "targetscan_combined_with_mir217.csv"
  ),

  row.names = FALSE
)


message("")
message(
  "[ok] complete TargetScan rows: ",
  nrow(combined)
)

message(
  "[ok] output: ",
  file.path(
    outdir,
    "targetscan_combined_with_mir217.csv"
  )
)

message("")
message("STAGE 3C COMPLETE")