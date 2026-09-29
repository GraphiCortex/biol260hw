#!/usr/bin/env Rscript

# ================================================================
# BIOL 260 — CDC14B circRNA project
# Stage 3: TargetScan 8 target extraction through multiMiR
#
# This script:
#   1. reads the Stage-2 sequence shortlist
#   2. takes the top N sequence-ranked miRNA families
#   3. queries ONLY the TargetScan table in multiMiR
#   4. keeps the top 20% of TargetScan predictions
#   5. exports raw TargetScan results for downstream pathway analysis
#
# INPUT:
#   stage2/sequence_shortlist.csv
#
# OUTPUT:
#   stage3_multimir/
#
# Usage:
#   Rscript 03_targetscan_multimir.R \
#       stage2/sequence_shortlist.csv \
#       stage3_multimir \
#       6
#
# Third argument = number of sequence-ranked candidates to query.
# Default = 6.
# ================================================================


# ------------------------------------------------
# 1. ARGUMENTS
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {

  stop(
    paste0(
      "\nUsage:\n",
      "Rscript 03_targetscan_multimir.R ",
      "stage2/sequence_shortlist.csv ",
      "stage3_multimir ",
      "[top_n]\n"
    )
  )
}


shortlist_file <- args[1]
outdir <- args[2]

top_n <- if (length(args) >= 3) {
  as.integer(args[3])
} else {
  6L
}


if (!file.exists(shortlist_file)) {
  stop("Input file does not exist: ", shortlist_file)
}


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


message("==============================================")
message("BIOL 260 — TARGETSCAN STAGE 3")
message("==============================================")
message("")


# ------------------------------------------------
# 2. INSTALL multiMiR
# ------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


if (!requireNamespace("multiMiR", quietly = TRUE)) {

  message("[install] multiMiR")

  BiocManager::install(
    "multiMiR",
    ask = FALSE,
    update = FALSE
  )
}


if (!requireNamespace("multiMiR", quietly = TRUE)) {

  stop(
    "Could not install multiMiR."
  )
}


suppressPackageStartupMessages(
  library(multiMiR)
)


message("[ok] multiMiR loaded")


# ------------------------------------------------
# 3. RECORD DATABASE PROVENANCE
# ------------------------------------------------

message("[info] retrieving multiMiR database information...")


db_info <- tryCatch(

  multimir_dbInfo(),

  error = function(e) {

    warning(
      "Could not retrieve database metadata: ",
      conditionMessage(e)
    )

    NULL
  }
)


if (!is.null(db_info)) {

  write.csv(
    db_info,
    file.path(
      outdir,
      "multimir_database_info.csv"
    ),
    row.names = FALSE
  )


  ts_info <- db_info[
    tolower(db_info$map_name) == "targetscan",
    ,
    drop = FALSE
  ]


  if (nrow(ts_info) > 0) {

    message("")
    message("[TargetScan database metadata]")

    print(ts_info)

    message("")
  }
}


# ------------------------------------------------
# 4. READ STAGE-2 SHORTLIST
# ------------------------------------------------

shortlist <- read.csv(
  shortlist_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required <- c(
  "seed_family_rank",
  "seed_2_8",
  "representative_miRNA",
  "weighted_seed_score"
)


missing <- setdiff(
  required,
  colnames(shortlist)
)


if (length(missing) > 0) {

  stop(
    "Missing columns: ",
    paste(
      missing,
      collapse = ", "
    )
  )
}


shortlist <- shortlist[
  order(
    shortlist$seed_family_rank
  ),
  ,
  drop = FALSE
]


top_n <- min(
  top_n,
  nrow(shortlist)
)


candidates <- shortlist[
  seq_len(top_n),
  ,
  drop = FALSE
]


write.csv(
  candidates,
  file.path(
    outdir,
    "candidates_queried.csv"
  ),
  row.names = FALSE
)


message(
  "[ok] sequence-ranked candidates selected: ",
  nrow(candidates)
)


print(
  candidates[
    ,
    c(
      "seed_family_rank",
      "representative_miRNA",
      "seed_2_8",
      "weighted_seed_score"
    )
  ]
)


# ------------------------------------------------
# 5. QUERY TARGETSCAN ONE miRNA AT A TIME
# ------------------------------------------------
#
# We deliberately query one miRNA at a time because:
#   - it is more robust to network timeouts
#   - each candidate gets its own reproducible CSV
#   - one failed query will not destroy the whole run
#
# We use:
#   table = "targetscan"
#   predicted.site = "all"
#   predicted.cutoff = 20
#   predicted.cutoff.type = "p"
#
# Therefore we retain the top 20% of TargetScan predictions.
# This is a pre-specified filter, not manual gene cherry-picking.
# ------------------------------------------------


all_results <- list()

failed_rows <- list()

result_counter <- 0L
failed_counter <- 0L


for (i in seq_len(nrow(candidates))) {

  mir <- candidates$representative_miRNA[i]

  rank <- candidates$seed_family_rank[i]

  seed <- candidates$seed_2_8[i]

  seed_score <- candidates$weighted_seed_score[i]


  message("")
  message("----------------------------------------------")
  message(
    "[query] rank ",
    rank,
    ": ",
    mir
  )
  message("----------------------------------------------")


  result <- tryCatch(

    get_multimir(

      org = "hsa",

      mirna = mir,

      table = "targetscan",

      predicted.cutoff = 20,

      predicted.cutoff.type = "p",

      predicted.site = "all",

      summary = FALSE,

      add.link = FALSE,

      use.tibble = FALSE
    ),

    error = function(e) {

      message(
        "[error] ",
        conditionMessage(e)
      )

      NULL
    }
  )


  if (is.null(result)) {

    failed_counter <- failed_counter + 1L

    failed_rows[[failed_counter]] <- data.frame(

      seed_family_rank = rank,

      representative_miRNA = mir,

      seed_2_8 = seed,

      reason = "query_error",

      stringsAsFactors = FALSE
    )

    next
  }


  dat <- tryCatch(

    as.data.frame(
      result@data
    ),

    error = function(e) {

      message(
        "[error] could not extract result@data: ",
        conditionMessage(e)
      )

      NULL
    }
  )


  if (
    is.null(dat) ||
    nrow(dat) == 0
  ) {

    message(
      "[warning] no TargetScan records returned"
    )


    failed_counter <- failed_counter + 1L

    failed_rows[[failed_counter]] <- data.frame(

      seed_family_rank = rank,

      representative_miRNA = mir,

      seed_2_8 = seed,

      reason = "no_TargetScan_records",

      stringsAsFactors = FALSE
    )

    next
  }


  # Add our circRNA-stage metadata.
  dat$sequence_rank <- rank

  dat$circRNA_seed_2_8 <- seed

  dat$circRNA_weighted_seed_score <- seed_score

  dat$queried_miRNA <- mir


  # Move our metadata to the front.
  metadata_cols <- c(
    "sequence_rank",
    "queried_miRNA",
    "circRNA_seed_2_8",
    "circRNA_weighted_seed_score"
  )


  dat <- dat[
    ,
    c(
      metadata_cols,
      setdiff(
        colnames(dat),
        metadata_cols
      )
    ),
    drop = FALSE
  ]


  safe_name <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    mir
  )


  outfile <- file.path(
    outdir,
    paste0(
      sprintf(
        "%02d",
        rank
      ),
      "_",
      safe_name,
      "_TargetScan.csv"
    )
  )


  write.csv(
    dat,
    outfile,
    row.names = FALSE
  )


  message(
    "[ok] TargetScan records: ",
    nrow(dat)
  )

  message(
    "[ok] wrote: ",
    outfile
  )


  result_counter <- result_counter + 1L

  all_results[[result_counter]] <- dat
}


# ------------------------------------------------
# 6. COMBINE ALL SUCCESSFUL TARGETSCAN RESULTS
# ------------------------------------------------

if (length(all_results) > 0) {

  all_columns <- unique(
    unlist(
      lapply(
        all_results,
        colnames
      )
    )
  )


  aligned <- lapply(

    all_results,

    function(x) {

      missing_cols <- setdiff(
        all_columns,
        colnames(x)
      )


      for (nm in missing_cols) {

        x[[nm]] <- NA
      }


      x[
        ,
        all_columns,
        drop = FALSE
      ]
    }
  )


  combined <- do.call(
    rbind,
    aligned
  )


  write.csv(

    combined,

    file.path(
      outdir,
      "targetscan_combined_raw.csv"
    ),

    row.names = FALSE
  )


  message("")
  message(
    "[ok] combined TargetScan rows: ",
    nrow(combined)
  )

} else {

  message("")
  message(
    "[warning] no successful TargetScan queries"
  )
}


# ------------------------------------------------
# 7. WRITE FAILED / UNMAPPED CANDIDATES
# ------------------------------------------------

if (length(failed_rows) > 0) {

  failed <- do.call(
    rbind,
    failed_rows
  )

} else {

  failed <- data.frame(

    seed_family_rank = integer(),

    representative_miRNA = character(),

    seed_2_8 = character(),

    reason = character(),

    stringsAsFactors = FALSE
  )
}


write.csv(

  failed,

  file.path(
    outdir,
    "targetscan_failed_or_unmapped.csv"
  ),

  row.names = FALSE
)


# ------------------------------------------------
# 8. DONE
# ------------------------------------------------

message("")
message("==============================================")
message("STAGE 3 COMPLETE")
message("==============================================")

message(
  "Candidates attempted: ",
  nrow(candidates)
)

message(
  "Successful TargetScan queries: ",
  length(all_results)
)

message(
  "Failed / unmapped candidates: ",
  nrow(failed)
)

message("")
message(
  "Output directory: ",
  normalizePath(
    outdir,
    winslash = "/",
    mustWork = FALSE
  )
)