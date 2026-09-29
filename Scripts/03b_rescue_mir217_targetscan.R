#!/usr/bin/env Rscript

# ================================================================
# BIOL 260 — CDC14B circRNA project
# Stage 3B: rescue / nomenclature check for miR-217 in TargetScan
#
# Purpose:
# Stage 3 returned no TargetScan records for hsa-miR-217-5p.
# Older TargetScan releases often use mature-miRNA names without
# the 5p/3p suffix, so we test several aliases explicitly.
# ================================================================


# ------------------------------------------------
# 1. Output directory
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

outdir <- if (length(args) >= 1) {
  args[1]
} else {
  "stage3b_mir217"
}

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------
# 2. Load multiMiR
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

  stop("Could not install/load multiMiR.")
}


suppressPackageStartupMessages(
  library(multiMiR)
)


message("==============================================")
message("BIOL 260 — miR-217 TargetScan rescue")
message("==============================================")


# ------------------------------------------------
# 3. Test common aliases
# ------------------------------------------------

aliases <- c(
  "hsa-miR-217-5p",
  "hsa-miR-217",
  "miR-217-5p",
  "miR-217"
)


all_results <- list()

status_rows <- list()


for (i in seq_along(aliases)) {

  alias <- aliases[i]


  message("")
  message("----------------------------------------------")
  message("[query] ", alias)
  message("----------------------------------------------")


  result <- tryCatch(

    get_multimir(

      org = "hsa",

      mirna = alias,

      table = "targetscan",

      predicted.cutoff = 300000,

      predicted.cutoff.type = "n",

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

    status_rows[[i]] <- data.frame(

      alias = alias,

      status = "query_error",

      n_records = 0,

      stringsAsFactors = FALSE
    )

    next
  }


  dat <- tryCatch(

    as.data.frame(
      result@data
    ),

    error = function(e) NULL
  )


  if (
    is.null(dat) ||
    nrow(dat) == 0
  ) {

    message(
      "[result] 0 TargetScan records"
    )


    status_rows[[i]] <- data.frame(

      alias = alias,

      status = "no_records",

      n_records = 0,

      stringsAsFactors = FALSE
    )

    next
  }


  message(
    "[result] ",
    nrow(dat),
    " TargetScan records"
  )


  dat$queried_alias <- alias


  all_results <- append(
  all_results,
  list(dat)
)


  safe_name <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    alias
  )


  write.csv(

    dat,

    file.path(
      outdir,
      paste0(
        safe_name,
        "_TargetScan_all.csv"
      )
    ),

    row.names = FALSE
  )


  status_rows[[i]] <- data.frame(

    alias = alias,

    status = "success",

    n_records = nrow(dat),

    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------
# 4. Write alias-status table
# ------------------------------------------------

status <- do.call(
  rbind,
  status_rows
)


write.csv(

  status,

  file.path(
    outdir,
    "mir217_alias_status.csv"
  ),

  row.names = FALSE
)


# ------------------------------------------------
# 5. Combine successful TargetScan queries
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

      missing_columns <- setdiff(
        all_columns,
        colnames(x)
      )


      for (nm in missing_columns) {

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


  # Remove duplicate records produced when multiple aliases
  # resolve to the same underlying TargetScan record.

  comparison_columns <- setdiff(
    colnames(combined),
    "queried_alias"
  )


  row_keys <- apply(

    combined[
      ,
      comparison_columns,
      drop = FALSE
    ],

    1,

    function(x) {

      paste(
        x,
        collapse = "|||"
      )
    }
  )


  combined_unique <- combined[
    !duplicated(row_keys),
    ,
    drop = FALSE
  ]


  write.csv(

    combined_unique,

    file.path(
      outdir,
      "mir217_TargetScan_combined_unique.csv"
    ),

    row.names = FALSE
  )


  message("")
  message(
    "[ok] unique rescued TargetScan records: ",
    nrow(combined_unique)
  )

} else {

  message("")
  message(
    "[warning] none of the miR-217 aliases returned TargetScan records"
  )
}


message("")
message("==============================================")
message("STAGE 3B COMPLETE")
message("==============================================")

message(
  "Alias-status file: ",
  file.path(
    outdir,
    "mir217_alias_status.csv"
  )
)