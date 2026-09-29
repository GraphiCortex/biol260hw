#!/usr/bin/env Rscript

# ================================================================
# BIOL 260 — CDC14B circRNA project
# Stage 3: conserved TargetScan target extraction via Bioconductor
#
# INPUT:
#   stage2/sequence_shortlist.csv
#
# OUTPUT:
#   stage3_bioc/
#       targetscan_family_map.csv
#       targetscan_conserved_targets_all.csv
#       targetscan_conserved_targets_summary.csv
#       unmapped_seed_families.csv
#
# IMPORTANT:
# targetscan.Hs.eg.db contains CONSERVED human TargetScan targets.
# A miRNA being absent here does NOT mean that it has no targets.
# ================================================================


# ------------------------------------------------
# 1. COMMAND-LINE ARGUMENTS
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript 03_targetscan_via_bioconductor.R ",
      "stage2/sequence_shortlist.csv stage3_bioc\n"
    )
  )
}

shortlist_file <- args[1]
outdir <- args[2]

if (!file.exists(shortlist_file)) {
  stop("Shortlist file not found: ", shortlist_file)
}

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("[start] TargetScan Stage 3")


# ------------------------------------------------
# 2. INSTALL REQUIRED PACKAGES
# ------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  message("[install] BiocManager")

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "targetscan.Hs.eg.db"
)


for (pkg in packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("[install] ", pkg)

    try(
      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      ),
      silent = TRUE
    )
  }
}


# Fallback for the older TargetScan annotation package
if (!requireNamespace("targetscan.Hs.eg.db", quietly = TRUE)) {

  message(
    "[fallback] Trying archived targetscan.Hs.eg.db package"
  )

  archive_url <- paste0(
    "https://bioconductor.posit.co/packages/3.19/",
    "data/annotation/src/contrib/",
    "targetscan.Hs.eg.db_0.6.1.tar.gz"
  )

  try(
    install.packages(
      archive_url,
      repos = NULL,
      type = "source"
    ),
    silent = TRUE
  )
}


for (pkg in packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      "\nCould not install/load package: ",
      pkg,
      "\nSend the full terminal output back to ChatGPT."
    )
  }
}


# ------------------------------------------------
# 3. LOAD PACKAGES
# ------------------------------------------------

suppressPackageStartupMessages({

  library(AnnotationDbi)

  library(org.Hs.eg.db)

  library(targetscan.Hs.eg.db)
})


message("[ok] packages loaded")


# ------------------------------------------------
# 4. HELPER FUNCTIONS
# ------------------------------------------------

clean_seed <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  chartr(
    "T",
    "U",
    x
  )
}


clean_mirna <- function(x) {

  tolower(
    gsub(
      "_",
      "-",
      trimws(
        as.character(x)
      )
    )
  )
}


collapse_values <- function(x) {

  x <- unique(
    x[
      !is.na(x) &
      nzchar(x)
    ]
  )

  paste(
    sort(x),
    collapse = ";"
  )
}


safe_lookup <- function(key, map_object) {

  result <- tryCatch(

    mget(
      key,
      map_object,
      ifnotfound = list(NA)
    )[[1]],

    error = function(e) NA
  )

  result
}


# ------------------------------------------------
# 5. READ STAGE-2 SHORTLIST
# ------------------------------------------------

shortlist <- read.csv(
  shortlist_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "seed_family_rank",
  "seed_2_8",
  "representative_miRNA",
  "family_members",
  "weighted_seed_score"
)


missing_columns <- setdiff(
  required_columns,
  colnames(shortlist)
)


if (length(missing_columns) > 0) {

  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


shortlist$seed_2_8 <- clean_seed(
  shortlist$seed_2_8
)


message(
  "[ok] Stage-2 seed families loaded: ",
  nrow(shortlist)
)

# ------------------------------------------------
# 6. BUILD HUMAN TARGETSCAN miRNA LOOKUP
# ------------------------------------------------

mir_ids <- AnnotationDbi::mappedkeys(
  targetscan.Hs.egMIRNA
)

message(
  "[info] Total TargetScan miRNA records available: ",
  length(mir_ids)
)

# Human miRBase identifiers use the hsa- prefix.
human_mir_ids <- mir_ids[
  grepl(
    "^hsa-",
    mir_ids,
    ignore.case = TRUE
  )
]

message(
  "[info] Human hsa-* records: ",
  length(human_mir_ids)
)

if (length(human_mir_ids) == 0) {
  stop(
    "No hsa-* TargetScan miRNA identifiers were found."
  )
}


lookup_rows <- list()

counter <- 0L


for (mid in human_mir_ids) {

  info <- safe_lookup(
    mid,
    targetscan.Hs.egMIRNA
  )


  if (
    length(info) == 1 &&
    all(is.na(info))
  ) {
    next
  }


  fam <- safe_lookup(
    mid,
    targetscan.Hs.egMIRBASE2FAMILY
  )


  if (
    length(fam) == 1 &&
    all(is.na(fam))
  ) {
    next
  }


  seed <- tryCatch(
    as.character(
      info$Seed.m8
    ),
    error = function(e) NA_character_
  )


  conservation <- tryCatch(
    as.character(
      info$Family.conservation
    ),
    error = function(e) NA_character_
  )


  mature_seq <- tryCatch(
    as.character(
      info$Mature.sequence
    ),
    error = function(e) NA_character_
  )


  # Skip malformed records without a usable seed.
  if (
    length(seed) == 0 ||
    is.na(seed[1]) ||
    !nzchar(seed[1])
  ) {
    next
  }


  counter <- counter + 1L


  lookup_rows[[counter]] <- data.frame(

    MiRBase_ID =
      mid,

    Seed_m8 =
      clean_seed(
        seed[1]
      ),

    TargetScan_family =
      as.character(
        fam[1]
      ),

    Family_conservation =
      if (
        length(conservation) > 0
      ) {
        conservation[1]
      } else {
        NA_character_
      },

    Mature_sequence =
      if (
        length(mature_seq) > 0
      ) {
        mature_seq[1]
      } else {
        NA_character_
      },

    stringsAsFactors = FALSE
  )
}


if (length(lookup_rows) == 0) {

  stop(
    "Human hsa-* records exist, but none could be parsed into ",
    "TargetScan families and seeds."
  )
}


mir_lookup <- do.call(
  rbind,
  lookup_rows
)


message(
  "[ok] human TargetScan miRNAs parsed: ",
  nrow(mir_lookup)
)


message(
  "[info] unique human TargetScan seeds: ",
  length(
    unique(
      mir_lookup$Seed_m8
    )
  )
)


message(
  "[info] example records:"
)

print(
  head(
    mir_lookup,
    5
  )
)
# ------------------------------------------------
# 7. MAP OUR SEED FAMILIES TO TARGETSCAN FAMILIES
# ------------------------------------------------

family_map_rows <- list()


for (i in seq_len(nrow(shortlist))) {

  row <- shortlist[i, ]


  seed <- clean_seed(
    row$seed_2_8
  )


  seed_hits <- mir_lookup[
    mir_lookup$Seed_m8 == seed,
    ,
    drop = FALSE
  ]


  member_names <- unlist(
    strsplit(
      row$family_members,
      ";",
      fixed = TRUE
    )
  )


  member_names <- unique(
    clean_mirna(
      member_names
    )
  )


  exact_hits <- mir_lookup[
    clean_mirna(
      mir_lookup$MiRBase_ID
    ) %in% member_names,
    ,
    drop = FALSE
  ]


  mapped_families <- unique(
    c(
      seed_hits$TargetScan_family,
      exact_hits$TargetScan_family
    )
  )


  mapped_families <- mapped_families[
    !is.na(mapped_families) &
    nzchar(mapped_families)
  ]


  mapped_ids <- unique(
    c(
      seed_hits$MiRBase_ID,
      exact_hits$MiRBase_ID
    )
  )


  conservation <- unique(
    c(
      seed_hits$Family_conservation,
      exact_hits$Family_conservation
    )
  )


  family_map_rows[[i]] <- data.frame(

    seed_family_rank =
      row$seed_family_rank,

    seed_2_8 =
      seed,

    representative_miRNA =
      row$representative_miRNA,

    family_members =
      row$family_members,

    weighted_seed_score =
      row$weighted_seed_score,

    TargetScan_families =
      collapse_values(
        mapped_families
      ),

    n_TargetScan_families =
      length(
        mapped_families
      ),

    mapped_miRBase_IDs =
      collapse_values(
        mapped_ids
      ),

    conservation_labels =
      collapse_values(
        conservation
      ),

    stringsAsFactors = FALSE
  )
}


family_map <- do.call(
  rbind,
  family_map_rows
)


write.csv(
  family_map,
  file.path(
    outdir,
    "targetscan_family_map.csv"
  ),
  row.names = FALSE
)


message(
  "[ok] seed families mapped to TargetScan: ",
  sum(
    family_map$n_TargetScan_families > 0
  )
)


# ------------------------------------------------
# 8. WRITE UNMAPPED FAMILIES
# ------------------------------------------------

unmapped <- family_map[
  family_map$n_TargetScan_families == 0,
  ,
  drop = FALSE
]


write.csv(
  unmapped,
  file.path(
    outdir,
    "unmapped_seed_families.csv"
  ),
  row.names = FALSE
)


message(
  "[ok] unmapped seed families: ",
  nrow(unmapped)
)


# ------------------------------------------------
# 9. REVERSE TARGETSCAN GENE → miRNA MAP
#
# targetscan.Hs.egTARGETS maps:
#
#     Entrez gene → TargetScan miRNA family
#
# revmap converts it to:
#
#     TargetScan miRNA family → Entrez genes
# ------------------------------------------------

family_to_genes <- AnnotationDbi::revmap(
  targetscan.Hs.egTARGETS
)


target_rows <- list()

target_counter <- 0L


# ------------------------------------------------
# 10. EXTRACT TARGET GENES
# ------------------------------------------------

for (i in seq_len(nrow(family_map))) {

  fm <- family_map[i, ]


  if (
    fm$n_TargetScan_families == 0
  ) {
    next
  }


  families <- unlist(
    strsplit(
      fm$TargetScan_families,
      ";",
      fixed = TRUE
    )
  )


  families <- families[
    nzchar(families)
  ]


  for (family in families) {


    entrez_ids <- safe_lookup(
      family,
      family_to_genes
    )


    entrez_ids <- unique(
      as.character(
        entrez_ids
      )
    )


    entrez_ids <- entrez_ids[
      !is.na(entrez_ids) &
      nzchar(entrez_ids)
    ]


    if (length(entrez_ids) == 0) {
      next
    }


    gene_symbols <- AnnotationDbi::mapIds(

      org.Hs.eg.db,

      keys = entrez_ids,

      keytype = "ENTREZID",

      column = "SYMBOL",

      multiVals = "first"
    )


    for (j in seq_along(entrez_ids)) {

      target_counter <- target_counter + 1L


      eid <- entrez_ids[j]


      target_rows[[target_counter]] <- data.frame(

        seed_family_rank =
          fm$seed_family_rank,

        seed_2_8 =
          fm$seed_2_8,

        representative_miRNA =
          fm$representative_miRNA,

        weighted_seed_score =
          fm$weighted_seed_score,

        TargetScan_family =
          family,

        Entrez_ID =
          eid,

        Gene_Symbol =
          unname(
            gene_symbols[eid]
          ),

        stringsAsFactors = FALSE
      )
    }
  }
}


# ------------------------------------------------
# 11. COMBINE TARGET RESULTS
# ------------------------------------------------

if (length(target_rows) > 0) {

  targets <- do.call(
    rbind,
    target_rows
  )

} else {

  targets <- data.frame(

    seed_family_rank =
      integer(),

    seed_2_8 =
      character(),

    representative_miRNA =
      character(),

    weighted_seed_score =
      numeric(),

    TargetScan_family =
      character(),

    Entrez_ID =
      character(),

    Gene_Symbol =
      character(),

    stringsAsFactors = FALSE
  )
}


# remove duplicate gene-family combinations
targets <- unique(
  targets
)


targets <- targets[
  order(
    targets$seed_family_rank,
    targets$TargetScan_family,
    targets$Gene_Symbol
  ),
  ,
  drop = FALSE
]


write.csv(
  targets,
  file.path(
    outdir,
    "targetscan_conserved_targets_all.csv"
  ),
  row.names = FALSE
)


message(
  "[ok] TargetScan conserved target rows: ",
  nrow(targets)
)


# ------------------------------------------------
# 12. SUMMARIZE NUMBER OF TARGETS PER SEED FAMILY
# ------------------------------------------------

summary_rows <- list()


for (i in seq_len(nrow(family_map))) {

  fm <- family_map[i, ]


  sub <- targets[
    targets$seed_family_rank ==
      fm$seed_family_rank,
    ,
    drop = FALSE
  ]


  unique_genes <- unique(
    sub$Gene_Symbol[
      !is.na(
        sub$Gene_Symbol
      )
    ]
  )


  summary_rows[[i]] <- data.frame(

    seed_family_rank =
      fm$seed_family_rank,

    seed_2_8 =
      fm$seed_2_8,

    representative_miRNA =
      fm$representative_miRNA,

    weighted_seed_score =
      fm$weighted_seed_score,

    TargetScan_families =
      fm$TargetScan_families,

    conservation_labels =
      fm$conservation_labels,

    n_conserved_target_genes =
      length(
        unique_genes
      ),

    stringsAsFactors = FALSE
  )
}


target_summary <- do.call(
  rbind,
  summary_rows
)


write.csv(
  target_summary,
  file.path(
    outdir,
    "targetscan_conserved_targets_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------
# 13. FINISHED
# ------------------------------------------------

message("")
message("==============================================")
message("STAGE 3 COMPLETE")
message("==============================================")

message(
  "Mapped TargetScan families: ",
  sum(
    family_map$n_TargetScan_families > 0
  )
)

message(
  "Unmapped seed families: ",
  nrow(unmapped)
)

message(
  "Conserved TargetScan target records: ",
  nrow(targets)
)

message("")
message(
  "Outputs written to: ",
  normalizePath(
    outdir,
    winslash = "/",
    mustWork = FALSE
  )
)

message("")
message(
  "NOTE: This is the conserved TargetScan layer only."
)

message(
  "Absence here does not imply absence of biological targets."
)