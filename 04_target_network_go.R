#!/usr/bin/env Rscript

# ================================================================
# BIOL 260 — CDC14B circRNA project
# Stage 4: gene-level TargetScan consolidation + GO-BP convergence
#
# INPUT
#   stage3_multimir/targetscan_combined_raw.csv
#
# OUTPUT
#   stage4/
#     targetscan_gene_level.csv
#     candidate_target_summary.csv
#     target_overlap_counts.csv
#     target_overlap_jaccard.csv
#     go_bp_enrichment_human_background.csv
#     go_bp_enrichment_targetscan_union_background.csv
#     go_bp_top_terms.csv
#     cancer_theme_summary.csv
#
# Biological logic
# ----------------
# 1) Collapse multiple TargetScan records for the same miRNA-gene pair.
# 2) Preserve the most favorable (lowest / most negative) TargetScan score.
# 3) Test GO Biological Process enrichment.
# 4) Use a PRE-SPECIFIED set of cancer/tumor-initiation themes rather than
#    choosing pathways after seeing the answer.
#
# IMPORTANT
# ---------
# TargetScan predictions are computational hypotheses, not proof of direct
# regulation. GO enrichment is descriptive pathway convergence, not causality.
# ================================================================


# ------------------------------------------------
# 1. Arguments
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript 04_target_network_go.R ",
      "stage3_multimir/targetscan_combined_raw.csv ",
      "stage4\n"
    )
  )
}

input_file <- args[1]
outdir <- args[2]

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("==============================================")
message("BIOL 260 — STAGE 4")
message("Target-network and GO convergence analysis")
message("==============================================")


# ------------------------------------------------
# 2. Required Bioconductor annotation packages
# ------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

needed <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db"
)

for (pkg in needed) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("[install] ", pkg)

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

for (pkg in needed) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Could not install/load package: ", pkg)
  }
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

message("[ok] annotation packages loaded")


# ------------------------------------------------
# 3. Helpers
# ------------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

collapse_values <- function(x) {

  x <- unique(
    as.character(
      x[
        !is.na(x) &
        nzchar(as.character(x))
      ]
    )
  )

  paste(
    sort(x),
    collapse = ";"
  )
}

first_nonmissing <- function(x) {

  x <- x[
    !is.na(x) &
    nzchar(as.character(x))
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  as.character(x[1])
}


# ------------------------------------------------
# 4. Read TargetScan results
# ------------------------------------------------

raw <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "sequence_rank",
  "queried_miRNA",
  "circRNA_seed_2_8",
  "circRNA_weighted_seed_score",
  "target_symbol",
  "target_entrez",
  "target_ensembl",
  "score"
)

missing_cols <- setdiff(
  required_cols,
  colnames(raw)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

raw$score <- safe_numeric(raw$score)

raw <- raw[
  !is.na(raw$target_symbol) &
  nzchar(raw$target_symbol),
  ,
  drop = FALSE
]

message("[ok] raw TargetScan records: ", nrow(raw))
message("[ok] successful miRNAs: ", length(unique(raw$queried_miRNA)))


# ------------------------------------------------
# 5. Collapse to one row per miRNA-gene pair
#
# If the same target gene occurs more than once, retain:
#   - lowest score = most favorable predicted TargetScan interaction
#   - mean score
#   - number of TargetScan records
# ------------------------------------------------

pair_key <- paste(
  raw$queried_miRNA,
  raw$target_symbol,
  sep = "|||"
)

raw_split <- split(
  raw,
  pair_key
)

gene_rows <- vector(
  "list",
  length(raw_split)
)

idx <- 0L

for (d in raw_split) {

  idx <- idx + 1L

  scores <- d$score[
    !is.na(d$score)
  ]

  best_score <- if (length(scores) > 0) {
    min(scores)
  } else {
    NA_real_
  }

  mean_score <- if (length(scores) > 0) {
    mean(scores)
  } else {
    NA_real_
  }

  gene_rows[[idx]] <- data.frame(

    sequence_rank =
      d$sequence_rank[1],

    queried_miRNA =
      d$queried_miRNA[1],

    circRNA_seed_2_8 =
      d$circRNA_seed_2_8[1],

    circRNA_weighted_seed_score =
      d$circRNA_weighted_seed_score[1],

    target_symbol =
      d$target_symbol[1],

    target_entrez =
      first_nonmissing(d$target_entrez),

    target_ensembl =
      first_nonmissing(d$target_ensembl),

    best_TargetScan_score =
      best_score,

    mean_TargetScan_score =
      mean_score,

    n_TargetScan_records =
      nrow(d),

    stringsAsFactors = FALSE
  )
}

gene_level <- do.call(
  rbind,
  gene_rows
)

gene_level$sequence_rank <- safe_numeric(
  gene_level$sequence_rank
)

gene_level <- gene_level[
  order(
    gene_level$sequence_rank,
    gene_level$best_TargetScan_score,
    gene_level$target_symbol
  ),
  ,
  drop = FALSE
]

write.csv(
  gene_level,
  file.path(
    outdir,
    "targetscan_gene_level.csv"
  ),
  row.names = FALSE
)

message("[ok] unique miRNA-gene pairs: ", nrow(gene_level))


# ------------------------------------------------
# 6. Map gene symbols to Entrez IDs
# ------------------------------------------------

symbols <- unique(
  gene_level$target_symbol
)

symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = symbols,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

gene_level$Entrez_ID <- unname(
  symbol_to_entrez[
    gene_level$target_symbol
  ]
)

mapping_rate <- mean(
  !is.na(gene_level$Entrez_ID)
)

message(
  "[ok] Entrez mapping rate: ",
  round(100 * mapping_rate, 1),
  "%"
)

# overwrite file including mapped Entrez IDs
write.csv(
  gene_level,
  file.path(
    outdir,
    "targetscan_gene_level.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------
# 7. Candidate target summary
# ------------------------------------------------

mirnas <- unique(
  gene_level$queried_miRNA[
    order(gene_level$sequence_rank)
  ]
)

summary_rows <- list()

for (i in seq_along(mirnas)) {

  mir <- mirnas[i]

  d <- gene_level[
    gene_level$queried_miRNA == mir,
    ,
    drop = FALSE
  ]

  scores <- d$best_TargetScan_score[
    !is.na(d$best_TargetScan_score)
  ]

  summary_rows[[i]] <- data.frame(

    sequence_rank =
      d$sequence_rank[1],

    miRNA =
      mir,

    circRNA_seed_2_8 =
      d$circRNA_seed_2_8[1],

    circRNA_weighted_seed_score =
      d$circRNA_weighted_seed_score[1],

    n_unique_TargetScan_genes =
      nrow(d),

    n_mapped_Entrez_genes =
      sum(!is.na(d$Entrez_ID)),

    best_score =
      if (length(scores) > 0) min(scores) else NA_real_,

    median_score =
      if (length(scores) > 0) median(scores) else NA_real_,

    n_score_le_minus_0_30 =
      sum(scores <= -0.30),

    n_score_le_minus_0_50 =
      sum(scores <= -0.50),

    stringsAsFactors = FALSE
  )
}

candidate_summary <- do.call(
  rbind,
  summary_rows
)

candidate_summary <- candidate_summary[
  order(candidate_summary$sequence_rank),
  ,
  drop = FALSE
]

write.csv(
  candidate_summary,
  file.path(
    outdir,
    "candidate_target_summary.csv"
  ),
  row.names = FALSE
)

message("[ok] candidate target summary written")


# ------------------------------------------------
# 8. Target-overlap matrices
# ------------------------------------------------

target_sets <- lapply(
  mirnas,
  function(mir) {

    unique(
      gene_level$target_symbol[
        gene_level$queried_miRNA == mir
      ]
    )
  }
)

names(target_sets) <- mirnas


overlap_counts <- matrix(
  0,
  nrow = length(mirnas),
  ncol = length(mirnas),
  dimnames = list(mirnas, mirnas)
)

jaccard <- matrix(
  0,
  nrow = length(mirnas),
  ncol = length(mirnas),
  dimnames = list(mirnas, mirnas)
)


for (i in seq_along(mirnas)) {

  for (j in seq_along(mirnas)) {

    a <- target_sets[[i]]
    b <- target_sets[[j]]

    intersection_n <- length(
      intersect(a, b)
    )

    union_n <- length(
      union(a, b)
    )

    overlap_counts[i, j] <- intersection_n

    jaccard[i, j] <- if (union_n > 0) {
      intersection_n / union_n
    } else {
      NA_real_
    }
  }
}


write.csv(
  data.frame(
    miRNA = rownames(overlap_counts),
    overlap_counts,
    check.names = FALSE
  ),
  file.path(
    outdir,
    "target_overlap_counts.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    miRNA = rownames(jaccard),
    jaccard,
    check.names = FALSE
  ),
  file.path(
    outdir,
    "target_overlap_jaccard.csv"
  ),
  row.names = FALSE
)

message("[ok] target-overlap matrices written")


# ------------------------------------------------
# 9. Build GO Biological Process annotation universe
# ------------------------------------------------

message("[info] building GO Biological Process map...")

go_ids <- AnnotationDbi::mappedkeys(
  org.Hs.egGO2ALLEGS
)

go_meta <- AnnotationDbi::select(
  GO.db,
  keys = go_ids,
  keytype = "GOID",
  columns = c(
    "TERM",
    "ONTOLOGY"
  )
)

go_meta <- go_meta[
  !is.na(go_meta$GOID) &
  !is.na(go_meta$ONTOLOGY) &
  go_meta$ONTOLOGY == "BP",
  ,
  drop = FALSE
]

go_meta <- go_meta[
  !duplicated(go_meta$GOID),
  ,
  drop = FALSE
]

bp_ids <- go_meta$GOID

message("[info] GO-BP terms available: ", length(bp_ids))


bp_maps <- mget(
  bp_ids,
  org.Hs.egGO2ALLEGS,
  ifnotfound = NA
)

names(bp_maps) <- bp_ids

# Convert missing GO mappings to empty vectors,
# otherwise retain unique Entrez IDs as character strings.
bp_maps <- lapply(
  bp_maps,
  function(x) {

    if (
      length(x) == 1 &&
      is.na(x)
    ) {
      return(character(0))
    }

    unique(
      as.character(x)
    )
  }
)

term_name <- setNames(
  go_meta$TERM,
  go_meta$GOID
)

human_bp_background <- unique(
  unlist(
    bp_maps,
    use.names = FALSE
  )
)

targetscan_union_background <- unique(
  gene_level$Entrez_ID[
    !is.na(gene_level$Entrez_ID)
  ]
)

message(
  "[info] human GO-BP background genes: ",
  length(human_bp_background)
)

message(
  "[info] TargetScan-union background genes: ",
  length(targetscan_union_background)
)


# ------------------------------------------------
# 10. Hypergeometric GO enrichment function
# ------------------------------------------------

run_go_enrichment <- function(
  candidate_genes,
  background_genes,
  miRNA_name,
  min_term_size = 10,
  max_term_size = 500,
  min_overlap = 3
) {

  candidate_genes <- unique(
    intersect(
      as.character(candidate_genes),
      as.character(background_genes)
    )
  )

  background_genes <- unique(
    as.character(background_genes)
  )

  n <- length(candidate_genes)
  N <- length(background_genes)

  if (n == 0 || N == 0) {
    return(data.frame())
  }

  rows <- list()
  krow <- 0L


  for (go_id in bp_ids) {

    term_genes <- intersect(
      bp_maps[[go_id]],
      background_genes
    )

    K <- length(term_genes)

    if (
      K < min_term_size ||
      K > max_term_size
    ) {
      next
    }

    overlap <- intersect(
      candidate_genes,
      term_genes
    )

    k <- length(overlap)

    if (k < min_overlap) {
      next
    }

    p <- phyper(
      k - 1,
      K,
      N - K,
      n,
      lower.tail = FALSE
    )

    fold_enrichment <- (
      (k / n) /
      (K / N)
    )

    symbols_overlap <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = overlap,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )

    krow <- krow + 1L

    rows[[krow]] <- data.frame(

      miRNA =
        miRNA_name,

      GO_ID =
        go_id,

      Term =
        unname(term_name[go_id]),

      overlap_genes =
        k,

      candidate_genes =
        n,

      term_genes =
        K,

      background_genes =
        N,

      fold_enrichment =
        fold_enrichment,

      p_value =
        p,

      overlap_symbols =
        collapse_values(
          unname(symbols_overlap)
        ),

      stringsAsFactors = FALSE
    )
  }


  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- do.call(
    rbind,
    rows
  )

  out$FDR <- p.adjust(
    out$p_value,
    method = "BH"
  )

  out <- out[
    order(
      out$FDR,
      out$p_value,
      -out$fold_enrichment
    ),
    ,
    drop = FALSE
  ]

  out
}


# ------------------------------------------------
# 11. Run GO-BP enrichment with TWO backgrounds
#
# A) all human genes annotated to GO-BP
# B) union of genes appearing in our successful TargetScan queries
#
# The second analysis acts as a sensitivity analysis for the fact that
# miRNA-targetable 3' UTR genes are not a random sample of all human genes.
# ------------------------------------------------

human_results <- list()
union_results <- list()

hcounter <- 0L
ucounter <- 0L


for (mir in mirnas) {

  candidate_entrez <- unique(
    gene_level$Entrez_ID[
      gene_level$queried_miRNA == mir &
      !is.na(gene_level$Entrez_ID)
    ]
  )

  message("[GO] ", mir)

  res_human <- run_go_enrichment(
    candidate_genes = candidate_entrez,
    background_genes = human_bp_background,
    miRNA_name = mir
  )

  if (nrow(res_human) > 0) {

    hcounter <- hcounter + 1L
    human_results[[hcounter]] <- res_human
  }


  res_union <- run_go_enrichment(
    candidate_genes = candidate_entrez,
    background_genes = targetscan_union_background,
    miRNA_name = mir
  )

  if (nrow(res_union) > 0) {

    ucounter <- ucounter + 1L
    union_results[[ucounter]] <- res_union
  }
}


go_human <- if (length(human_results) > 0) {
  do.call(rbind, human_results)
} else {
  data.frame()
}

go_union <- if (length(union_results) > 0) {
  do.call(rbind, union_results)
} else {
  data.frame()
}


write.csv(
  go_human,
  file.path(
    outdir,
    "go_bp_enrichment_human_background.csv"
  ),
  row.names = FALSE
)


write.csv(
  go_union,
  file.path(
    outdir,
    "go_bp_enrichment_targetscan_union_background.csv"
  ),
  row.names = FALSE
)


message("[ok] GO enrichment complete")


# ------------------------------------------------
# 12. Top GO terms for compact review
# ------------------------------------------------

top_terms_list <- list()

tcounter <- 0L

for (mir in mirnas) {

  d <- go_human[
    go_human$miRNA == mir,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    next
  }

  # Prefer significant terms; if fewer than 15, still report the best 15
  d <- d[
    order(
      d$FDR,
      d$p_value,
      -d$fold_enrichment
    ),
    ,
    drop = FALSE
  ]

  d <- head(
    d,
    15
  )

  tcounter <- tcounter + 1L

  top_terms_list[[tcounter]] <- d
}


top_terms <- if (length(top_terms_list) > 0) {
  do.call(rbind, top_terms_list)
} else {
  data.frame()
}


write.csv(
  top_terms,
  file.path(
    outdir,
    "go_bp_top_terms.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------
# 13. Pre-specified cancer / tumor-initiation themes
# ------------------------------------------------

themes <- list(

  "Cell cycle / mitosis" =
    paste0(
      "cell cycle|mitotic|mitosis|chromosome segregation|",
      "DNA replication|spindle|G1|G2|S phase"
    ),

  "Apoptosis / cell death" =
    paste0(
      "apopt|programmed cell death|cell death"
    ),

  "DNA damage / repair / checkpoints" =
    paste0(
      "DNA damage|DNA repair|double.strand break|",
      "checkpoint|p53"
    ),

  "Polarity / EMT / migration / adhesion" =
    paste0(
      "cell polarity|epithelial.*mesenchymal|",
      "cell migration|cell motility|cell adhesion|",
      "cell junction|epithelial morphogenesis"
    ),

  "PI3K / AKT / mTOR" =
    paste0(
      "phosphatidylinositol 3.kinase|PI3K|",
      "AKT|mTOR|TOR signaling"
    ),

  "Wnt / beta-catenin" =
    paste0(
      "Wnt|beta.catenin"
    ),

  "Stem / progenitor / differentiation" =
    paste0(
      "stem cell|progenitor|cell differentiation|",
      "mammary|branching morphogenesis"
    ),

  "Proliferation / growth" =
    paste0(
      "cell proliferation|regulation of growth|",
      "cell growth|tissue growth"
    )
)


theme_rows <- list()

theme_counter <- 0L


for (mir in mirnas) {

  d <- go_human[
    go_human$miRNA == mir,
    ,
    drop = FALSE
  ]


  for (theme_name in names(themes)) {

    pattern <- themes[[theme_name]]

    matched <- d[
      grepl(
        pattern,
        d$Term,
        ignore.case = TRUE,
        perl = TRUE
      ),
      ,
      drop = FALSE
    ]


    significant <- matched[
      !is.na(matched$FDR) &
      matched$FDR < 0.05,
      ,
      drop = FALSE
    ]


    if (nrow(matched) > 0) {

      matched <- matched[
        order(
          matched$FDR,
          matched$p_value
        ),
        ,
        drop = FALSE
      ]

      best_term <- matched$Term[1]
      best_FDR <- matched$FDR[1]
      best_fold <- matched$fold_enrichment[1]

      top_sig_terms <- if (nrow(significant) > 0) {

        significant <- significant[
          order(
            significant$FDR,
            significant$p_value
          ),
          ,
          drop = FALSE
        ]

        paste(
          head(
            significant$Term,
            5
          ),
          collapse = " | "
        )

      } else {
        ""
      }

    } else {

      best_term <- NA_character_
      best_FDR <- NA_real_
      best_fold <- NA_real_
      top_sig_terms <- ""
    }


    theme_counter <- theme_counter + 1L

    theme_rows[[theme_counter]] <- data.frame(

      miRNA =
        mir,

      Theme =
        theme_name,

      n_matching_GO_terms =
        nrow(matched),

      n_FDR_lt_0_05 =
        nrow(significant),

      best_GO_term =
        best_term,

      best_FDR =
        best_FDR,

      best_fold_enrichment =
        best_fold,

      top_significant_terms =
        top_sig_terms,

      stringsAsFactors = FALSE
    )
  }
}


theme_summary <- do.call(
  rbind,
  theme_rows
)


write.csv(
  theme_summary,
  file.path(
    outdir,
    "cancer_theme_summary.csv"
  ),
  row.names = FALSE
)


message("[ok] pre-specified cancer-theme summary written")


# ------------------------------------------------
# 14. Finished
# ------------------------------------------------

message("")
message("==============================================")
message("STAGE 4 COMPLETE")
message("==============================================")
message("Output directory: ",
        normalizePath(outdir, winslash = "/", mustWork = FALSE))
message("")
message("Files of greatest interest:")
message("  candidate_target_summary.csv")
message("  go_bp_top_terms.csv")
message("  cancer_theme_summary.csv")
message("  target_overlap_jaccard.csv")
message("")
message("Interpretation rule:")
message("Do NOT choose a final miRNA from pathway enrichment alone.")
message("Stage 4 is one evidence layer to integrate with:")
message("  - circRNA seed architecture")
message("  - breast-cancer literature")
message("  - TargetScan prediction")
message("  - later Kaplan-Meier evidence")
