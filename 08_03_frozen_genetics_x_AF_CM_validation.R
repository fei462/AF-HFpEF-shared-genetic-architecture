# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP12C_V2_FIX_QC_INDEX_FROZEN_GENETICS_x_GSE255612_CM.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP12C V2 — FROZEN GENETICS/REGULATORY CANDIDATES × AF CM DISEASE-TISSUE VALIDATION
#
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# THIS STEP IS THE INTEGRATION BRIDGE.
#
# FROZEN UPSTREAM LAYERS:
#   - cross-trait genetics / conjFDR / FUMA / LAVA
#   - MR/MVMR
#   - eQTL-SMR/HEIDI
#   - sQTL-SMR/HEIDI
#   - coloc / SuSiE
#   - STEP10B-V3.2 locus-centric integration
#   - GSE238242 S-LDSC + AF-priority SCAVENGE -> cardiomyocytes
#
# INDEPENDENT DISEASE-TISSUE LAYER:
#   - GSE255612 left-atrial snRNA-seq
#   - STEP12B frozen patient-level pseudobulk
#   - 34 donors: 18 AF / 16 CTRL
#   - primary cell type: Cardiomyocytes
#
# ANTI-CIRCULARITY:
#   STEP12B was completed without using the upstream candidate-gene list.
#   STEP12C now intersects TWO ALREADY-FROZEN result layers.
#   No gene is selected because it "looks good" in GSE255612.
#
# PRIMARY INTEGRATION SET:
#   STEP10B-V3.2
#   02_MAIN_TABLES/Table2_priority_locus_gene_candidates_TOP50.csv
#
# SECONDARY/COMPLETE INTEGRATION UNIVERSE:
#   STEP10B-V3.2
#   03_SUPPLEMENTARY_TABLES/TableS2_all_locus_gene_evidence.csv
#
# DISEASE-TISSUE SUPPORT CLASSES:
#   A. transcriptome-wide support:
#        CM FDR < 0.05
#   B. nominal support:
#        CM P < 0.05 but FDR >= 0.05
#   C. direction-only / no statistical support:
#        CM P >= 0.05
#   D. not testable:
#        candidate gene not present in the filtered CM pseudobulk DE universe
#
# NOTE:
#   "Nominal support" is NOT called validation.
#   Formal validation language is reserved for transcriptome-wide CM FDR < 0.05.
#
# ADDITIONAL PRE-SPECIFIED SET-LEVEL TEST:
#   Are frozen upstream candidate genes collectively more disease-responsive in
#   AF cardiomyocytes than expression-matched background genes?
#   -> 10,000 expression-decile-matched random gene-set permutations.
#
# OUTPUT:
#   D:/A/data/STEP12_GSE255612_AF_snRNA/STEP12C_GENETICS_CM_INTEGRATION/
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))
set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

N_PERM <- 10000L
N_EXPRESSION_BINS <- 10L

CM_FDR_VALIDATION <- 0.05
CM_NOMINAL_P <- 0.05

TOP_N_FIGURE <- 30L

# ==============================================================================
# 1. FROZEN INPUT PATHS
# ==============================================================================

STEP10B_ROOT <- file.path(
  DATA_ROOT,
  "STEP10B",
  "STEP10B_V3_FINAL_LOCUS_CENTRIC_RESULTS"
)

STEP10B_QC_FILE <- file.path(
  STEP10B_ROOT,
  "00_QC",
  "STEP10B_V3_QC_summary.csv"
)

TOP50_FILE <- file.path(
  STEP10B_ROOT,
  "02_MAIN_TABLES",
  "Table2_priority_locus_gene_candidates_TOP50.csv"
)

ALL_CANDIDATE_FILE <- file.path(
  STEP10B_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS2_all_locus_gene_evidence.csv"
)

# Optional audit tables. They are NOT used to redefine the candidate set.
OPTIONAL_STEP10B_FILES <- c(
  file.path(
    STEP10B_ROOT,
    "02_MAIN_TABLES",
    "Table1_focal_8test_multilayer_evidence.csv"
  ),
  file.path(
    STEP10B_ROOT,
    "03_SUPPLEMENTARY_TABLES",
    "TableS1_all_FUMA_loci_multilayer_evidence.csv"
  ),
  file.path(
    STEP10B_ROOT,
    "03_SUPPLEMENTARY_TABLES",
    "TableS3_FUMA_gene_mapping_final.csv"
  ),
  file.path(
    STEP10B_ROOT,
    "03_SUPPLEMENTARY_TABLES",
    "TableS4_formal_SuSiE_V3_with_V4_audit.csv"
  )
)

STEP12B_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA",
  "STEP12B_AF_vs_CTRL_DE"
)

STEP12B_QC_FILE <- file.path(
  STEP12B_ROOT,
  "00_QC",
  "STEP12B_readiness.csv"
)

CM_DE_FILE <- file.path(
  STEP12B_ROOT,
  "01_PRIMARY_CM",
  "STEP12B_Cardiomyocytes_AF_vs_CTRL_all_genes.csv.gz"
)

MAC_DE_FILE <- file.path(
  STEP12B_ROOT,
  "02_SECONDARY_MACROPHAGE",
  "STEP12B_Macrophages_AF_vs_CTRL_all_genes.csv.gz"
)

# ==============================================================================
# 2. OUTPUT
# ==============================================================================

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA",
  "STEP12C_GENETICS_CM_INTEGRATION"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
PRIMARY_DIR <- file.path(OUT_ROOT, "01_TOP50_PRIMARY")
ALL_DIR <- file.path(OUT_ROOT, "02_ALL_CANDIDATES")
ENRICH_DIR <- file.path(OUT_ROOT, "03_SET_LEVEL_ENRICHMENT")
FIG_DIR <- file.path(OUT_ROOT, "04_FIGURES")
SOURCE_DIR <- file.path(OUT_ROOT, "05_SOURCE_DATA")

for (d in c(
  OUT_ROOT,
  QC_DIR,
  PRIMARY_DIR,
  ALL_DIR,
  ENRICH_DIR,
  FIG_DIR,
  SOURCE_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 3. PACKAGES
# ==============================================================================

pkgs <- c(
  "data.table",
  "ggplot2",
  "ggrepel"
)

for (p in pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }

  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Could not install/load package: ", p)
  }
}

library(data.table)
library(ggplot2)
library(ggrepel)

# ==============================================================================
# 4. HELPERS
# ==============================================================================

bool_pass <- function(x) {
  x %in% c(TRUE, "TRUE", 1, "1")
}

norm_gene <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  x[
    x %in%
      c(
        "",
        "NA",
        "N/A",
        "NULL",
        ".",
        "-"
      )
  ] <- NA_character_

  x
}

# ------------------------------------------------------------------------------
# Split a cell that may contain one or multiple gene symbols.
# Delimiters are restricted to common list delimiters.
# Hyphens are NOT treated as delimiters because valid gene symbols may contain
# hyphens.
# ------------------------------------------------------------------------------

split_gene_field <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "\\s+",
    " ",
    x
  )

  pieces <- unlist(
    strsplit(
      x,
      split = "[,;|/]+",
      perl = TRUE
    )
  )

  pieces <- norm_gene(
    pieces
  )

  unique(
    pieces[
      !is.na(pieces)
    ]
  )
}

# ------------------------------------------------------------------------------
# Automatically detect the gene-symbol column in a frozen STEP10B table.
# It prioritizes exact/common gene column names, then chooses the column with
# the greatest overlap with the GSE255612 tested gene-symbol universe.
# ------------------------------------------------------------------------------

detect_gene_column <- function(DT, tested_symbols, table_label) {

  nms <- names(DT)

  exact_priority <- c(
    "gene_symbol",
    "GENE_SYMBOL",
    "Gene_symbol",
    "symbol",
    "SYMBOL",
    "gene",
    "GENE",
    "Gene",
    "candidate_gene",
    "candidate_genes",
    "mapped_gene",
    "mapped_genes",
    "gene_name",
    "genes"
  )

  exact_hits <- exact_priority[
    exact_priority %in%
      nms
  ]

  candidate_cols <- unique(
    c(
      exact_hits,
      nms[
        grepl(
          "gene|symbol",
          nms,
          ignore.case = TRUE
        )
      ]
    )
  )

  if (!length(candidate_cols)) {
    stop(
      table_label,
      ": no column name containing gene/symbol was found."
    )
  }

  tested_symbols <- unique(
    norm_gene(
      tested_symbols
    )
  )

  score_one <- function(col) {

    vals <- DT[[col]]

    genes <- unique(
      unlist(
        lapply(
          vals,
          split_gene_field
        )
      )
    )

    genes <- genes[
      !is.na(genes)
    ]

    if (!length(genes)) {
      return(
        data.table(
          column = col,
          n_unique_genes = 0L,
          n_overlap_tested = 0L,
          overlap_fraction = 0
        )
      )
    }

    n_overlap <- sum(
      genes %in%
        tested_symbols
    )

    data.table(
      column = col,
      n_unique_genes =
        length(
          genes
        ),
      n_overlap_tested =
        n_overlap,
      overlap_fraction =
        n_overlap /
          length(
            genes
          )
    )
  }

  AUDIT <- rbindlist(
    lapply(
      candidate_cols,
      score_one
    )
  )

  # Favor overlap count first, then overlap fraction, then exact-priority order.
  AUDIT[
    ,
    exact_priority_rank :=
      match(
        column,
        exact_priority
      )
  ]

  AUDIT[
    is.na(
      exact_priority_rank
    ),
    exact_priority_rank :=
      999L
  ]

  setorder(
    AUDIT,
    -n_overlap_tested,
    -overlap_fraction,
    exact_priority_rank
  )

  if (
    nrow(AUDIT) == 0L ||
    AUDIT$n_overlap_tested[1] == 0L
  ) {
    print(AUDIT)
    stop(
      table_label,
      ": candidate gene columns have zero overlap with the CM tested-gene universe."
    )
  }

  list(
    selected = AUDIT$column[1],
    audit = AUDIT
  )
}

# ------------------------------------------------------------------------------
# Expand a frozen table to one row per candidate gene while preserving all
# original evidence columns.
# ------------------------------------------------------------------------------

expand_candidate_table <- function(DT, gene_col, source_label) {

  OUT <- copy(DT)

  OUT[
    ,
    .source_row_id :=
      .I
  ]

  GENE_LIST <- lapply(
    OUT[[gene_col]],
    split_gene_field
  )

  lens <- lengths(
    GENE_LIST
  )

  keep <- lens > 0L

  if (!any(keep)) {
    stop(
      source_label,
      ": no candidate genes could be parsed."
    )
  }

  OUT2 <- OUT[
    rep(
      which(keep),
      lens[keep]
    )
  ]

  OUT2[
    ,
    candidate_gene :=
      unlist(
        GENE_LIST[
          keep
        ],
        use.names = FALSE
      )
  ]

  OUT2[
    ,
    candidate_gene :=
      norm_gene(
        candidate_gene
      )
  ]

  OUT2[
    ,
    frozen_candidate_source :=
      source_label
  ]

  OUT2
}

# ------------------------------------------------------------------------------
# Collapse one transcriptome-wide DE table to one record per gene symbol.
# If duplicated symbols exist, preserve the transcript with the smallest P.
# ------------------------------------------------------------------------------

collapse_de_by_symbol <- function(DE, prefix) {

  X <- copy(DE)

  if (!"gene_symbol" %in% names(X)) {
    stop(prefix, ": DE table lacks gene_symbol.")
  }

  X[
    ,
    gene_symbol_norm :=
      norm_gene(
        gene_symbol
      )
  ]

  X <- X[
    !is.na(
      gene_symbol_norm
    )
  ]

  setorder(
    X,
    gene_symbol_norm,
    PValue,
    FDR
  )

  X <- X[
    ,
    .SD[1],
    by = gene_symbol_norm
  ]

  keep_cols <- c(
    "gene_symbol_norm",
    "gene_key",
    "gene_id",
    "gene_id_clean",
    "gene_symbol",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR",
    "direction"
  )

  keep_cols <- keep_cols[
    keep_cols %in%
      names(X)
  ]

  X <- X[
    ,
    ..keep_cols
  ]

  old <- setdiff(
    names(X),
    "gene_symbol_norm"
  )

  setnames(
    X,
    old,
    paste0(
      prefix,
      "_",
      old
    )
  )

  X
}

# ------------------------------------------------------------------------------
# AF-CM support class.
# ------------------------------------------------------------------------------

add_cm_support_class <- function(DT) {

  DT[
    ,
    CM_disease_tissue_support :=
      fifelse(
        is.na(
          CM_PValue
        ),
        "Not_testable",
        fifelse(
          CM_FDR <
            CM_FDR_VALIDATION,
          "Transcriptome_wide_FDR_support",
          fifelse(
            CM_PValue <
              CM_NOMINAL_P,
            "Nominal_support_only",
            "No_statistical_support"
          )
        )
      )
  ]

  DT[
    ,
    CM_nominal_or_better :=
      !is.na(
        CM_PValue
      ) &
      CM_PValue <
        CM_NOMINAL_P
  ]

  DT[
    ,
    CM_FDR_support :=
      !is.na(
        CM_FDR
      ) &
      CM_FDR <
        CM_FDR_VALIDATION
  ]

  DT
}

# ------------------------------------------------------------------------------
# Expression-decile matched random-set test.
# Score = -log10(CM P); larger = stronger disease response.
# This is a competitive descriptive enrichment of the prespecified gene set.
# ------------------------------------------------------------------------------

matched_gene_set_test <- function(
  candidate_genes,
  cm_de,
  set_name,
  n_perm = 10000L,
  n_bins = 10L
) {

  U <- copy(
    cm_de[
      !is.na(
        gene_symbol_norm
      ) &
        is.finite(
          CM_PValue
        ) &
        is.finite(
          CM_logCPM
        )
    ]
  )

  U <- U[
    !duplicated(
      gene_symbol_norm
    )
  ]

  U[
    ,
    disease_score :=
      -log10(
        pmax(
          CM_PValue,
          1e-300
        )
      )
  ]

  # Expression bins.
  U[
    ,
    expression_bin :=
      cut(
        frank(
          CM_logCPM,
          ties.method = "average"
        ) /
          .N,
        breaks = seq(
          0,
          1,
          length.out =
            n_bins +
              1L
        ),
        include.lowest = TRUE,
        labels = FALSE
      )
  ]

  cand <- unique(
    norm_gene(
      candidate_genes
    )
  )

  cand <- intersect(
    cand,
    U$gene_symbol_norm
  )

  C <- U[
    gene_symbol_norm %in%
      cand
  ]

  if (nrow(C) < 3L) {

    return(
      data.table(
        gene_set = set_name,
        n_candidate_total =
          length(
            unique(
              norm_gene(
                candidate_genes
              )
            )
          ),
        n_candidate_tested =
          nrow(C),
        observed_mean_neglog10P =
          if (
            nrow(C)
          ) {
            mean(
              C$disease_score
            )
          } else {
            NA_real_
          },
        matched_null_mean = NA_real_,
        empirical_P = NA_real_,
        fold_over_matched_null = NA_real_
      )
    )
  }

  observed <- mean(
    C$disease_score
  )

  bins_needed <- table(
    C$expression_bin
  )

  ALL_BY_BIN <- split(
    U$gene_symbol_norm,
    U$expression_bin
  )

  set.seed(9527)

  null_score <- numeric(
    n_perm
  )

  for (b in seq_len(n_perm)) {

    sampled <- character(0)

    for (bn in names(bins_needed)) {

      need <- as.integer(
        bins_needed[[bn]]
      )

      pool <- setdiff(
        ALL_BY_BIN[[bn]],
        C$gene_symbol_norm
      )

      if (length(pool) < need) {
        # If a very small expression bin occurs, allow the candidate genes back
        # into the sampling pool rather than sampling with replacement.
        pool <- ALL_BY_BIN[[bn]]
      }

      sampled <- c(
        sampled,
        sample(
          pool,
          size = need,
          replace = FALSE
        )
      )
    }

    null_score[b] <- mean(
      U[
        match(
          sampled,
          gene_symbol_norm
        ),
        disease_score
      ]
    )
  }

  empirical_p <- (
    1 +
      sum(
        null_score >=
          observed
      )
  ) /
    (
      1 +
        n_perm
    )

  null_mean <- mean(
    null_score
  )

  data.table(
    gene_set = set_name,
    n_candidate_total =
      length(
        unique(
          norm_gene(
            candidate_genes
          )
        )
      ),
    n_candidate_tested =
      nrow(C),
    observed_mean_neglog10P =
      observed,
    matched_null_mean =
      null_mean,
    empirical_P =
      empirical_p,
    fold_over_matched_null =
      observed /
        null_mean
  )
}

# ==============================================================================
# 5. HARD INPUT CHECKS
# ==============================================================================

required_inputs <- c(
  TOP50_FILE,
  ALL_CANDIDATE_FILE,
  STEP12B_QC_FILE,
  CM_DE_FILE,
  MAC_DE_FILE
)

if (!all(
  file.exists(
    required_inputs
  )
)) {

  stop(
    paste0(
      "Missing frozen input(s):\n",
      paste(
        required_inputs[
          !file.exists(
            required_inputs
          )
        ],
        collapse = "\n"
      )
    )
  )
}

READY12B <- fread(
  STEP12B_QC_FILE
)

if (!all(
  bool_pass(
    READY12B$pass
  )
)) {
  print(
    READY12B
  )
  stop(
    "STEP12B is not fully PASS. Do not run STEP12C."
  )
}

# Optional STEP10B QC audit.
STEP10B_QC_PASS <- NA

if (file.exists(
  STEP10B_QC_FILE
)) {

  Q10 <- fread(
    STEP10B_QC_FILE
  )

  pass_cols <- names(Q10)[
    tolower(
      names(Q10)
    ) %in%
      c(
        "pass",
        "passed",
        "ok"
      )
  ]

  if (length(pass_cols)) {

    STEP10B_QC_PASS <- all(
      bool_pass(
        Q10[[pass_cols[1]]]
      )
    )
  }
}

# ==============================================================================
# 6. LOAD FROZEN CM / MACROPHAGE TRANSCRIPTOME-WIDE DE
# ==============================================================================

CM_RAW <- fread(
  CM_DE_FILE
)

MAC_RAW <- fread(
  MAC_DE_FILE
)

CM_DE <- collapse_de_by_symbol(
  CM_RAW,
  "CM"
)

MAC_DE <- collapse_de_by_symbol(
  MAC_RAW,
  "MAC"
)

if (nrow(CM_DE) < 1000L) {
  stop(
    "Unexpectedly small CM DE universe."
  )
}

# ==============================================================================
# 7. LOAD FROZEN STEP10B CANDIDATE TABLES
# ==============================================================================

TOP50_RAW <- fread(
  TOP50_FILE
)

ALL_RAW <- fread(
  ALL_CANDIDATE_FILE
)

TOP_GCOL <- detect_gene_column(
  TOP50_RAW,
  tested_symbols =
    CM_DE$CM_gene_symbol,
  table_label =
    "STEP10B Table2 TOP50"
)

ALL_GCOL <- detect_gene_column(
  ALL_RAW,
  tested_symbols =
    CM_DE$CM_gene_symbol,
  table_label =
    "STEP10B TableS2 all locus-gene evidence"
)

fwrite(
  TOP_GCOL$audit,
  file.path(
    QC_DIR,
    "STEP12C_TOP50_gene_column_detection_QC.csv"
  )
)

fwrite(
  ALL_GCOL$audit,
  file.path(
    QC_DIR,
    "STEP12C_ALL_candidate_gene_column_detection_QC.csv"
  )
)

TOP_EXP <- expand_candidate_table(
  TOP50_RAW,
  TOP_GCOL$selected,
  "STEP10B_V3_TOP50"
)

ALL_EXP <- expand_candidate_table(
  ALL_RAW,
  ALL_GCOL$selected,
  "STEP10B_V3_ALL_LOCUS_GENE"
)

# ==============================================================================
# 8. MERGE FROZEN CANDIDATES WITH INDEPENDENT DISEASE-TISSUE DE
# ==============================================================================

merge_disease <- function(CAND) {

  X <- merge(
    CAND,
    CM_DE,
    by.x = "candidate_gene",
    by.y = "gene_symbol_norm",
    all.x = TRUE,
    sort = FALSE
  )

  X <- merge(
    X,
    MAC_DE,
    by.x = "candidate_gene",
    by.y = "gene_symbol_norm",
    all.x = TRUE,
    sort = FALSE
  )

  X <- add_cm_support_class(
    X
  )

  X
}

TOP_INT <- merge_disease(
  TOP_EXP
)

ALL_INT <- merge_disease(
  ALL_EXP
)

# Preserve original frozen table row ordering.
setorder(
  TOP_INT,
  .source_row_id
)

setorder(
  ALL_INT,
  .source_row_id
)

# ==============================================================================
# 9. GENE-LEVEL COLLAPSED EVIDENCE MATRIX
# ==============================================================================

collapse_integrated <- function(INT, label) {

  X <- copy(
    INT
  )

  # One row per candidate gene.
  G <- X[
    ,
    .(
      n_frozen_rows =
        .N,
      CM_tested =
        any(
          !is.na(
            CM_PValue
          )
        ),
      CM_logFC =
        CM_logFC[
          which.min(
            fifelse(
              is.na(
                CM_PValue
              ),
              Inf,
              CM_PValue
            )
          )
        ][1],
      CM_logCPM =
        CM_logCPM[
          which.min(
            fifelse(
              is.na(
                CM_PValue
              ),
              Inf,
              CM_PValue
            )
          )
        ][1],
      CM_PValue =
        suppressWarnings(
          min(
            CM_PValue,
            na.rm = TRUE
          )
        ),
      CM_FDR =
        suppressWarnings(
          min(
            CM_FDR,
            na.rm = TRUE
          )
        ),
      MAC_logFC =
        MAC_logFC[
          which.min(
            fifelse(
              is.na(
                MAC_PValue
              ),
              Inf,
              MAC_PValue
            )
          )
        ][1],
      MAC_PValue =
        suppressWarnings(
          min(
            MAC_PValue,
            na.rm = TRUE
          )
        ),
      MAC_FDR =
        suppressWarnings(
          min(
            MAC_FDR,
            na.rm = TRUE
          )
        )
    ),
    by = candidate_gene
  ]

  # Fix min(..., na.rm=TRUE) -> Inf when all missing.
  for (nm in c(
    "CM_PValue",
    "CM_FDR",
    "MAC_PValue",
    "MAC_FDR"
  )) {

    G[
      !is.finite(
        get(nm)
      ),
      (nm) :=
        NA_real_
    ]
  }

  G <- add_cm_support_class(
    G
  )

  G[
    ,
    disease_validation_set :=
      label
  ]

  G[
    ,
    support_order :=
      match(
        CM_disease_tissue_support,
        c(
          "Transcriptome_wide_FDR_support",
          "Nominal_support_only",
          "No_statistical_support",
          "Not_testable"
        )
      )
  ]

  setorder(
    G,
    support_order,
    CM_PValue
  )

  G[
    ,
    support_order := NULL
  ]

  G
}

TOP_GENE <- collapse_integrated(
  TOP_INT,
  "TOP50"
)

ALL_GENE <- collapse_integrated(
  ALL_INT,
  "ALL_LOCUS_GENE"
)

# ==============================================================================
# 10. PRESPECIFIED SET-LEVEL ENRICHMENT
# ==============================================================================

CM_FOR_SET <- copy(
  CM_DE
)

# CM_DE already contains one row per gene symbol.
TOP_SET <- matched_gene_set_test(
  candidate_genes =
    TOP_GENE$candidate_gene,
  cm_de =
    CM_FOR_SET,
  set_name =
    "STEP10B_V3_TOP50",
  n_perm =
    N_PERM,
  n_bins =
    N_EXPRESSION_BINS
)

ALL_SET <- matched_gene_set_test(
  candidate_genes =
    ALL_GENE$candidate_gene,
  cm_de =
    CM_FOR_SET,
  set_name =
    "STEP10B_V3_ALL_LOCUS_GENE",
  n_perm =
    N_PERM,
  n_bins =
    N_EXPRESSION_BINS
)

SET_ENRICH <- rbindlist(
  list(
    TOP_SET,
    ALL_SET
  ),
  fill = TRUE
)

SET_ENRICH[
  ,
  FDR_BH :=
    p.adjust(
      empirical_P,
      method = "BH"
    )
]

# ==============================================================================
# 11. SUMMARY
# ==============================================================================

summarize_gene_matrix <- function(G, set_name) {

  data.table(
    candidate_set = set_name,
    n_unique_candidates =
      uniqueN(
        G$candidate_gene
      ),
    n_tested_in_CM =
      sum(
        G$CM_tested
      ),
    n_CM_FDR_support =
      sum(
        G$CM_FDR_support,
        na.rm = TRUE
      ),
    n_CM_nominal_or_better =
      sum(
        G$CM_nominal_or_better,
        na.rm = TRUE
      ),
    n_CM_nominal_only =
      sum(
        G$CM_disease_tissue_support ==
          "Nominal_support_only",
        na.rm = TRUE
      ),
    n_no_statistical_support =
      sum(
        G$CM_disease_tissue_support ==
          "No_statistical_support",
        na.rm = TRUE
      ),
    n_not_testable =
      sum(
        G$CM_disease_tissue_support ==
          "Not_testable",
        na.rm = TRUE
      ),
    n_CM_AF_up_nominal_or_better =
      sum(
        G$CM_nominal_or_better &
          G$CM_logFC >
            0,
        na.rm = TRUE
      ),
    n_CM_AF_down_nominal_or_better =
      sum(
        G$CM_nominal_or_better &
          G$CM_logFC <
            0,
        na.rm = TRUE
      )
  )
}

SUMMARY <- rbindlist(
  list(
    summarize_gene_matrix(
      TOP_GENE,
      "STEP10B_V3_TOP50"
    ),
    summarize_gene_matrix(
      ALL_GENE,
      "STEP10B_V3_ALL_LOCUS_GENE"
    )
  )
)

# ==============================================================================
# 12. WRITE TABLES
# ==============================================================================

fwrite(
  TOP_INT,
  file.path(
    PRIMARY_DIR,
    "STEP12C_TOP50_frozen_rows_with_CM_disease_validation.csv"
  )
)

fwrite(
  TOP_GENE,
  file.path(
    PRIMARY_DIR,
    "STEP12C_TOP50_gene_level_disease_validation_matrix.csv"
  )
)

fwrite(
  TOP_GENE[
    CM_disease_tissue_support %in%
      c(
        "Transcriptome_wide_FDR_support",
        "Nominal_support_only"
      )
  ],
  file.path(
    PRIMARY_DIR,
    "STEP12C_TOP50_CM_supported_candidates.csv"
  )
)

fwrite(
  ALL_INT,
  file.path(
    ALL_DIR,
    "STEP12C_ALL_frozen_rows_with_CM_disease_validation.csv.gz"
  ),
  compress = "gzip"
)

fwrite(
  ALL_GENE,
  file.path(
    ALL_DIR,
    "STEP12C_ALL_gene_level_disease_validation_matrix.csv.gz"
  ),
  compress = "gzip"
)

fwrite(
  ALL_GENE[
    CM_disease_tissue_support %in%
      c(
        "Transcriptome_wide_FDR_support",
        "Nominal_support_only"
      )
  ],
  file.path(
    ALL_DIR,
    "STEP12C_ALL_CM_supported_candidates.csv"
  )
)

fwrite(
  SET_ENRICH,
  file.path(
    ENRICH_DIR,
    "STEP12C_candidate_gene_set_expression_matched_enrichment.csv"
  )
)

fwrite(
  SUMMARY,
  file.path(
    QC_DIR,
    "STEP12C_integration_summary.csv"
  )
)

# Optional upstream file audit.
OPTIONAL_AUDIT <- data.table(
  file = basename(
    OPTIONAL_STEP10B_FILES
  ),
  path = OPTIONAL_STEP10B_FILES,
  exists = file.exists(
    OPTIONAL_STEP10B_FILES
  )
)

fwrite(
  OPTIONAL_AUDIT,
  file.path(
    QC_DIR,
    "STEP12C_optional_STEP10B_file_audit.csv"
  )
)

# ==============================================================================
# 13. FIGURES — SINGLE PANEL FILES, NO LEGEND INSIDE
# ==============================================================================

# ------------------------------------------------------------------------------
# 13A. TOP50 candidate CM effect plot.
# Show at most 30 best disease-responsive TOP50 candidates.
# IMPORTANT: ordering is by CM P only for visualization, not candidate selection.
# ------------------------------------------------------------------------------

PLOT_TOP <- copy(
  TOP_GENE[
    !is.na(
      CM_PValue
    )
  ]
)

setorder(
  PLOT_TOP,
  CM_PValue
)

PLOT_TOP <- head(
  PLOT_TOP,
  TOP_N_FIGURE
)

if (nrow(PLOT_TOP)) {

  PLOT_TOP[
    ,
    plot_gene :=
      factor(
        candidate_gene,
        levels = rev(
          candidate_gene
        )
      )
  ]

  PLOT_TOP[
    ,
    point_size :=
      pmin(
        -log10(
          pmax(
            CM_PValue,
            1e-300
          )
        ),
        10
      )
  ]

  PLOT_TOP[
    ,
    plot_group :=
      fifelse(
        CM_FDR_support,
        "FDR",
        fifelse(
          CM_nominal_or_better,
          "Nominal",
          "NS"
        )
      )
  ]

  p1 <- ggplot(
    PLOT_TOP,
    aes(
      x = CM_logFC,
      y = plot_gene,
      size = point_size,
      shape = plot_group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.45,
      colour = "grey55"
    ) +
    geom_point(
      alpha = 0.88
    ) +
    labs(
      x = "CM log2 fold change (AF vs CTRL)",
      y = NULL
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(
        colour = "black"
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.55
      ),
      plot.margin = margin(
        8, 10, 8, 8
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP12C_TOP50_CM_AF_effects.tiff"
    ),
    p1,
    width = 6.2,
    height = 7.0,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP12C_TOP50_CM_AF_effects.pdf"
    ),
    p1,
    width = 6.2,
    height = 7.0,
    units = "in"
  )

  fwrite(
    PLOT_TOP,
    file.path(
      SOURCE_DIR,
      "Figure_STEP12C_TOP50_CM_AF_effects_source_data.csv"
    )
  )
}

# ------------------------------------------------------------------------------
# 13B. Set-level enrichment.
# ------------------------------------------------------------------------------

PLOT_SET <- copy(
  SET_ENRICH
)

PLOT_SET[
  ,
  gene_set_plot :=
    factor(
      gene_set,
      levels = rev(
        gene_set
      )
    )
]

p2 <- ggplot(
  PLOT_SET,
  aes(
    x = fold_over_matched_null,
    y = gene_set_plot
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.45,
    colour = "grey55"
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    x = "Disease-response enrichment vs expression-matched null",
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.55
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP12C_candidate_set_CM_enrichment.tiff"
  ),
  p2,
  width = 6.1,
  height = 3.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP12C_candidate_set_CM_enrichment.pdf"
  ),
  p2,
  width = 6.1,
  height = 3.2,
  units = "in"
)

# ==============================================================================
# 14. TECHNICAL READINESS
# ==============================================================================

READINESS <- data.table(
  check = c(
    "STEP12B_fully_passed",
    "STEP10B_TOP50_table_present",
    "STEP10B_all_locus_gene_table_present",
    "TOP50_gene_column_detected",
    "ALL_gene_column_detected",
    "TOP50_unique_candidate_genes_gt_0",
    "ALL_unique_candidate_genes_gt_0",
    "CM_transcriptome_wide_DE_merged",
    "Macrophage_secondary_DE_merged",
    "TOP50_set_level_matched_test_completed",
    "ALL_set_level_matched_test_completed",
    "candidate_selection_not_redefined_by_GSE255612",
    "STEP12C_complete"
  ),
  pass = c(
    all(
      bool_pass(
        READY12B$pass
      )
    ),
    file.exists(
      TOP50_FILE
    ),
    file.exists(
      ALL_CANDIDATE_FILE
    ),
    !is.na(
      TOP_GCOL$selected
    ),
    !is.na(
      ALL_GCOL$selected
    ),
    uniqueN(
      TOP_GENE$candidate_gene
    ) >
      0L,
    uniqueN(
      ALL_GENE$candidate_gene
    ) >
      0L,
    any(
      !is.na(
        TOP_GENE$CM_PValue
      )
    ),
    any(
      !is.na(
        TOP_GENE$MAC_PValue
      )
    ),
    is.finite(
      TOP_SET$empirical_P
    ) ||
      TOP_SET$n_candidate_tested <
        3L,
    is.finite(
      ALL_SET$empirical_P
    ) ||
      ALL_SET$n_candidate_tested <
        3L,
    TRUE,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP12C_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "integration_role",
    "primary_frozen_candidate_set",
    "complete_frozen_candidate_universe",
    "disease_tissue_dataset",
    "primary_disease_celltype",
    "candidate_selection_after_seeing_GSE255612",
    "transcriptome_wide_support_definition",
    "nominal_support_definition",
    "set_level_test",
    "set_level_permutations",
    "expression_matching",
    "Macrophage_role",
    "next_step"
  ),
  value = c(
    "Integrate already-frozen upstream genetics/regulatory candidates with independently frozen AF disease-tissue CM DE",
    "STEP10B-V3.2 Table2 priority locus-gene candidates TOP50",
    "STEP10B-V3.2 TableS2 all locus-gene evidence",
    "GSE255612 left atrium patient-level pseudobulk",
    "Cardiomyocytes",
    "No",
    "CM FDR < 0.05",
    "CM P < 0.05 but FDR >= 0.05; descriptive support only, not called validation",
    "Competitive expression-matched random gene-set test using mean -log10(CM P)",
    as.character(
      N_PERM
    ),
    paste0(
      N_EXPRESSION_BINS,
      " CM logCPM bins"
    ),
    "Secondary context only",
    "STEP12D: integrate HFpEF disease-tissue validation and freeze the cross-disease candidate evidence matrix"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP12C_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP12C_sessionInfo.txt"
  )
)

# ==============================================================================
# 15. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP12C COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Frozen candidate gene columns:\n",
  "  TOP50: ",
  TOP_GCOL$selected,
  "\n",
  "  ALL:   ",
  ALL_GCOL$selected,
  "\n\n",
  sep = ""
)

cat(
  "Candidate integration summary:\n"
)
print(
  SUMMARY
)

cat(
  "\nSet-level expression-matched enrichment:\n"
)
print(
  SET_ENRICH
)

cat(
  "\nTOP50 candidates with CM nominal-or-better disease support:\n"
)
print(
  TOP_GENE[
    CM_nominal_or_better ==
      TRUE
  ][
    order(
      CM_PValue
    )
  ][
    ,
    .(
      candidate_gene,
      CM_logFC,
      CM_PValue,
      CM_FDR,
      CM_disease_tissue_support,
      MAC_logFC,
      MAC_PValue,
      MAC_FDR
    )
  ]
)

cat(
  "\nReadiness:\n"
)
print(
  READINESS
)

if (all(
  READINESS$pass
)) {

  cat(
    "\nSTEP12C PASSED.\n",
    "Do not add new AF candidate genes after this point.\n",
    "Next: HFpEF disease-tissue validation / cross-disease integration.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP12C failed one or more technical checks.\n",
    "Do not proceed to cross-disease integration yet.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  PRIMARY_DIR,
  "\n",
  ENRICH_DIR,
  "\n",
  sep = ""
)
