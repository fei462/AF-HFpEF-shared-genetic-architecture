# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP13C_V2_FIX_PERMUTATION_NULL_AF_FDR.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP13C — FINAL CROSS-DISEASE INTEGRATION
# AF-CM × HFpEF-CM × FROZEN GENETICS / REGULATORY CANDIDATES
#
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# ROLE
#   This is the LAST core statistical integration block.
#
# FROZEN INPUT LAYERS
#   1) STEP10B-V3.2 upstream genetics/regulatory candidates
#   2) STEP12B/12C GSE255612 AF left-atrial cardiomyocyte disease-tissue result
#   3) STEP13B SCP3342 HFpEF cardiomyocyte disease-tissue result
#
# ANTI-CIRCULARITY
#   - Candidate genes were frozen before disease-tissue integration.
#   - AF-CM DE was frozen without using the candidate list.
#   - HFpEF-CM DE was frozen without using the candidate list.
#   - STEP13C only INTERSECTS already-frozen layers.
#   - No new candidate is created because it "looks good" in AF or HFpEF.
#
# PRIMARY QUESTIONS
#   A) Which frozen candidate genes show disease-tissue support in AF-CM?
#   B) Which frozen candidate genes show disease-tissue support in HFpEF-CM?
#   C) Which frozen candidates show convergent support across BOTH diseases?
#   D) Are the frozen candidate sets collectively more HFpEF-responsive than
#      expression-matched background genes?
#   E) How similar are AF-CM and HFpEF-CM transcriptome-wide directions?
#
# HFpEF SUPPORT DEFINITIONS
#   - Transcriptome-wide FDR support: adj.P.Val < 0.05
#   - Robust FDR support: adj.P.Val < 0.05 AND low_background == TRUE
#   - Nominal support: P.Value < 0.05 but adj.P.Val >= 0.05
#   - FDR support with background flag is retained but explicitly labeled
#
# AF SUPPORT DEFINITIONS
#   Reuse the already-frozen STEP12C classification:
#     Transcriptome_wide_FDR_support
#     Nominal_support_only
#     No_statistical_support
#     Not_testable
#
# SET-LEVEL TEST
#   10,000 expression-decile-matched permutations using
#   mean[-log10(HFpEF P)] as the prespecified gene-set statistic.
#
# OUTPUT
#   D:/A/data/STEP13_CROSS_DISEASE_INTEGRATION/
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))
set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

AF_FDR <- 0.05
HFPEF_FDR <- 0.05
NOMINAL_P <- 0.05

N_PERM <- 10000L
N_EXPR_BINS <- 10L

TOP_N_FIGURE <- 30L

# ==============================================================================
# 1. INPUT PATHS
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP12C — frozen genetics × AF-CM integration
# ------------------------------------------------------------------------------

STEP12C_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA",
  "STEP12C_GENETICS_CM_INTEGRATION"
)

STEP12C_READY_FILE <- file.path(
  STEP12C_ROOT,
  "00_QC",
  "STEP12C_readiness.csv"
)

AF_TOP50_FILE <- file.path(
  STEP12C_ROOT,
  "01_TOP50_PRIMARY",
  "STEP12C_TOP50_gene_level_disease_validation_matrix.csv"
)

AF_ALL_CAND_FILE <- file.path(
  STEP12C_ROOT,
  "02_ALL_CANDIDATES",
  "STEP12C_ALL_gene_level_disease_validation_matrix.csv.gz"
)

# Raw AF-CM transcriptome for transcriptome-wide cross-disease audit.
STEP12B_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA",
  "STEP12B_AF_vs_CTRL_DE"
)

AF_CM_ALL_FILE <- file.path(
  STEP12B_ROOT,
  "01_PRIMARY_CM",
  "STEP12B_Cardiomyocytes_AF_vs_CTRL_all_genes.csv.gz"
)

STEP12B_READY_FILE <- file.path(
  STEP12B_ROOT,
  "00_QC",
  "STEP12B_readiness.csv"
)

# ------------------------------------------------------------------------------
# STEP13B — frozen HFpEF-CM disease-tissue result
# ------------------------------------------------------------------------------

STEP13B_ROOT <- file.path(
  DATA_ROOT,
  "STEP13_HFpEF_SCP3342_snRNA",
  "STEP13B_HFpEF_vs_CTRL_CM_DE"
)

STEP13B_READY_FILE <- file.path(
  STEP13B_ROOT,
  "00_QC",
  "STEP13B_readiness.csv"
)

HFPEF_CM_ALL_FILE <- file.path(
  STEP13B_ROOT,
  "02_PRIMARY_DE",
  "STEP13B_Cardiomyocyte_HFpEF_vs_CTRL_all_genes.csv.gz"
)

HFPEF_SENS_FILE <- file.path(
  STEP13B_ROOT,
  "03_SENSITIVITY",
  "STEP13B_primary_vs_sensitivity_gene_comparison.csv.gz"
)

HFPEF_SUMMARY_FILE <- file.path(
  STEP13B_ROOT,
  "00_QC",
  "STEP13B_DE_summary.csv"
)

# ------------------------------------------------------------------------------
# Optional original frozen upstream candidate tables — audit only.
# ------------------------------------------------------------------------------

STEP10B_ROOT <- file.path(
  DATA_ROOT,
  "STEP10B",
  "STEP10B_V3_FINAL_LOCUS_CENTRIC_RESULTS"
)

STEP10B_TOP50_FILE <- file.path(
  STEP10B_ROOT,
  "02_MAIN_TABLES",
  "Table2_priority_locus_gene_candidates_TOP50.csv"
)

STEP10B_ALL_FILE <- file.path(
  STEP10B_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS2_all_locus_gene_evidence.csv"
)

# ==============================================================================
# 2. OUTPUT
# ==============================================================================

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP13_CROSS_DISEASE_INTEGRATION"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
TOP50_DIR <- file.path(OUT_ROOT, "01_TOP50_CROSS_DISEASE")
ALL_DIR <- file.path(OUT_ROOT, "02_ALL_CANDIDATES_CROSS_DISEASE")
SET_DIR <- file.path(OUT_ROOT, "03_SET_LEVEL")
TX_DIR <- file.path(OUT_ROOT, "04_TRANSCRIPTOME_WIDE")
FIG_DIR <- file.path(OUT_ROOT, "05_FIGURES")
SOURCE_DIR <- file.path(OUT_ROOT, "06_SOURCE_DATA")

for (d in c(
  OUT_ROOT,
  QC_DIR,
  TOP50_DIR,
  ALL_DIR,
  SET_DIR,
  TX_DIR,
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

first_existing_col <- function(
  DT,
  candidates,
  required = TRUE,
  label = ""
) {

  hit <- candidates[
    candidates %in%
      names(DT)
  ]

  if (!length(hit)) {

    if (required) {
      stop(
        "Could not find required column for ",
        label,
        ". Tried: ",
        paste(
          candidates,
          collapse = ", "
        )
      )
    }

    return(
      NA_character_
    )
  }

  hit[1]
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

safe_log10p <- function(p) {

  p <- safe_num(p)

  p[
    is.finite(p) &
      p <= 0
  ] <- .Machine$double.xmin

  -log10(p)
}

# ------------------------------------------------------------------------------
# Standardize AF transcriptome-wide STEP12B result.
# ------------------------------------------------------------------------------

standardize_af_de <- function(DT) {

  GENE_COL <- first_existing_col(
    DT,
    c(
      "gene_symbol",
      "gene_key",
      "symbol",
      "gene"
    ),
    TRUE,
    "AF gene"
  )

  LOGFC_COL <- first_existing_col(
    DT,
    c(
      "logFC",
      "log2FC",
      "estimate"
    ),
    TRUE,
    "AF logFC"
  )

  EXPR_COL <- first_existing_col(
    DT,
    c(
      "logCPM",
      "AveExpr",
      "baseMean"
    ),
    FALSE,
    "AF expression"
  )

  P_COL <- first_existing_col(
    DT,
    c(
      "PValue",
      "P.Value",
      "pvalue",
      "p_value",
      "P"
    ),
    TRUE,
    "AF P"
  )

  FDR_COL <- first_existing_col(
    DT,
    c(
      "FDR",
      "adj.P.Val",
      "padj",
      "qvalue"
    ),
    TRUE,
    "AF FDR"
  )

  OUT <- data.table(
    gene = norm_gene(
      DT[[GENE_COL]]
    ),
    AF_logFC = safe_num(
      DT[[LOGFC_COL]]
    ),
    AF_expr =
      if (!is.na(EXPR_COL)) {
        safe_num(
          DT[[EXPR_COL]]
        )
      } else {
        NA_real_
      },
    AF_P = safe_num(
      DT[[P_COL]]
    ),
    AF_FDR = safe_num(
      DT[[FDR_COL]]
    )
  )

  OUT <- OUT[
    !is.na(gene)
  ]

  # One row per symbol: retain the smallest-P representation.
  setorder(
    OUT,
    gene,
    AF_P
  )

  OUT <- OUT[
    !duplicated(gene)
  ]

  OUT
}

# ------------------------------------------------------------------------------
# Standardize HFpEF transcriptome-wide STEP13B result.
# ------------------------------------------------------------------------------

standardize_hfpef_de <- function(DT) {

  GENE_COL <- first_existing_col(
    DT,
    c(
      "gene_symbol",
      "gene_key",
      "symbol",
      "gene"
    ),
    TRUE,
    "HFpEF gene"
  )

  LOGFC_COL <- first_existing_col(
    DT,
    c(
      "logFC",
      "log2FC",
      "estimate"
    ),
    TRUE,
    "HFpEF logFC"
  )

  EXPR_COL <- first_existing_col(
    DT,
    c(
      "AveExpr",
      "logCPM",
      "baseMean"
    ),
    TRUE,
    "HFpEF expression"
  )

  P_COL <- first_existing_col(
    DT,
    c(
      "P.Value",
      "PValue",
      "pvalue",
      "p_value",
      "P"
    ),
    TRUE,
    "HFpEF P"
  )

  FDR_COL <- first_existing_col(
    DT,
    c(
      "adj.P.Val",
      "FDR",
      "padj",
      "qvalue"
    ),
    TRUE,
    "HFpEF FDR"
  )

  BKG_COL <- first_existing_col(
    DT,
    c(
      "background_heuristic"
    ),
    FALSE,
    "HFpEF background heuristic"
  )

  LOW_BKG_COL <- first_existing_col(
    DT,
    c(
      "low_background",
      "FDR05_low_background"
    ),
    FALSE,
    "HFpEF low-background flag"
  )

  OUT <- data.table(
    gene = norm_gene(
      DT[[GENE_COL]]
    ),
    HFpEF_logFC =
      safe_num(
        DT[[LOGFC_COL]]
      ),
    HFpEF_AveExpr =
      safe_num(
        DT[[EXPR_COL]]
      ),
    HFpEF_P =
      safe_num(
        DT[[P_COL]]
      ),
    HFpEF_FDR =
      safe_num(
        DT[[FDR_COL]]
      ),
    HFpEF_background_heuristic =
      if (!is.na(BKG_COL)) {
        safe_num(
          DT[[BKG_COL]]
        )
      } else {
        NA_real_
      },
    HFpEF_low_background =
      if (!is.na(LOW_BKG_COL)) {

        z <- DT[[LOW_BKG_COL]]

        if (is.logical(z)) {
          z
        } else {
          toupper(
            as.character(z)
          ) %in%
            c(
              "TRUE",
              "T",
              "1",
              "YES"
            )
        }

      } else {
        NA
      }
  )

  OUT <- OUT[
    !is.na(gene)
  ]

  setorder(
    OUT,
    gene,
    HFpEF_P
  )

  OUT <- OUT[
    !duplicated(gene)
  ]

  OUT
}

# ------------------------------------------------------------------------------
# Add HFpEF support class.
# ------------------------------------------------------------------------------

add_hfpef_support <- function(DT) {

  X <- copy(DT)

  X[
    ,
    HFpEF_tested :=
      is.finite(
        HFpEF_P
      ) &
      is.finite(
        HFpEF_FDR
      )
  ]

  X[
    ,
    HFpEF_FDR_support :=
      HFpEF_tested &
      HFpEF_FDR <
        HFPEF_FDR
  ]

  X[
    ,
    HFpEF_FDR_robust_support :=
      HFpEF_FDR_support &
      (
        is.na(
          HFpEF_low_background
        ) |
        HFpEF_low_background
      )
  ]

  X[
    ,
    HFpEF_nominal_or_better :=
      HFpEF_tested &
      HFpEF_P <
        NOMINAL_P
  ]

  X[
    ,
    HFpEF_disease_tissue_support :=
      fcase(
        !HFpEF_tested,
        "Not_testable",

        HFpEF_FDR_support &
          !is.na(
            HFpEF_low_background
          ) &
          !HFpEF_low_background,
        "Transcriptome_wide_FDR_support_background_flagged",

        HFpEF_FDR_support,
        "Transcriptome_wide_FDR_robust_support",

        HFpEF_P <
          NOMINAL_P,
        "Nominal_support_only",

        default =
          "No_statistical_support"
      )
  ]

  X
}

# ------------------------------------------------------------------------------
# Add cross-disease support class and direction.
# ------------------------------------------------------------------------------

add_cross_disease_class <- function(DT) {

  X <- copy(DT)

  # AF support is frozen from STEP12C.
  X[
    ,
    AF_tested :=
      is.finite(
        CM_PValue
      ) &
      is.finite(
        CM_FDR
      )
  ]

  X[
    ,
    AF_FDR_support :=
      AF_tested &
      CM_FDR <
        AF_FDR
  ]

  X[
    ,
    AF_nominal_or_better :=
      AF_tested &
      CM_PValue <
        NOMINAL_P
  ]

  X <- add_hfpef_support(
    X
  )

  X[
    ,
    direction_relation :=
      fcase(
        !AF_tested |
          !HFpEF_tested |
          !is.finite(
            CM_logFC
          ) |
          !is.finite(
            HFpEF_logFC
          ),
        "Not_evaluable",

        CM_logFC ==
          0 |
          HFpEF_logFC ==
            0,
        "Zero_effect",

        sign(
          CM_logFC
        ) ==
          sign(
            HFpEF_logFC
          ),
        "Same_direction",

        default =
          "Opposite_direction"
      )
  ]

  X[
    ,
    cross_disease_support_class :=
      fcase(
        !AF_tested |
          !HFpEF_tested,
        "Not_testable_in_both",

        AF_FDR_support &
          HFpEF_FDR_robust_support,
        "Dual_FDR_robust_support",

        AF_FDR_support &
          HFpEF_FDR_support,
        "Dual_FDR_HFpEF_background_flagged",

        AF_FDR_support &
          HFpEF_nominal_or_better,
        "AF_FDR_HFpEF_nominal",

        HFpEF_FDR_support &
          AF_nominal_or_better,
        "HFpEF_FDR_AF_nominal",

        AF_nominal_or_better &
          HFpEF_nominal_or_better,
        "Dual_nominal_only",

        AF_FDR_support,
        "AF_FDR_only",

        HFpEF_FDR_support,
        "HFpEF_FDR_only",

        AF_nominal_or_better |
          HFpEF_nominal_or_better,
        "Single_nominal_only",

        default =
          "No_statistical_support"
      )
  ]

  X[
    ,
    dual_nominal_or_better :=
      AF_nominal_or_better &
      HFpEF_nominal_or_better
  ]

  X[
    ,
    dual_FDR :=
      AF_FDR_support &
      HFpEF_FDR_support
  ]

  X[
    ,
    dual_FDR_robust :=
      AF_FDR_support &
      HFpEF_FDR_robust_support
  ]

  X[
    ,
    convergent_dual_FDR_robust :=
      dual_FDR_robust &
      direction_relation ==
        "Same_direction"
  ]

  X
}

# ------------------------------------------------------------------------------
# Expression-decile-matched candidate set test.
# Statistic = mean[-log10(P)].
# ------------------------------------------------------------------------------

matched_gene_set_test <- function(
  candidate_genes,
  de,
  set_name,
  n_perm = 10000L,
  n_bins = 10L
) {

  D <- copy(
    de
  )

  D <- D[
    is.finite(
      HFpEF_P
    ) &
      is.finite(
        HFpEF_AveExpr
      ) &
      !is.na(
        gene
      )
  ]

  D[
    ,
    score :=
      safe_log10p(
        HFpEF_P
      )
  ]

  if (
    nrow(D) <
      100L
  ) {
    stop(
      "HFpEF DE universe too small for matched test."
    )
  }

  # Robust quantile bins. unique() handles tied quantiles.
  brks <- unique(
    as.numeric(
      quantile(
        D$HFpEF_AveExpr,
        probs =
          seq(
            0,
            1,
            length.out =
              n_bins +
              1L
          ),
        na.rm = TRUE,
        names = FALSE
      )
    )
  )

  if (
    length(brks) <
      3L
  ) {
    stop(
      "Could not construct expression bins."
    )
  }

  brks[1] <- -Inf
  brks[length(brks)] <- Inf

  D[
    ,
    expr_bin :=
      cut(
        HFpEF_AveExpr,
        breaks = brks,
        include.lowest = TRUE,
        labels = FALSE
      )
  ]

  CANDS <- unique(
    norm_gene(
      candidate_genes
    )
  )

  CANDS <- CANDS[
    !is.na(CANDS)
  ]

  C <- D[
    gene %in%
      CANDS
  ]

  if (
    nrow(C) <
      3L
  ) {

    return(
      data.table(
        candidate_set =
          set_name,
        n_candidate_input =
          length(CANDS),
        n_candidate_tested =
          nrow(C),
        observed_mean_neglog10P =
          if (
            nrow(C)
          ) {
            mean(
              C$score
            )
          } else {
            NA_real_
          },
        empirical_P =
          NA_real_,
        permutations =
          n_perm,
        note =
          "Fewer than 3 candidate genes testable"
      )
    )
  }

  NEED <- table(
    C$expr_bin
  )

  OBS <- mean(
    C$score
  )

  NULL_STAT <- numeric(
    n_perm
  )

  set.seed(
    9527
  )

  for (b in seq_len(
    n_perm
  )) {

    sampled <- character(
      0
    )

    for (bin_name in names(
      NEED
    )) {

      n_need <- as.integer(
        NEED[
          bin_name
        ]
      )

      pool <- D[
        expr_bin ==
          as.integer(
            bin_name
          ),
        gene
      ]

      # Prefer noncandidate background.
      pool_non_candidate <- setdiff(
        pool,
        CANDS
      )

      if (
        length(
          pool_non_candidate
        ) >=
          n_need
      ) {

        pool_use <-
          pool_non_candidate

      } else {

        pool_use <-
          pool
      }

      sampled <- c(
        sampled,
        sample(
          pool_use,
          size =
            n_need,
          replace =
            length(
              pool_use
            ) <
            n_need
        )
      )
    }

    NULL_STAT[b] <- mean(
      D[
        match(
          sampled,
          gene
        ),
        score
      ],
      na.rm = TRUE
    )
  }

  EMP <- (
    1 +
      sum(
        NULL_STAT >=
          OBS,
        na.rm = TRUE
      )
  ) /
    (
      1 +
        n_perm
    )

  data.table(
    candidate_set =
      set_name,
    n_candidate_input =
      length(CANDS),
    n_candidate_tested =
      nrow(C),
    observed_mean_neglog10P =
      OBS,
    null_mean =
      mean(
        NULL_STAT
      ),
    null_sd =
      sd(
        NULL_STAT
      ),
    empirical_P =
      EMP,
    permutations =
      n_perm,
    note =
      "HFpEF expression-decile matched"
  )
}

# ==============================================================================
# 5. HARD INPUT CHECK / READINESS
# ==============================================================================

required_files <- c(
  STEP12C_READY_FILE,
  AF_TOP50_FILE,
  AF_ALL_CAND_FILE,
  STEP12B_READY_FILE,
  AF_CM_ALL_FILE,
  STEP13B_READY_FILE,
  HFPEF_CM_ALL_FILE,
  HFPEF_SENS_FILE
)

missing <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing
  )
) {

  stop(
    paste0(
      "Missing required STEP13C input(s):\n",
      paste(
        missing,
        collapse = "\n"
      )
    )
  )
}

READY12C <- fread(
  STEP12C_READY_FILE
)

READY12B <- fread(
  STEP12B_READY_FILE
)

READY13B <- fread(
  STEP13B_READY_FILE
)

if (
  !all(
    bool_pass(
      READY12C$pass
    )
  )
) {
  stop(
    "STEP12C is not fully PASS."
  )
}

if (
  !all(
    bool_pass(
      READY12B$pass
    )
  )
) {
  stop(
    "STEP12B is not fully PASS."
  )
}

if (
  !all(
    bool_pass(
      READY13B$pass
    )
  )
) {
  stop(
    "STEP13B is not fully PASS."
  )
}

# ==============================================================================
# 6. LOAD FROZEN AF-CM CANDIDATE MATRICES
# ==============================================================================

TOP50_AF <- fread(
  AF_TOP50_FILE
)

ALL_AF <- fread(
  AF_ALL_CAND_FILE
)

required_af_candidate_cols <- c(
  "candidate_gene",
  "CM_logFC",
  "CM_PValue",
  "CM_FDR",
  "CM_disease_tissue_support"
)

for (
  nm in required_af_candidate_cols
) {

  if (
    !nm %in%
      names(
        TOP50_AF
      )
  ) {
    stop(
      "TOP50 AF integration missing column: ",
      nm
    )
  }

  if (
    !nm %in%
      names(
        ALL_AF
      )
  ) {
    stop(
      "ALL-candidate AF integration missing column: ",
      nm
    )
  }
}

TOP50_AF[
  ,
  candidate_gene :=
    norm_gene(
      candidate_gene
    )
]

ALL_AF[
  ,
  candidate_gene :=
    norm_gene(
      candidate_gene
    )
]

TOP50_AF <- TOP50_AF[
  !is.na(
    candidate_gene
  )
]

ALL_AF <- ALL_AF[
  !is.na(
    candidate_gene
  )
]

# Guarantee one row per candidate gene.
setorder(
  TOP50_AF,
  candidate_gene,
  CM_PValue
)

TOP50_AF <- TOP50_AF[
  !duplicated(
    candidate_gene
  )
]

setorder(
  ALL_AF,
  candidate_gene,
  CM_PValue
)

ALL_AF <- ALL_AF[
  !duplicated(
    candidate_gene
  )
]

# ==============================================================================
# 7. LOAD / STANDARDIZE TRANSCRIPTOME-WIDE AF AND HFpEF CM DE
# ==============================================================================

AF_RAW <- fread(
  AF_CM_ALL_FILE
)

HF_RAW <- fread(
  HFPEF_CM_ALL_FILE
)

AF_DE <- standardize_af_de(
  AF_RAW
)

HF_DE <- standardize_hfpef_de(
  HF_RAW
)

if (
  nrow(
    AF_DE
  ) <
    1000L
) {
  stop(
    "AF-CM transcriptome-wide DE universe unexpectedly small."
  )
}

if (
  nrow(
    HF_DE
  ) <
    1000L
) {
  stop(
    "HFpEF-CM transcriptome-wide DE universe unexpectedly small."
  )
}

if (
  !all(
    is.finite(
      AF_DE$AF_P
    )
  ) ||
  !all(
    is.finite(
      AF_DE$AF_FDR
    )
  )
) {
  stop(
    "AF-CM transcriptome-wide P/FDR contains non-finite values."
  )
}

if (
  !all(
    is.finite(
      HF_DE$HFpEF_P
    )
  ) ||
  !all(
    is.finite(
      HF_DE$HFpEF_FDR
    )
  )
) {
  stop(
    "HFpEF-CM transcriptome-wide P/FDR contains non-finite values."
  )
}

# ==============================================================================
# 8. MERGE FROZEN CANDIDATES WITH HFpEF-CM
# ==============================================================================

merge_hfpef <- function(AF_CAND) {

  X <- merge(
    AF_CAND,
    HF_DE,
    by.x =
      "candidate_gene",
    by.y =
      "gene",
    all.x = TRUE,
    sort = FALSE
  )

  X <- add_cross_disease_class(
    X
  )

  X
}

TOP50 <- merge_hfpef(
  TOP50_AF
)

ALLC <- merge_hfpef(
  ALL_AF
)

# ==============================================================================
# 9. EVIDENCE ORDERING — NO CANDIDATE RESELECTION
# ==============================================================================

support_order <- c(
  "Dual_FDR_robust_support",
  "Dual_FDR_HFpEF_background_flagged",
  "AF_FDR_HFpEF_nominal",
  "HFpEF_FDR_AF_nominal",
  "Dual_nominal_only",
  "AF_FDR_only",
  "HFpEF_FDR_only",
  "Single_nominal_only",
  "No_statistical_support",
  "Not_testable_in_both"
)

for (nm in c(
  "TOP50",
  "ALLC"
)) {

  X <- get(
    nm
  )

  X[
    ,
    support_rank :=
      match(
        cross_disease_support_class,
        support_order
      )
  ]

  setorder(
    X,
    support_rank,
    HFpEF_P,
    CM_PValue
  )

  X[
    ,
    support_rank := NULL
  ]

  assign(
    nm,
    X
  )
}

# ==============================================================================
# 10. HFpEF EXPRESSION-MATCHED SET-LEVEL TEST
# ==============================================================================

TOP_SET <- matched_gene_set_test(
  candidate_genes =
    TOP50$candidate_gene,
  de =
    HF_DE,
  set_name =
    "STEP10B_V3_TOP50",
  n_perm =
    N_PERM,
  n_bins =
    N_EXPR_BINS
)

ALL_SET <- matched_gene_set_test(
  candidate_genes =
    ALLC$candidate_gene,
  de =
    HF_DE,
  set_name =
    "STEP10B_V3_ALL_LOCUS_GENE",
  n_perm =
    N_PERM,
  n_bins =
    N_EXPR_BINS
)

SET_ENRICH <- rbindlist(
  list(
    TOP_SET,
    ALL_SET
  ),
  use.names = TRUE,
  fill = TRUE
)

# CODE RELEASE PATCH: downstream section replaced by validated STEP13C V3 scope-safe implementation.
# ==============================================================================
# 11. TRANSCRIPTOME-WIDE AF-CM vs HFpEF-CM AUDIT
# ==============================================================================

TX <- merge(
  AF_DE,
  HF_DE,
  by = "gene",
  all = FALSE
)

# Explicit vector assignments avoid data.table name masking between
# the AF_FDR column and the AF FDR threshold constant.
TX$AF_FDR05 <- TX$AF_FDR < 0.05
TX$HFpEF_FDR05 <- TX$HFpEF_FDR < HFPEF_FDR

# Explicit vector assignments avoid data.table evaluation/scope ambiguity.
TX$same_direction <-
  sign(TX$AF_logFC) ==
  sign(TX$HFpEF_logFC)

TX$dual_FDR05 <-
  as.logical(TX$AF_FDR05) &
  as.logical(TX$HFpEF_FDR05)

if (!"dual_FDR05" %in% names(TX)) {
  stop("Failed to create TX$dual_FDR05.")
}

GLOBAL_SPEARMAN <- suppressWarnings(
  cor(
    TX$AF_logFC,
    TX$HFpEF_logFC,
    method = "spearman",
    use = "complete.obs"
  )
)

DUAL <- TX[
  which(TX$dual_FDR05 %in% TRUE)
]

DUAL_SAME_DIRECTION_FRAC <-
  if (
    nrow(DUAL)
  ) {
    mean(
      DUAL$same_direction,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

# Fisher association of transcriptome-wide significance status.
TAB_SIG <- table(
  AF =
    TX$AF_FDR05,
  HFpEF =
    TX$HFpEF_FDR05
)

FISHER_SIG <- fisher.test(
  TAB_SIG
)

TX_SUMMARY <- data.table(
  metric = c(
    "common_CM_genes",
    "AF_FDR05_genes",
    "HFpEF_FDR05_genes",
    "dual_FDR05_genes",
    "dual_FDR05_same_direction_genes",
    "dual_FDR05_opposite_direction_genes",
    "dual_FDR05_same_direction_fraction",
    "AF_HFpEF_logFC_Spearman_all_common_genes",
    "Fisher_significance_overlap_OR",
    "Fisher_significance_overlap_P"
  ),
  value = c(
    nrow(TX),
    sum(
      TX$AF_FDR05
    ),
    sum(
      TX$HFpEF_FDR05
    ),
    nrow(DUAL),
    sum(
      DUAL$same_direction,
      na.rm = TRUE
    ),
    sum(
      !DUAL$same_direction,
      na.rm = TRUE
    ),
    DUAL_SAME_DIRECTION_FRAC,
    GLOBAL_SPEARMAN,
    unname(
      FISHER_SIG$estimate
    ),
    FISHER_SIG$p.value
  )
)

# ==============================================================================
# 12. CANDIDATE-LEVEL SUMMARY
# ==============================================================================

summarize_candidate_set <- function(
  X,
  set_name
) {

  data.table(
    candidate_set =
      set_name,
    n_candidate_genes =
      uniqueN(
        X$candidate_gene
      ),
    n_AF_tested =
      sum(
        X$AF_tested,
        na.rm = TRUE
      ),
    n_HFpEF_tested =
      sum(
        X$HFpEF_tested,
        na.rm = TRUE
      ),
    n_tested_in_both =
      sum(
        X$AF_tested &
          X$HFpEF_tested,
        na.rm = TRUE
      ),
    n_AF_FDR =
      sum(
        X$AF_FDR_support,
        na.rm = TRUE
      ),
    n_HFpEF_FDR =
      sum(
        X$HFpEF_FDR_support,
        na.rm = TRUE
      ),
    n_HFpEF_FDR_robust =
      sum(
        X$HFpEF_FDR_robust_support,
        na.rm = TRUE
      ),
    n_dual_FDR =
      sum(
        X$dual_FDR,
        na.rm = TRUE
      ),
    n_dual_FDR_robust =
      sum(
        X$dual_FDR_robust,
        na.rm = TRUE
      ),
    n_dual_nominal_or_better =
      sum(
        X$dual_nominal_or_better,
        na.rm = TRUE
      ),
    n_convergent_dual_FDR_robust =
      sum(
        X$convergent_dual_FDR_robust,
        na.rm = TRUE
      ),
    n_dual_nominal_same_direction =
      sum(
        X$dual_nominal_or_better &
          X$direction_relation ==
            "Same_direction",
        na.rm = TRUE
      ),
    n_dual_nominal_opposite_direction =
      sum(
        X$dual_nominal_or_better &
          X$direction_relation ==
            "Opposite_direction",
        na.rm = TRUE
      )
  )
}

CAND_SUMMARY <- rbindlist(
  list(
    summarize_candidate_set(
      TOP50,
      "TOP50"
    ),
    summarize_candidate_set(
      ALLC,
      "ALL_LOCUS_GENE"
    )
  )
)

# ==============================================================================
# 13. WRITE FINAL TABLES
# ==============================================================================

TOP50_OUT <- file.path(
  TOP50_DIR,
  "STEP13C_TOP50_cross_disease_evidence_matrix.csv"
)

ALL_OUT <- file.path(
  ALL_DIR,
  "STEP13C_ALL_cross_disease_evidence_matrix.csv.gz"
)

fwrite(
  TOP50,
  TOP50_OUT
)

fwrite(
  ALLC,
  ALL_OUT,
  compress = "gzip"
)

fwrite(
  TOP50[
    dual_nominal_or_better ==
      TRUE
  ],
  file.path(
    TOP50_DIR,
    "STEP13C_TOP50_dual_nominal_or_better.csv"
  )
)

fwrite(
  TOP50[
    dual_FDR ==
      TRUE
  ],
  file.path(
    TOP50_DIR,
    "STEP13C_TOP50_dual_FDR_support.csv"
  )
)

fwrite(
  TOP50[
    dual_FDR_robust ==
      TRUE
  ],
  file.path(
    TOP50_DIR,
    "STEP13C_TOP50_dual_FDR_robust_support.csv"
  )
)

fwrite(
  TOP50[
    convergent_dual_FDR_robust ==
      TRUE
  ],
  file.path(
    TOP50_DIR,
    "STEP13C_TOP50_convergent_dual_FDR_robust.csv"
  )
)

fwrite(
  ALLC[
    dual_FDR_robust ==
      TRUE
  ],
  file.path(
    ALL_DIR,
    "STEP13C_ALL_dual_FDR_robust_support.csv"
  )
)

fwrite(
  SET_ENRICH,
  file.path(
    SET_DIR,
    "STEP13C_HFpEF_expression_matched_candidate_set_enrichment.csv"
  )
)

fwrite(
  TX_SUMMARY,
  file.path(
    TX_DIR,
    "STEP13C_AF_HFpEF_CM_transcriptome_overlap_summary.csv"
  )
)

fwrite(
  TX,
  file.path(
    TX_DIR,
    "STEP13C_AF_HFpEF_CM_common_transcriptome.csv.gz"
  ),
  compress = "gzip"
)

fwrite(
  CAND_SUMMARY,
  file.path(
    QC_DIR,
    "STEP13C_candidate_cross_disease_summary.csv"
  )
)

# ==============================================================================
# 14. OPTIONAL AUDIT — ORIGINAL STEP10B FILES PRESENT
# ==============================================================================

UPSTREAM_AUDIT <- data.table(
  file = c(
    STEP10B_TOP50_FILE,
    STEP10B_ALL_FILE
  ),
  exists = file.exists(
    c(
      STEP10B_TOP50_FILE,
      STEP10B_ALL_FILE
    )
  )
)

fwrite(
  UPSTREAM_AUDIT,
  file.path(
    QC_DIR,
    "STEP13C_STEP10B_original_candidate_file_audit.csv"
  )
)

# ==============================================================================
# 15. FIGURES — SINGLE FIGURE FILES, NO LEGEND INSIDE
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure A: TOP50 AF-CM vs HFpEF-CM effect scatter.
# ------------------------------------------------------------------------------

PLOT_TOP <- TOP50[
  is.finite(
    CM_logFC
  ) &
    is.finite(
      HFpEF_logFC
    )
]

PLOT_TOP[
  ,
  plot_class :=
    fifelse(
      dual_FDR_robust,
      "Dual FDR robust",
      fifelse(
        dual_nominal_or_better,
        "Dual nominal+",
        "Other"
      )
    )
]

# Label strongest cross-disease evidence genes.
LABEL_TOP <- PLOT_TOP[
  dual_nominal_or_better ==
    TRUE
][
  order(
    pmax(
      CM_PValue,
      HFpEF_P,
      na.rm = TRUE
    )
  )
]

LABEL_TOP <- LABEL_TOP[
  seq_len(
    min(
      nrow(LABEL_TOP),
      15L
    )
  )
]

p1 <- ggplot(
  PLOT_TOP,
  aes(
    x = CM_logFC,
    y = HFpEF_logFC
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    linetype = 2
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    linetype = 2
  ) +
  geom_point(
    aes(
      shape = plot_class,
      size = pmin(
        8,
        safe_log10p(
          pmax(
            CM_PValue,
            HFpEF_P,
            na.rm = TRUE
          )
        )
      )
    ),
    alpha = 0.80
  ) +
  ggrepel::geom_text_repel(
    data = LABEL_TOP,
    aes(
      label =
        candidate_gene
    ),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  labs(
    x = "AF cardiomyocyte log2FC",
    y = "HFpEF cardiomyocyte log2FC"
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
    "Figure_STEP13C_TOP50_AF_vs_HFpEF_CM_effects.tiff"
  ),
  p1,
  width = 5.6,
  height = 5.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13C_TOP50_AF_vs_HFpEF_CM_effects.pdf"
  ),
  p1,
  width = 5.6,
  height = 5.0,
  units = "in"
)

fwrite(
  PLOT_TOP,
  file.path(
    SOURCE_DIR,
    "Figure_STEP13C_TOP50_AF_vs_HFpEF_CM_effects_source.csv"
  )
)

# ------------------------------------------------------------------------------
# Figure B: evidence matrix for the strongest TOP50 candidates.
# Two rows per candidate: AF and HFpEF.
# ------------------------------------------------------------------------------

PLOT_EVID <- copy(
  TOP50
)

PLOT_EVID[
  ,
  score_for_plot :=
    pmin(
      safe_log10p(
        CM_PValue
      ),
      safe_log10p(
        HFpEF_P
      ),
      na.rm = TRUE
    )
]

PLOT_EVID[
  !is.finite(
    score_for_plot
  ),
  score_for_plot :=
    -Inf
]

setorder(
  PLOT_EVID,
  -dual_FDR_robust,
  -dual_FDR,
  -dual_nominal_or_better,
  -score_for_plot
)

PLOT_EVID <- PLOT_EVID[
  seq_len(
    min(
      nrow(PLOT_EVID),
      TOP_N_FIGURE
    )
  )
]

LONG <- rbindlist(
  list(
    PLOT_EVID[
      ,
      .(
        candidate_gene,
        disease = "AF",
        logFC =
          CM_logFC,
        P =
          CM_PValue,
        FDR =
          CM_FDR
      )
    ],
    PLOT_EVID[
      ,
      .(
        candidate_gene,
        disease = "HFpEF",
        logFC =
          HFpEF_logFC,
        P =
          HFpEF_P,
        FDR =
          HFpEF_FDR
      )
    ]
  )
)

LONG[
  ,
  neglog10P :=
    safe_log10p(
      P
    )
]

LONG[
  ,
  candidate_gene :=
    factor(
      candidate_gene,
      levels =
        rev(
          unique(
            PLOT_EVID$candidate_gene
          )
        )
    )
]

LONG[
  ,
  disease :=
    factor(
      disease,
      levels = c(
        "AF",
        "HFpEF"
      )
    )
]

p2 <- ggplot(
  LONG,
  aes(
    x = disease,
    y = candidate_gene
  )
) +
  geom_point(
    aes(
      size =
        pmin(
          neglog10P,
          10
        ),
      fill =
        logFC
    ),
    shape = 21,
    colour = "black",
    stroke = 0.25
  ) +
  scale_fill_gradient2(
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13C_TOP50_cross_disease_evidence_matrix.tiff"
  ),
  p2,
  width = 4.4,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13C_TOP50_cross_disease_evidence_matrix.pdf"
  ),
  p2,
  width = 4.4,
  height = 7.2,
  units = "in"
)

fwrite(
  LONG,
  file.path(
    SOURCE_DIR,
    "Figure_STEP13C_TOP50_cross_disease_evidence_matrix_source.csv"
  )
)

# ------------------------------------------------------------------------------
# Figure C: TOP50 support-category counts.
# ------------------------------------------------------------------------------

CLASS_COUNTS <- TOP50[
  ,
  .N,
  by =
    cross_disease_support_class
][
  order(
    -N
  )
]

CLASS_COUNTS[
  ,
  support_plot :=
    factor(
      cross_disease_support_class,
      levels =
        rev(
          cross_disease_support_class
        )
    )
]

p3 <- ggplot(
  CLASS_COUNTS,
  aes(
    x = N,
    y = support_plot
  )
) +
  geom_col(
    width = 0.72
  ) +
  labs(
    x = "Number of frozen TOP50 candidate genes",
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13C_TOP50_support_class_counts.tiff"
  ),
  p3,
  width = 6.4,
  height = 4.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

fwrite(
  CLASS_COUNTS,
  file.path(
    SOURCE_DIR,
    "Figure_STEP13C_TOP50_support_class_counts_source.csv"
  )
)

# ==============================================================================
# 16. FINAL QC / READINESS
# ==============================================================================

SENS <- fread(
  HFPEF_SENS_FILE
)

SENS_LOGFC_COR <- NA_real_
SENS_DIR_AGREE <- NA_real_

if (
  all(
    c(
      "primary_logFC",
      "sensitivity_logFC"
    ) %in%
      names(
        SENS
      )
  )
) {

  SENS_LOGFC_COR <- suppressWarnings(
    cor(
      SENS$primary_logFC,
      SENS$sensitivity_logFC,
      method = "spearman",
      use = "complete.obs"
    )
  )

  SENS_DIR_AGREE <- mean(
    sign(
      SENS$primary_logFC
    ) ==
      sign(
        SENS$sensitivity_logFC
      ),
    na.rm = TRUE
  )
}

HFPEF_SENS_QC <- data.table(
  metric = c(
    "HFpEF_primary_vs_sensitivity_logFC_Spearman",
    "HFpEF_primary_vs_sensitivity_direction_agreement"
  ),
  value = c(
    SENS_LOGFC_COR,
    SENS_DIR_AGREE
  )
)

fwrite(
  HFPEF_SENS_QC,
  file.path(
    QC_DIR,
    "STEP13C_HFpEF_sensitivity_carryforward_QC.csv"
  )
)

READINESS <- data.table(
  check = c(
    "STEP12B_fully_passed",
    "STEP12C_fully_passed",
    "STEP13B_fully_passed",
    "AF_TOP50_frozen_matrix_present",
    "AF_ALL_frozen_matrix_present",
    "AF_CM_transcriptome_gt1000_genes",
    "HFpEF_CM_transcriptome_gt1000_genes",
    "TOP50_unique_candidate_genes_gt0",
    "ALL_unique_candidate_genes_gt0",
    "TOP50_HFpEF_merge_completed",
    "ALL_HFpEF_merge_completed",
    "TOP50_expression_matched_test_completed",
    "ALL_expression_matched_test_completed",
    "transcriptome_wide_common_gene_audit_completed",
    "candidate_selection_not_redefined",
    "final_TOP50_matrix_written",
    "final_ALL_matrix_written",
    "STEP13C_complete"
  ),
  pass = c(
    all(
      bool_pass(
        READY12B$pass
      )
    ),
    all(
      bool_pass(
        READY12C$pass
      )
    ),
    all(
      bool_pass(
        READY13B$pass
      )
    ),
    file.exists(
      AF_TOP50_FILE
    ),
    file.exists(
      AF_ALL_CAND_FILE
    ),
    nrow(
      AF_DE
    ) >
      1000L,
    nrow(
      HF_DE
    ) >
      1000L,
    uniqueN(
      TOP50$candidate_gene
    ) >
      0L,
    uniqueN(
      ALLC$candidate_gene
    ) >
      0L,
    any(
      TOP50$HFpEF_tested,
      na.rm = TRUE
    ),
    any(
      ALLC$HFpEF_tested,
      na.rm = TRUE
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
    nrow(
      TX
    ) >
      1000L,
    TRUE,
    file.exists(
      TOP50_OUT
    ),
    file.exists(
      ALL_OUT
    ),
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP13C_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "integration_role",
    "upstream_candidate_status",
    "AF_disease_tissue_layer",
    "HFpEF_disease_tissue_layer",
    "primary_cell_type_both_diseases",
    "AF_support_definition",
    "HFpEF_support_definition",
    "HFpEF_robust_support_definition",
    "nominal_support_definition",
    "candidate_reselection_after_disease_tissue",
    "HFpEF_set_level_test",
    "HFpEF_set_level_permutations",
    "expression_matching",
    "transcriptome_wide_cross_disease_audit",
    "next_step"
  ),
  value = c(
    "Final cross-disease evidence integration",
    "Frozen before AF/HFpEF disease-tissue inspection",
    "GSE255612 left atrial patient-level cardiomyocyte pseudobulk",
    "SCP3342 ventricular-septal patient-level cardiomyocyte pseudobulk",
    "Cardiomyocyte",
    "AF CM FDR < 0.05 from frozen STEP12C",
    "HFpEF CM adj.P.Val < 0.05",
    "HFpEF adj.P.Val < 0.05 and low_background == TRUE when available",
    "P < 0.05 but FDR >= 0.05; descriptive support only",
    "No",
    "Competitive expression-matched random gene-set test using mean -log10(HFpEF P)",
    as.character(
      N_PERM
    ),
    paste0(
      N_EXPR_BINS,
      " HFpEF-CM expression bins"
    ),
    "Descriptive cross-disease concordance only; does not redefine candidates",
    "STEP14 final evidence freeze: Figure 1–6, main/supplementary tables, Methods/Results/Discussion"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP13C_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP13C_sessionInfo.txt"
  )
)

# ==============================================================================
# 17. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP13C COMPLETE — FINAL CROSS-DISEASE INTEGRATION\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Candidate-level cross-disease summary:\n"
)

print(
  CAND_SUMMARY
)

cat(
  "\nHFpEF expression-matched candidate-set enrichment:\n"
)

print(
  SET_ENRICH
)

cat(
  "\nAF-CM vs HFpEF-CM transcriptome-wide audit:\n"
)

print(
  TX_SUMMARY
)

cat(
  "\nTOP50 dual-FDR robust genes:\n"
)

print(
  TOP50[
    dual_FDR_robust ==
      TRUE,
    .(
      candidate_gene,
      CM_logFC,
      CM_PValue,
      CM_FDR,
      HFpEF_logFC,
      HFpEF_P,
      HFpEF_FDR,
      HFpEF_background_heuristic,
      direction_relation,
      cross_disease_support_class
    )
  ]
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

if (
  all(
    READINESS$pass
  )
) {

  cat(
    "\nSTEP13C PASSED.\n",
    "Core statistical analysis is complete.\n",
    "Next and final board: STEP14 final evidence/figure/table/manuscript freeze.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP13C has failed technical checks.\n",
    "Do not freeze the manuscript evidence matrix yet.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  TOP50_DIR,
  "\n",
  SET_DIR,
  "\n",
  TX_DIR,
  "\n",
  sep = ""
)
