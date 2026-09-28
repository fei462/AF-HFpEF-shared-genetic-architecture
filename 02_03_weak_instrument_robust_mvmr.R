# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP3D_QHET_WEAK_INSTRUMENT_ROBUST_MVMR_V2_FIXED.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP 3D — Weak-instrument robust MVMR sensitivity analysis
# Project: AF – HFpEF – BMI – OSA
#
# PURPOSE
#   STEP3C showed that the best pre-specified instrument set was:
#       BMI GWS (P<5e-8) + OSA P<5e-6
#   but OSA conditional F remained <10.
#
#   STEP3D therefore performs Q-statistic minimisation (qhet_mvmr) as a
#   WEAK-INSTRUMENT / PLEIOTROPY-ROBUST SENSITIVITY ANALYSIS.
#
# IMPORTANT SCIENTIFIC RULES
#   1. This script DOES NOT further relax SNP P-value thresholds.
#   2. This script DOES NOT choose a phenotypic correlation (rho) by looking
#      at outcome P-values.
#   3. qhet_mvmr requires a phenotypic correlation matrix between BMI and OSA.
#      Genetic correlation (LDSC rg) must NOT be substituted for phenotypic rho.
#   4. If no externally justified rho is available, this script runs a broad
#      rho sensitivity grid and reports robustness across rho values.
#      Such grid results are sensitivity analyses, NOT a single definitive
#      causal estimate.
#   5. If an externally justified rho later becomes available, set
#      RHO_PRIMARY to that value. The script will then run bootstrap 95% CIs
#      at that prespecified rho.
#
# INPUT
#   D:/A/data/STEP3_MR_STEP3C
#     - 02_harmonized/
#     - 03_mvmr_results/
#     - 05_reports/
#
# OUTPUT
#   D:/A/data/STEP3_MR_STEP3D
#
# MAIN OUTPUTS
#   01_results/STEP3D_qhet_rho_grid_results.csv
#   01_results/STEP3D_qhet_stability_summary.csv
#   01_results/STEP3D_qhet_primary_rho_results.csv          (only if set)
#   02_figures/*.tiff / *.pdf
#   03_reports/STEP3D_FINAL_DECISION.txt
#   03_reports/STEP3D_master_console.log
#
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

DATA_DIR <- "D:/A/data"

STEP3C_DIR <- file.path(
  DATA_DIR,
  "STEP3_MR_STEP3C"
)

OUT_DIR <- file.path(
  DATA_DIR,
  "STEP3_MR_STEP3D"
)

# ------------------------------------------------------------------
# PHENOTYPIC CORRELATION SETTINGS
# ------------------------------------------------------------------

# IMPORTANT:
# If you later obtain a defensible BMI–OSA PHENOTYPIC correlation from:
#   - individual-level data, OR
#   - a suitable external study using comparable phenotypes,
# put it here.
#
# Example only:
# RHO_PRIMARY <- 0.30
#
# Default = NA, meaning:
#   no single rho is treated as primary;
#   only sensitivity-grid analyses are performed.
RHO_PRIMARY <- NA_real_

# Optional source description if you specify RHO_PRIMARY later.
RHO_PRIMARY_SOURCE <- ""

# Broad diagnostic sensitivity grid.
# This grid is NOT a statement that every rho here is biologically plausible.
# It is used only to determine whether qhet MVMR conclusions are sensitive
# to uncertainty in the required phenotypic correlation.
RHO_GRID <- seq(
  from = -0.50,
  to   =  0.80,
  by   =  0.05
)

# qhet bootstrap only runs if RHO_PRIMARY is supplied.
QHET_BOOTSTRAP_ITERATIONS <- 1000L

# Windows qhet bootstrap is effectively single-core.
QHET_NCORES <- 1L

# Conventional conditional-F reference threshold.
F_REFERENCE <- 10

# Figures
FIG_DPI <- 600
BASE_FAMILY <- "Arial"

# Include strict BMI-GWS + OSA-GWS set as a secondary sensitivity comparator.
RUN_STRICT_GWS_COMPARATOR <- TRUE


# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

DIR_RES <- file.path(
  OUT_DIR,
  "01_results"
)

DIR_FIG <- file.path(
  OUT_DIR,
  "02_figures"
)

DIR_REP <- file.path(
  OUT_DIR,
  "03_reports"
)

for (d in c(
  OUT_DIR,
  DIR_RES,
  DIR_FIG,
  DIR_REP
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 2. PACKAGES
# ==============================================================================

options(
  repos = c(
    MRCIEU = "https://mrcieu.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)

install_if_missing <- function(pkg) {

  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    install.packages(pkg)
  }
}

for (p in c(
  "data.table",
  "ggplot2",
  "MVMR"
)) {

  try(
    install_if_missing(p),
    silent = TRUE
  )
}

if (!requireNamespace(
  "MVMR",
  quietly = TRUE
)) {

  stop(
    "MVMR package is required."
  )
}

# qhet_mvmr was corrected in MVMR 0.4.7;
# strhet_mvmr was corrected in MVMR 0.4.8.
if (
  utils::packageVersion("MVMR") <
    "0.4.8"
) {

  message(
    "Updating MVMR to >=0.4.8 ..."
  )

  try(
    install.packages(
      "MVMR",
      repos = c(
        "https://mrcieu.r-universe.dev",
        "https://cloud.r-project.org"
      )
    ),
    silent = TRUE
  )
}

if (
  utils::packageVersion("MVMR") <
    "0.4.8"
) {

  stop(
    "STEP3D requires MVMR >=0.4.8.\n",
    "Current version: ",
    as.character(
      utils::packageVersion("MVMR")
    ),
    "\nPlease update MVMR, restart R, and rerun."
  )
}

library(data.table)
library(ggplot2)


# ==============================================================================
# 3. LOGGING
# ==============================================================================

LOG_FILE <- file.path(
  DIR_REP,
  "STEP3D_master_console.log"
)

if (file.exists(LOG_FILE)) {
  unlink(LOG_FILE)
}

log_msg <- function(...) {

  x <- paste0(...)

  cat(
    x,
    "\n"
  )

  cat(
    paste0(
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      " | ",
      x,
      "\n"
    ),
    file = LOG_FILE,
    append = TRUE
  )
}

log_msg(
  "============================================================"
)

log_msg(
  "STEP 3D START"
)

log_msg(
  "Weak-instrument robust MVMR using qhet_mvmr"
)

log_msg(
  "MVMR package version = ",
  as.character(
    utils::packageVersion("MVMR")
  )
)

log_msg(
  "============================================================"
)


# ==============================================================================
# 4. INPUT CHECKS
# ==============================================================================

if (!dir.exists(STEP3C_DIR)) {

  stop(
    "STEP3C directory not found:\n",
    STEP3C_DIR
  )
}

SELECTED_FILE <- file.path(
  STEP3C_DIR,
  "05_reports",
  "STEP3C_SELECTED_strength_optimized_tiers.csv"
)

EFFECT_FILE <- file.path(
  STEP3C_DIR,
  "03_mvmr_results",
  "STEP3C_SELECTED_MVMR_effects.csv"
)

HARM_DIR <- file.path(
  STEP3C_DIR,
  "02_harmonized"
)

required_paths <- c(
  SELECTED_FILE,
  EFFECT_FILE,
  HARM_DIR
)

if (
  any(!file.exists(
    required_paths[1:2]
  )) ||
  !dir.exists(HARM_DIR)
) {

  stop(
    "Required STEP3C outputs are missing.\n",
    paste(
      required_paths,
      collapse = "\n"
    )
  )
}


# ==============================================================================
# 5. READ STEP3C SELECTION
# ==============================================================================

selected <- fread(
  SELECTED_FILE
)

step3c_ivw <- fread(
  EFFECT_FILE
)

required_sel_cols <- c(
  "combo_id",
  "outcome",
  "F_BMI",
  "F_OSA",
  "selection_status"
)

missing_sel <- setdiff(
  required_sel_cols,
  names(selected)
)

if (
  length(missing_sel) > 0
) {

  stop(
    "Selected-tier file is missing columns: ",
    paste(
      missing_sel,
      collapse = ", "
    )
  )
}

log_msg(
  "STEP3C selected tiers:"
)

for (
  i in seq_len(
    nrow(selected)
  )
) {

  log_msg(
    "  ",
    selected$outcome[i],
    ": ",
    selected$combo_id[i],
    " | F_BMI=",
    signif(
      selected$F_BMI[i],
      5
    ),
    " | F_OSA=",
    signif(
      selected$F_OSA[i],
      5
    ),
    " | ",
    selected$selection_status[i]
  )
}


# ==============================================================================
# 6. DEFINE ANALYSIS SETS
# ==============================================================================

# STEP3C strength-selected set for each outcome
analysis_sets <- selected[
  ,
  .(
    outcome,
    combo_id,
    analysis_set =
      "STEP3C_SELECTED",
    F_BMI_step3c = F_BMI,
    F_OSA_step3c = F_OSA
  )
]

# Optional strict comparator
if (
  RUN_STRICT_GWS_COMPARATOR
) {

  strict_rows <- data.table(
    outcome = c(
      "HFpEF",
      "AF"
    ),
    combo_id =
      "BMI_GWS__OSA_GWS",
    analysis_set =
      "STRICT_GWS_COMPARATOR",
    F_BMI_step3c = NA_real_,
    F_OSA_step3c = NA_real_
  )

  analysis_sets <- rbindlist(
    list(
      analysis_sets,
      strict_rows
    ),
    fill = TRUE
  )
}

# Remove duplicate set if selected == strict
analysis_sets <- unique(
  analysis_sets,
  by = c(
    "outcome",
    "combo_id",
    "analysis_set"
  )
)


# ==============================================================================
# 7. LOCATE HARMONIZED FILES
# ==============================================================================

harmonized_path <- function(
  combo_id,
  outcome
) {

  file.path(
    HARM_DIR,
    paste0(
      combo_id,
      "__to__",
      outcome,
      "__harmonized.csv"
    )
  )
}

analysis_sets[
  ,
  harmonized_file :=
    mapply(
      harmonized_path,
      combo_id,
      outcome
    )
]

analysis_sets[
  ,
  file_exists :=
    file.exists(
      harmonized_file
    )
]

fwrite(
  analysis_sets,
  file.path(
    DIR_REP,
    "STEP3D_analysis_sets.csv"
  )
)

if (
  any(!analysis_sets$file_exists)
) {

  missing_files <-
    analysis_sets[
      file_exists == FALSE,
      harmonized_file
    ]

  stop(
    "Missing STEP3C harmonized files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ==============================================================================
# 8. BUILD MVMR INPUT
# ==============================================================================

prepare_rinput <- function(
  harmonized_file
) {

  d <- fread(
    harmonized_file
  )

  required <- c(
    "SNP",
    "beta_BMI",
    "se_BMI",
    "beta_OSA",
    "se_OSA",
    "beta_OUT",
    "se_OUT"
  )

  missing <- setdiff(
    required,
    names(d)
  )

  if (
    length(missing) > 0
  ) {

    stop(
      "Harmonized file is missing columns: ",
      paste(
        missing,
        collapse = ", "
      ),
      "\nFile: ",
      harmonized_file
    )
  }

  d <- d[
    complete.cases(
      d[
        ,
        ..required
      ]
    )
  ]

  if (
    nrow(d) < 3
  ) {

    stop(
      "Too few complete variants in:\n",
      harmonized_file
    )
  }

  BX <- as.matrix(
    d[
      ,
      .(
        beta_BMI,
        beta_OSA
      )
    ]
  )

  colnames(BX) <- c(
    "BMI",
    "OSA"
  )

  seBX <- as.matrix(
    d[
      ,
      .(
        se_BMI,
        se_OSA
      )
    ]
  )

  colnames(seBX) <- c(
    "BMI",
    "OSA"
  )

  rin <- MVMR::format_mvmr(
    BXGs = BX,
    BYG = d$beta_OUT,
    seBXGs = seBX,
    seBYG = d$se_OUT,
    RSID = d$SNP
  )

  list(
    data = d,
    r_input = rin
  )
}


# ==============================================================================
# 9. PHENOTYPIC CORRELATION MATRIX
# ==============================================================================

make_pcor <- function(rho) {

  if (
    !is.finite(rho) ||
    abs(rho) >= 1
  ) {

    stop(
      "rho must be finite and satisfy |rho| < 1."
    )
  }

  pcor <- matrix(
    c(
      1,
      rho,
      rho,
      1
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(pcor) <- c(
    "BMI",
    "OSA"
  )

  colnames(pcor) <- c(
    "BMI",
    "OSA"
  )

  eig <- eigen(
    pcor,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (
    any(eig <= 0)
  ) {

    stop(
      "pcor is not positive definite at rho=",
      rho
    )
  }

  pcor
}


# ==============================================================================
# 10. ROBUST qhet RESULT PARSER
# ==============================================================================

parse_qhet <- function(
  qh,
  rho,
  outcome,
  combo_id,
  analysis_set,
  n_snp,
  ci_requested
) {

  if (
    inherits(
      qh,
      "try-error"
    )
  ) {

    return(
      data.table(
        rho = rho,
        outcome = outcome,
        combo_id = combo_id,
        analysis_set = analysis_set,
        exposure = c(
          "BMI",
          "OSA"
        ),
        n_snp = n_snp,
        beta = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        OR = NA_real_,
        OR_low = NA_real_,
        OR_high = NA_real_,
        CI_requested = ci_requested,
        status = "QHET_ERROR",
        error_message =
          as.character(qh)
      )
    )
  }

  z <- as.data.frame(
    qh,
    check.names = FALSE
  )

  # qhet_mvmr normally returns one row per exposure.
  # We avoid depending on exact column names because they changed
  # across package versions.

  numeric_cols <- names(z)[
    vapply(
      z,
      is.numeric,
      logical(1)
    )
  ]

  if (
    length(numeric_cols) < 1
  ) {

    return(
      data.table(
        rho = rho,
        outcome = outcome,
        combo_id = combo_id,
        analysis_set = analysis_set,
        exposure = c(
          "BMI",
          "OSA"
        ),
        n_snp = n_snp,
        beta = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        OR = NA_real_,
        OR_low = NA_real_,
        OR_high = NA_real_,
        CI_requested = ci_requested,
        status = "QHET_PARSE_ERROR",
        error_message =
          "No numeric columns returned."
      )
    )
  }

  # Prefer an estimate/effect column if named.
  est_candidates <- numeric_cols[
    grepl(
      "effect|estimate|beta",
      numeric_cols,
      ignore.case = TRUE
    )
  ]

  if (
    length(est_candidates) > 0
  ) {

    est_col <- est_candidates[1]

  } else {

    est_col <- numeric_cols[1]
  }

  beta <- as.numeric(
    z[[est_col]]
  )

  # Expect 2 exposures; handle unexpected structures defensively.
  if (
    length(beta) < 2
  ) {

    beta <- c(
      beta,
      rep(
        NA_real_,
        2 - length(beta)
      )
    )
  }

  beta <- beta[1:2]

  ci_low <- rep(
    NA_real_,
    2
  )

  ci_high <- rep(
    NA_real_,
    2
  )

  if (
    ci_requested &&
    length(numeric_cols) >= 3
  ) {

    low_candidates <- numeric_cols[
      grepl(
        "lower|2.5|lci",
        numeric_cols,
        ignore.case = TRUE
      )
    ]

    high_candidates <- numeric_cols[
      grepl(
        "upper|97.5|uci",
        numeric_cols,
        ignore.case = TRUE
      )
    ]

    if (
      length(low_candidates) > 0 &&
      length(high_candidates) > 0
    ) {

      ci_low <- as.numeric(
        z[[low_candidates[1]]]
      )[1:2]

      ci_high <- as.numeric(
        z[[high_candidates[1]]]
      )[1:2]

    } else {

      # Fallback: after estimate, assume next two numeric columns
      # correspond to CI limits if present.
      other <- setdiff(
        numeric_cols,
        est_col
      )

      if (
        length(other) >= 2
      ) {

        a <- as.numeric(
          z[[other[1]]]
        )[1:2]

        b <- as.numeric(
          z[[other[2]]]
        )[1:2]

        ci_low <- pmin(
          a,
          b
        )

        ci_high <- pmax(
          a,
          b
        )
      }
    }
  }

  data.table(
    rho = rho,
    outcome = outcome,
    combo_id = combo_id,
    analysis_set = analysis_set,
    exposure = c(
      "BMI",
      "OSA"
    ),
    n_snp = n_snp,
    beta = beta,
    ci_low = ci_low,
    ci_high = ci_high,
    OR = exp(beta),
    OR_low = exp(ci_low),
    OR_high = exp(ci_high),
    CI_requested = ci_requested,
    status = "OK",
    error_message = NA_character_
  )
}


# ==============================================================================
# 11. RUN qhet_mvmr ACROSS THE RHO GRID
# ==============================================================================

RHO_GRID <- sort(
  unique(
    RHO_GRID[
      is.finite(
        RHO_GRID
      ) &
      abs(
        RHO_GRID
      ) < 1
    ]
  )
)

if (
  length(RHO_GRID) == 0
) {

  stop(
    "RHO_GRID contains no valid values."
  )
}

grid_results <- list()

for (
  i in seq_len(
    nrow(
      analysis_sets
    )
  )
) {

  aa <- analysis_sets[i]

  pp <- prepare_rinput(
    aa$harmonized_file
  )

  rin <- pp$r_input

  n_snp <- nrow(
    pp$data
  )

  log_msg(
    "qhet rho-grid | ",
    aa$analysis_set,
    " | ",
    aa$combo_id,
    " | ",
    aa$outcome,
    " | nSNP=",
    n_snp
  )

  for (
    rho in RHO_GRID
  ) {

    pcor <- make_pcor(
      rho
    )

    qh <- try(
      MVMR::qhet_mvmr(
        r_input = rin,
        pcor = pcor,
        CI = FALSE,
        iterations = 100,
        ncores = QHET_NCORES
      ),
      silent = TRUE
    )

    rr <- parse_qhet(
      qh = qh,
      rho = rho,
      outcome = aa$outcome,
      combo_id = aa$combo_id,
      analysis_set =
        aa$analysis_set,
      n_snp = n_snp,
      ci_requested = FALSE
    )

    grid_results[[length(grid_results) + 1]] <- rr
  }
}

grid_results <- rbindlist(
  grid_results,
  fill = TRUE
)

fwrite(
  grid_results,
  file.path(
    DIR_RES,
    "STEP3D_qhet_rho_grid_results.csv"
  )
)


# ==============================================================================
# 12. STABILITY SUMMARY ACROSS RHO
# ==============================================================================

stability <- grid_results[
  status == "OK" &
  is.finite(beta),
  .(
    n_rho_success =
      uniqueN(rho),

    rho_min =
      min(rho),

    rho_max =
      max(rho),

    beta_min =
      min(beta),

    beta_median =
      median(beta),

    beta_max =
      max(beta),

    OR_min =
      min(OR),

    OR_median =
      median(OR),

    OR_max =
      max(OR),

    all_positive =
      all(beta > 0),

    all_negative =
      all(beta < 0),

    crosses_null =
      min(beta) <= 0 &
      max(beta) >= 0,

    max_abs_beta_deviation =
      max(
        abs(
          beta -
            median(beta)
        )
      )
  ),
  by = .(
    analysis_set,
    combo_id,
    outcome,
    exposure,
    n_snp
  )
]

stability[
  ,
  direction_stable :=
    all_positive |
    all_negative
]

stability[
  ,
  stability_class :=
    fifelse(
      direction_stable,
      "DIRECTION_STABLE_ACROSS_RHO_GRID",
      "DIRECTION_SENSITIVE_TO_RHO"
    )
]

fwrite(
  stability,
  file.path(
    DIR_RES,
    "STEP3D_qhet_stability_summary.csv"
  )
)


# ==============================================================================
# 13. COMPARE qhet GRID WITH STEP3C CONVENTIONAL IVW
# ==============================================================================

ivw_selected <- copy(
  step3c_ivw
)

needed_ivw <- c(
  "exposure",
  "outcome",
  "combo_id",
  "b",
  "OR"
)

if (
  all(
    needed_ivw %in%
      names(
        ivw_selected
      )
  )
) {

  comparison <- merge(
    stability[
      analysis_set ==
        "STEP3C_SELECTED"
    ],
    ivw_selected[
      ,
      .(
        exposure,
        outcome,
        combo_id,
        IVW_beta = b,
        IVW_OR = OR,
        IVW_p = pval,
        F_BMI = F_BMI,
        F_OSA = F_OSA,
        IVW_analysis_class =
          analysis_class
      )
    ],
    by = c(
      "exposure",
      "outcome",
      "combo_id"
    ),
    all.x = TRUE
  )

  comparison[
    ,
    IVW_direction_agrees_with_qhet_grid :=
      (
        IVW_beta > 0 &
          all_positive
      ) |
      (
        IVW_beta < 0 &
          all_negative
      )
  ]

  fwrite(
    comparison,
    file.path(
      DIR_RES,
      "STEP3D_IVW_vs_qhet_comparison.csv"
    )
  )
}


# ==============================================================================
# 14. OPTIONAL PRIMARY RHO WITH BOOTSTRAP 95% CI
# ==============================================================================

primary_results <- data.table()

if (
  is.finite(
    RHO_PRIMARY
  )
) {

  if (
    abs(
      RHO_PRIMARY
    ) >= 1
  ) {

    stop(
      "RHO_PRIMARY must satisfy |rho| < 1."
    )
  }

  log_msg(
    "Primary prespecified rho supplied: ",
    RHO_PRIMARY
  )

  if (
    !nzchar(
      trimws(
        RHO_PRIMARY_SOURCE
      )
    )
  ) {

    warning(
      "RHO_PRIMARY was supplied but RHO_PRIMARY_SOURCE is blank. ",
      "Document the source before manuscript use."
    )
  }

  primary_list <- list()

  selected_sets <- analysis_sets[
    analysis_set ==
      "STEP3C_SELECTED"
  ]

  for (
    i in seq_len(
      nrow(
        selected_sets
      )
    )
  ) {

    aa <- selected_sets[i]

    pp <- prepare_rinput(
      aa$harmonized_file
    )

    pcor <- make_pcor(
      RHO_PRIMARY
    )

    log_msg(
      "Bootstrap qhet | outcome=",
      aa$outcome,
      " | rho=",
      RHO_PRIMARY,
      " | iterations=",
      QHET_BOOTSTRAP_ITERATIONS
    )

    qh <- try(
      MVMR::qhet_mvmr(
        r_input =
          pp$r_input,
        pcor = pcor,
        CI = TRUE,
        iterations =
          QHET_BOOTSTRAP_ITERATIONS,
        ncores =
          QHET_NCORES
      ),
      silent = TRUE
    )

    rr <- parse_qhet(
      qh = qh,
      rho = RHO_PRIMARY,
      outcome =
        aa$outcome,
      combo_id =
        aa$combo_id,
      analysis_set =
        "PRIMARY_RHO_QHET",
      n_snp =
        nrow(
          pp$data
        ),
      ci_requested = TRUE
    )

    rr[
      ,
      rho_source :=
        RHO_PRIMARY_SOURCE
    ]

    primary_list[[length(primary_list) + 1]] <- rr

    capture.output(
      qh,
      file = file.path(
        DIR_RES,
        paste0(
          "STEP3D_qhet_primary_",
          aa$outcome,
          "_rho_",
          gsub(
            "\\.",
            "p",
            format(
              RHO_PRIMARY,
              trim = TRUE
            )
          ),
          ".txt"
        )
      )
    )
  }

  primary_results <- rbindlist(
    primary_list,
    fill = TRUE
  )

  fwrite(
    primary_results,
    file.path(
      DIR_RES,
      "STEP3D_qhet_primary_rho_results.csv"
    )
  )

} else {

  writeLines(
    c(
      "No RHO_PRIMARY was supplied.",
      "",
      "Therefore STEP3D does NOT designate any single qhet_mvmr estimate as primary.",
      "Only rho-grid sensitivity results are produced.",
      "",
      "For manuscript-level use of a single qhet estimate, supply a defensible BMI–OSA PHENOTYPIC correlation from individual-level data or a suitable external study.",
      "Do NOT substitute LDSC genetic correlation."
    ),
    file.path(
      DIR_REP,
      "STEP3D_NO_PRIMARY_RHO.txt"
    )
  )
}


# ==============================================================================
# 15. RHO-SENSITIVITY FIGURES — SELECTED COMBO ONLY
# ==============================================================================

selected_grid <- grid_results[
  analysis_set ==
    "STEP3C_SELECTED" &
  status == "OK" &
  is.finite(OR)
]

# User preference: single figures, no combined panels, no internal legend.
# One plot for each exposure-outcome combination.

make_rho_plot <- function(
  exposure_name,
  outcome_name
) {

  d <- selected_grid[
    exposure ==
      exposure_name &
    outcome ==
      outcome_name
  ]

  if (
    nrow(d) == 0
  ) {

    return(
      invisible(
        NULL
      )
    )
  }

  # STEP3C conventional IVW estimate for reference
  ivw_row <- step3c_ivw[
    exposure ==
      exposure_name &
    outcome ==
      outcome_name
  ]

  ivw_or <- if (
    nrow(ivw_row) > 0
  ) {

    ivw_row$OR[1]

  } else {

    NA_real_
  }

  p <- ggplot(
    d,
    aes(
      x = rho,
      y = OR
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_line(
      linewidth = 0.75
    ) +
    geom_point(
      size = 2.1
    )

  if (
    is.finite(
      ivw_or
    )
  ) {

    p <- p +
      geom_hline(
        yintercept =
          ivw_or,
        linetype =
          "dotted",
        linewidth =
          0.50
      )
  }

  if (
    is.finite(
      RHO_PRIMARY
    )
  ) {

    p <- p +
      geom_vline(
        xintercept =
          RHO_PRIMARY,
        linetype =
          "dotdash",
        linewidth =
          0.50
      )
  }

  p <- p +
    labs(
      x =
        "Assumed BMI–OSA phenotypic correlation (rho)",
      y =
        paste0(
          "Q-minimization OR for ",
          exposure_name,
          " → ",
          outcome_name
        )
    ) +
    theme_classic(
      base_size = 11,
      base_family =
        BASE_FAMILY
    ) +
    theme(
      legend.position =
        "none",
      plot.margin =
        margin(
          8,
          12,
          8,
          8
        )
    )

  stub <- paste0(
    "STEP3D_qhet_rho_",
    exposure_name,
    "_to_",
    outcome_name
  )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(
        stub,
        ".tiff"
      )
    ),
    p,
    width = 6.5,
    height = 4.2,
    units = "in",
    dpi = FIG_DPI,
    compression =
      "lzw"
  )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(
        stub,
        ".pdf"
      )
    ),
    p,
    width = 6.5,
    height = 4.2,
    units = "in",
    device =
      cairo_pdf
  )

  fwrite(
    d,
    file.path(
      DIR_FIG,
      paste0(
        stub,
        "_source_data.csv"
      )
    )
  )
}

for (
  ex in c(
    "BMI",
    "OSA"
  )
) {

  for (
    oy in c(
      "HFpEF",
      "AF"
    )
  ) {

    make_rho_plot(
      exposure_name =
        ex,
      outcome_name =
        oy
    )
  }
}

writeLines(
  c(
    "STEP3D rho-sensitivity figure legend",
    "",
    "Each figure shows the Q-statistic-minimization MVMR effect estimate across a grid of assumed BMI–OSA phenotypic correlations.",
    "The solid curve represents qhet_mvmr estimates; the horizontal dashed line denotes OR=1.",
    "The horizontal dotted line, where present, denotes the corresponding conventional IVW MVMR point estimate from STEP3C.",
    "If a prespecified RHO_PRIMARY is supplied, a vertical dot-dash line marks that correlation.",
    "The rho grid is a sensitivity analysis and should not be interpreted as a set of empirically observed phenotypic correlations."
  ),
  file.path(
    DIR_FIG,
    "STEP3D_qhet_rho_figures_legend.txt"
  )
)


# ==============================================================================
# 16. OPTIONAL PRIMARY-RHO FOREST FIGURE
# ==============================================================================

if (
  nrow(
    primary_results
  ) > 0
) {

  d <- copy(
    primary_results
  )

  d[
    ,
    pair := paste0(
      exposure,
      " → ",
      outcome
    )
  ]

  d[
    ,
    pair := factor(
      pair,
      levels = rev(
        c(
          "BMI → HFpEF",
          "OSA → HFpEF",
          "BMI → AF",
          "OSA → AF"
        )
      )
    )
  ]

  # Only draw CIs if qhet returned them successfully.
  have_ci <- any(
    is.finite(
      d$OR_low
    ) &
    is.finite(
      d$OR_high
    )
  )

  p <- ggplot(
    d,
    aes(
      y = pair,
      x = OR
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45
    )

  if (
    have_ci
  ) {

    p <- p +
      geom_errorbarh(
        aes(
          xmin = OR_low,
          xmax = OR_high
        ),
        height = 0,
        linewidth = 0.70
      )
  }

  p <- p +
    geom_point(
      shape = 15,
      size = 3.2
    ) +
    scale_x_log10() +
    labs(
      x =
        "Weak-instrument robust direct-effect OR",
      y = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family =
        BASE_FAMILY
    ) +
    theme(
      legend.position = "none",
      axis.line.y =
        element_blank(),
      axis.ticks.y =
        element_blank()
    )

  ggsave(
    file.path(
      DIR_FIG,
      "STEP3D_qhet_PRIMARY_RHO_forest.tiff"
    ),
    p,
    width = 6.6,
    height = 4.2,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  ggsave(
    file.path(
      DIR_FIG,
      "STEP3D_qhet_PRIMARY_RHO_forest.pdf"
    ),
    p,
    width = 6.6,
    height = 4.2,
    units = "in",
    device = cairo_pdf
  )

  writeLines(
    c(
      "STEP3D primary-rho qhet MVMR legend",
      "",
      paste0(
        "Phenotypic correlation rho = ",
        RHO_PRIMARY,
        "."
      ),
      paste0(
        "Source: ",
        RHO_PRIMARY_SOURCE
      ),
      "Squares indicate qhet_mvmr Q-statistic-minimization point estimates; horizontal lines indicate bootstrap 95% confidence intervals where returned by the installed MVMR package.",
      "These estimates are intended as weak-instrument/pleiotropy-robust sensitivity analyses and should be interpreted alongside conventional MVMR conditional F-statistics and UVMR results."
    ),
    file.path(
      DIR_FIG,
      "STEP3D_qhet_PRIMARY_RHO_forest_legend.txt"
    )
  )
}


# ==============================================================================
# 17. AUTOMATED DECISION SUMMARY
# ==============================================================================

selected_stability <- stability[
  analysis_set ==
    "STEP3C_SELECTED"
]

decision <- c(
  "STEP 3D FINAL DECISION",
  "",
  "Purpose:",
  "Weak-instrument robust MVMR sensitivity analysis using Q-statistic minimisation (qhet_mvmr).",
  "",
  "STEP3C conventional instrument-strength context:"
)

for (
  i in seq_len(
    nrow(
      selected
    )
  )
) {

  decision <- c(
    decision,
    paste0(
      "- ",
      selected$outcome[i],
      ": ",
      selected$combo_id[i],
      "; F_BMI=",
      signif(
        selected$F_BMI[i],
        5
      ),
      "; F_OSA=",
      signif(
        selected$F_OSA[i],
        5
      ),
      "; both F>=10 = ",
      selected$F_BMI[i] >=
        F_REFERENCE &&
        selected$F_OSA[i] >=
        F_REFERENCE
    )
  )
}

decision <- c(
  decision,
  "",
  "qhet rho-grid direction stability:"
)

for (
  i in seq_len(
    nrow(
      selected_stability
    )
  )
) {

  decision <- c(
    decision,
    paste0(
      "- ",
      selected_stability$exposure[i],
      " -> ",
      selected_stability$outcome[i],
      ": ",
      selected_stability$stability_class[i],
      "; OR range=",
      signif(
        selected_stability$OR_min[i],
        4
      ),
      " to ",
      signif(
        selected_stability$OR_max[i],
        4
      ),
      " across rho ",
      selected_stability$rho_min[i],
      " to ",
      selected_stability$rho_max[i]
    )
  )
}

decision <- c(
  decision,
  ""
)

if (
  is.finite(
    RHO_PRIMARY
  )
) {

  decision <- c(
    decision,
    paste0(
      "Prespecified phenotypic rho = ",
      RHO_PRIMARY,
      "."
    ),
    paste0(
      "Source = ",
      RHO_PRIMARY_SOURCE
    ),
    "Bootstrap qhet estimates were generated at this prespecified rho."
  )

} else {

  decision <- c(
    decision,
    "No externally justified RHO_PRIMARY was supplied.",
    "Therefore no single qhet estimate should be presented as a definitive primary MVMR estimate.",
    "Use the rho-grid results only as robustness/sensitivity evidence."
  )
}

decision <- c(
  decision,
  "",
  "Interpretation rule:",
  "- Direction stable across a broad rho grid strengthens the sensitivity argument but does not repair the conventional conditional-F weakness.",
  "- A rho-sensitive direction means qhet conclusions depend materially on the assumed phenotypic correlation.",
  "- Conventional STEP3C IVW MVMR should remain exploratory when conditional F<10.",
  "- STEP3D must not be used to select a rho because it produces a more significant outcome association."
)

writeLines(
  decision,
  file.path(
    DIR_REP,
    "STEP3D_FINAL_DECISION.txt"
  )
)


# ==============================================================================
# 18. METHODS NOTE
# ==============================================================================

methods_note <- c(
  "STEP3D METHODS NOTE",
  "",
  "Because conventional MVMR showed inadequate conditional instrument strength for OSA, Q-statistic-minimization MVMR was performed as a weak-instrument/pleiotropy-robust sensitivity analysis.",
  "The qhet_mvmr method requires a phenotypic correlation matrix between the exposures.",
  "When no externally justified BMI–OSA phenotypic correlation was available, estimates were evaluated across a prespecified correlation grid rather than selecting a correlation from the outcome results.",
  "Genetic correlation estimates from LDSC were not substituted for phenotypic correlation.",
  "The genome-wide-significant BMI instrument set and STEP3C-selected OSA instrument tier were not further relaxed in STEP3D.",
  "Results from STEP3D are interpreted as sensitivity evidence and do not override the conventional conditional-F assessment."
)

writeLines(
  methods_note,
  file.path(
    DIR_REP,
    "STEP3D_METHODS_NOTE.txt"
  )
)


# ==============================================================================
# 19. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_REP,
    "STEP3D_sessionInfo.txt"
  )
)


# ==============================================================================
# 20. FINAL CONSOLE
# ==============================================================================

log_msg(
  "============================================================"
)

log_msg(
  "STEP 3D COMPLETE"
)

log_msg(
  "============================================================"
)

log_msg(
  "Output: ",
  OUT_DIR
)

log_msg(
  ""
)

log_msg(
  "Please send back:"
)

log_msg(
  "1) 01_results/STEP3D_qhet_rho_grid_results.csv"
)

log_msg(
  "2) 01_results/STEP3D_qhet_stability_summary.csv"
)

log_msg(
  "3) 01_results/STEP3D_IVW_vs_qhet_comparison.csv"
)

log_msg(
  "4) 03_reports/STEP3D_FINAL_DECISION.txt"
)

log_msg(
  "5) 03_reports/STEP3D_master_console.log"
)

if (
  is.finite(
    RHO_PRIMARY
  )
) {

  log_msg(
    "6) 01_results/STEP3D_qhet_primary_rho_results.csv"
  )
}

log_msg(
  ""
)

log_msg(
  "Figures are in 02_figures."
)
