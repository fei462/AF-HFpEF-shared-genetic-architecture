# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP8C_V4_1_SER_RELIABILITY_AUDIT_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP8C-V4.1 — SER RELIABILITY AUDIT (manifest lookup fix)
# Project: AF / HFpEF / BMI / OSA
#
# PURPOSE
# -------
# Audit SuSiE-RSS fits from STEP8C-V3 for which
# R_reliability_flag == TRUE.
#
# IMPORTANT:
#   - NO re-fitting.
#   - NO region/SNP/prior changes.
#   - NO new LD calculation.
#   - Uses the official diagnostic model already stored in:
#
#       fit$R_finite_diagnostics$ser_model
#
# The current susieR mismatch-aware workflow recommends comparing
# a flagged EB fit with this one-effect SER diagnostic model and
# interpreting the result conservatively.
#
# This script performs FOUR audits:
#
#   A) EB versus SER top SNP / PIP comparison.
#   B) EB credible set versus SER credible set overlap.
#   C) Variant-level EB PIP versus SER PIP comparison.
#   D) For each flagged trait, compare its SER credible set with
#      the PARTNER phenotype's EB credible set(s) from the same
#      coloc test.
#
# The partner comparison is useful for asking whether the safer
# SER fallback still localizes to the same signal cluster that
# supported cross-trait colocalization.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"

V3_ROOT <- file.path(
  DATA_ROOT,
  "STEP8_COLOC",
  "07_coloc_susie_1000G_EUR_mismatch_aware"
)

FIT_DIR <- file.path(
  V3_ROOT,
  "01_susie_fits"
)

V3_TABLE_DIR <- file.path(
  V3_ROOT,
  "03_tables"
)

AUDIT_ROOT <- file.path(
  V3_ROOT,
  "05_SER_reliability_audit"
)

AUDIT_TABLE_DIR <- file.path(
  AUDIT_ROOT,
  "01_tables"
)

AUDIT_VARIANT_DIR <- file.path(
  AUDIT_ROOT,
  "02_variant_level"
)

AUDIT_LOG_DIR <- file.path(
  AUDIT_ROOT,
  "03_logs"
)

for (d in c(
  AUDIT_ROOT,
  AUDIT_TABLE_DIR,
  AUDIT_VARIANT_DIR,
  AUDIT_LOG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 1. PACKAGE CHECK
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)

if (!requireNamespace("susieR", quietly = TRUE)) {
  stop(
    "susieR is not installed. STEP8C-V4 only reads existing ",
    "V3 fit objects, but the package namespace is still required."
  )
}

SUSIER_VERSION <- as.character(
  packageVersion("susieR")
)

cat(
  "\nsusieR version: ",
  SUSIER_VERSION,
  "\n",
  sep = ""
)

# ============================================================
# 2. LOCATE REQUIRED V3 TABLES ROBUSTLY
# ============================================================

locate_one <- function(
  dir,
  exact_name = NULL,
  regex = NULL,
  label = "file"
) {

  if (
    !is.null(exact_name)
  ) {

    f <- file.path(
      dir,
      exact_name
    )

    if (file.exists(f)) {
      return(
        normalizePath(
          f,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  if (
    !is.null(regex)
  ) {

    hits <- list.files(
      dir,
      pattern = regex,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(hits) == 1L) {
      return(
        normalizePath(
          hits,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }

    if (length(hits) > 1L) {
      stop(
        "Multiple candidate ",
        label,
        " files found:\n",
        paste(
          hits,
          collapse = "\n"
        )
      )
    }
  }

  stop(
    "Could not locate ",
    label,
    " under:\n",
    dir
  )
}

RELIABILITY_FILE <- locate_one(
  V3_TABLE_DIR,
  exact_name = "STEP8C_V3_susie_reliability_QC.csv",
  regex = "V3.*susie.*reliability.*QC.*\\.csv$",
  label = "V3 reliability QC"
)

MANIFEST_FILE <- locate_one(
  V3_TABLE_DIR,
  exact_name = "STEP8C_V3_test_manifest.csv",
  regex = "V3.*test.*manifest.*\\.csv$",
  label = "V3 test manifest"
)

FINEMAP_FILE <- locate_one(
  V3_TABLE_DIR,
  exact_name = "STEP8C_V3_susie_finemap_summary.csv",
  regex = "V3.*susie.*finemap.*summary.*\\.csv$",
  label = "V3 finemap summary"
)

if (!dir.exists(FIT_DIR)) {
  stop(
    "V3 SuSiE fit directory does not exist:\n",
    FIT_DIR
  )
}

# ============================================================
# 3. READ V3 RESULTS
# ============================================================

REL <- fread(
  RELIABILITY_FILE
)

MANIFEST <- fread(
  MANIFEST_FILE
)

FINEMAP_V3 <- fread(
  FINEMAP_FILE
)

required_rel <- c(
  "test_id",
  "trait",
  "converged",
  "R_sensitivity_flag",
  "R_reliability_flag"
)

if (!all(
  required_rel %in%
  names(REL)
)) {
  stop(
    "Reliability QC lacks columns: ",
    paste(
      setdiff(
        required_rel,
        names(REL)
      ),
      collapse = ", "
    )
  )
}

required_manifest <- c(
  "test_id",
  "trait1",
  "trait2"
)

if (!all(
  required_manifest %in%
  names(MANIFEST)
)) {
  stop(
    "Test manifest lacks columns: ",
    paste(
      setdiff(
        required_manifest,
        names(MANIFEST)
      ),
      collapse = ", "
    )
  )
}

# Hard QC: each test_id must occur exactly once in the V3 manifest.
manifest_ids_qc <- as.character(
  MANIFEST[["test_id"]]
)

if (
  anyNA(manifest_ids_qc) ||
  any(!nzchar(manifest_ids_qc))
) {
  stop(
    "V3 test manifest contains missing/blank test_id values."
  )
}

dup_manifest_ids <- unique(
  manifest_ids_qc[
    duplicated(manifest_ids_qc)
  ]
)

if (length(dup_manifest_ids)) {
  stop(
    "Duplicated test_id values found in V3 manifest: ",
    paste(
      dup_manifest_ids,
      collapse = ";"
    )
  )
}

cat(
  "\nManifest QC PASS: ",
  length(manifest_ids_qc),
  " unique test IDs.\n",
  sep = ""
)

# Robust TRUE/FALSE parser for CSV values.
as_flag <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in% c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
  ] <- TRUE

  out[
    y %in% c(
      "FALSE",
      "F",
      "0",
      "NO",
      "N"
    )
  ] <- FALSE

  out
}

REL[
  ,
  R_reliability_flag_bool :=
    as_flag(
      R_reliability_flag
    )
]

REL[
  ,
  R_sensitivity_flag_bool :=
    as_flag(
      R_sensitivity_flag
    )
]

FLAGGED <- REL[
  R_reliability_flag_bool == TRUE
]

if (!nrow(FLAGGED)) {

  out <- file.path(
    AUDIT_TABLE_DIR,
    "STEP8C_V4_no_flagged_fits.txt"
  )

  writeLines(
    c(
      "No STEP8C-V3 fit had R_reliability_flag == TRUE.",
      "SER reliability audit is therefore not required."
    ),
    out
  )

  cat(
    "\nNo reliability-flagged V3 fits were found.\n"
  )

  quit(
    save = "no",
    status = 0
  )
}

cat(
  "\nReliability-flagged V3 fits:\n"
)

print(
  FLAGGED[
    ,
    .(
      test_id,
      trait,
      R_sensitivity_flag,
      R_reliability_flag
    )
  ]
)

# ============================================================
# 4. HELPERS
# ============================================================

safe_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_\\-]+",
    "_",
    x
  )
}

locate_fit <- function(
  test_id,
  trait
) {

  exact <- file.path(
    FIT_DIR,
    paste0(
      safe_filename(test_id),
      "__",
      trait,
      "_susie.rds"
    )
  )

  if (file.exists(exact)) {
    return(
      normalizePath(
        exact,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  # Fallback if a slightly different filename was used.
  files <- list.files(
    FIT_DIR,
    pattern = "\\.rds$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  tid_safe <- safe_filename(test_id)

  hits <- files[
    grepl(
      tid_safe,
      basename(files),
      fixed = TRUE
    ) &
    grepl(
      paste0(
        "__",
        trait,
        "_susie"
      ),
      basename(files),
      fixed = TRUE
    )
  ]

  if (length(hits) != 1L) {
    stop(
      "Could not uniquely locate fit for:\n",
      "test_id = ", test_id, "\n",
      "trait   = ", trait, "\n",
      "Candidates:\n",
      paste(
        hits,
        collapse = "\n"
      )
    )
  }

  normalizePath(
    hits,
    winslash = "/",
    mustWork = TRUE
  )
}

infer_snp_names <- function(
  obj,
  fallback = NULL
) {

  n <- NULL

  if (!is.null(obj$pip)) {
    n <- length(
      obj$pip
    )
  } else if (!is.null(obj$alpha)) {
    if (is.matrix(obj$alpha)) {
      n <- ncol(
        obj$alpha
      )
    } else {
      n <- length(
        obj$alpha
      )
    }
  }

  if (is.null(n)) {
    stop(
      "Cannot infer number of variants from fit object."
    )
  }

  candidates <- list(
    if (!is.null(obj$pip)) {
      names(obj$pip)
    } else NULL,
    if (
      !is.null(obj$alpha) &&
      is.matrix(obj$alpha)
    ) {
      colnames(obj$alpha)
    } else NULL,
    if (!is.null(obj$variable_names)) {
      as.character(
        obj$variable_names
      )
    } else NULL,
    fallback
  )

  for (x in candidates) {

    if (
      !is.null(x) &&
      length(x) == n &&
      all(
        !is.na(x)
      ) &&
      all(
        nzchar(
          as.character(x)
        )
      )
    ) {
      return(
        as.character(x)
      )
    }
  }

  paste0(
    "SNP_",
    seq_len(n)
  )
}

extract_pip_vector <- function(
  obj,
  fallback_names = NULL
) {

  if (!is.null(obj$pip)) {

    pip <- as.numeric(
      obj$pip
    )

  } else if (!is.null(obj$alpha)) {

    if (is.matrix(obj$alpha)) {

      # SER has one alpha row. If more than one row exists, use
      # the standard marginal inclusion probability formula.
      if (nrow(obj$alpha) == 1L) {

        pip <- as.numeric(
          obj$alpha[1, ]
        )

      } else {

        pip <- 1 -
          apply(
            1 - obj$alpha,
            2,
            prod
          )
      }

    } else {

      pip <- as.numeric(
        obj$alpha
      )
    }

  } else {

    stop(
      "Object contains neither pip nor alpha."
    )
  }

  snps <- infer_snp_names(
    obj,
    fallback = fallback_names
  )

  # Some SER implementations may append a null variable.
  # Remove it when null_index points to an extra element.
  if (
    !is.null(obj$null_index)
  ) {

    ni <- suppressWarnings(
      as.integer(
        obj$null_index[1]
      )
    )

    if (
      is.finite(ni) &&
      ni >= 1L &&
      ni <= length(pip) &&
      !is.null(fallback_names) &&
      length(pip) ==
        length(fallback_names) + 1L
    ) {

      pip <- pip[-ni]
      snps <- snps[-ni]
    }
  }

  if (
    length(pip) !=
    length(snps)
  ) {
    stop(
      "PIP/SNP length mismatch."
    )
  }

  data.table(
    SNP = as.character(
      snps
    ),
    PIP = as.numeric(
      pip
    )
  )
}

extract_main_cs <- function(
  fit,
  snp_names
) {

  cs_list <- fit$sets$cs

  if (
    is.null(cs_list) ||
    !length(cs_list)
  ) {
    return(
      list()
    )
  }

  out <- list()

  for (
    k in seq_along(
      cs_list
    )
  ) {

    idx <- as.integer(
      cs_list[[k]]
    )

    idx <- idx[
      is.finite(idx) &
      idx >= 1L &
      idx <= length(
        snp_names
      )
    ]

    out[[
      paste0(
        "CS",
        k
      )
    ]] <- unique(
      snp_names[
        idx
      ]
    )
  }

  out
}

extract_ser_cs <- function(
  ser,
  ser_pip_dt,
  coverage = 0.95
) {

  # Prefer the official SER credible set saved by susieR.
  if (
    !is.null(
      ser$sets$cs
    ) &&
    length(
      ser$sets$cs
    ) >= 1L
  ) {

    snp_names <- ser_pip_dt$SNP

    idx <- as.integer(
      ser$sets$cs[[1]]
    )

    idx <- idx[
      is.finite(idx) &
      idx >= 1L &
      idx <= length(
        snp_names
      )
    ]

    if (length(idx)) {

      return(
        list(
          snps = unique(
            snp_names[
              idx
            ]
          ),
          method = "official_ser_sets_cs",
          achieved_pip = sum(
            ser_pip_dt[
              SNP %in%
                snp_names[
                  idx
                ],
              PIP
            ],
            na.rm = TRUE
          )
        )
      )
    }
  }

  # Fallback: construct the shortest cumulative-PIP set reaching
  # 95% using the SER PIP vector.
  x <- copy(
    ser_pip_dt
  )

  x <- x[
    is.finite(PIP) &
    PIP >= 0
  ]

  setorder(
    x,
    -PIP
  )

  if (!nrow(x)) {
    return(
      list(
        snps = character(0),
        method = "no_valid_ser_pip",
        achieved_pip = NA_real_
      )
    )
  }

  x[
    ,
    cumPIP := cumsum(
      PIP
    )
  ]

  idx_end <- which(
    x$cumPIP >= coverage
  )[1]

  if (
    is.na(idx_end)
  ) {
    idx_end <- nrow(x)
  }

  keep <- seq_len(
    idx_end
  )

  list(
    snps = x$SNP[
      keep
    ],
    method = paste0(
      "constructed_cumulative_",
      coverage
    ),
    achieved_pip = sum(
      x$PIP[
        keep
      ],
      na.rm = TRUE
    )
  )
}

set_overlap_metrics <- function(
  A,
  B
) {

  A <- unique(
    as.character(A)
  )

  B <- unique(
    as.character(B)
  )

  A <- A[
    !is.na(A) &
    nzchar(A)
  ]

  B <- B[
    !is.na(B) &
    nzchar(B)
  ]

  inter <- intersect(
    A,
    B
  )

  uni <- union(
    A,
    B
  )

  nA <- length(A)
  nB <- length(B)
  nI <- length(inter)

  data.table(
    n_A = nA,
    n_B = nB,
    n_overlap = nI,
    jaccard = if (
      length(uni)
    ) {
      nI /
        length(uni)
    } else {
      NA_real_
    },
    overlap_coefficient = if (
      min(nA, nB) > 0
    ) {
      nI /
        min(nA, nB)
    } else {
      NA_real_
    },
    A_fraction_overlapped = if (
      nA > 0
    ) {
      nI / nA
    } else {
      NA_real_
    },
    B_fraction_overlapped = if (
      nB > 0
    ) {
      nI / nB
    } else {
      NA_real_
    },
    overlap_snps = paste(
      inter,
      collapse = ";"
    )
  )
}

get_diag_scalar <- function(
  fit,
  name
) {

  d <- fit$R_finite_diagnostics

  if (
    is.null(d) ||
    is.null(
      d[[name]]
    )
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    as.numeric(
      d[[name]][1]
    )
  )
}

get_diag_flag <- function(
  fit,
  name
) {

  d <- fit$R_finite_diagnostics

  if (
    is.null(d) ||
    is.null(
      d[[name]]
    )
  ) {
    return(
      NA
    )
  }

  as.logical(
    d[[name]][1]
  )
}

get_partner_trait <- function(
  test_id_value,
  trait_value
) {

  # IMPORTANT:
  # Do NOT use data.table expressions such as:
  #   MANIFEST[test_id == test_id_value]
  # when an external variable may share a name with a column.
  # Use base-vector indexing explicitly so there is no scope ambiguity.
  manifest_ids <- as.character(
    MANIFEST[["test_id"]]
  )

  target_id <- as.character(
    test_id_value
  )

  idx <- which(
    manifest_ids == target_id
  )

  if (length(idx) != 1L) {
    stop(
      "Manifest lookup failed for test_id: ",
      target_id,
      " | matches found = ",
      length(idx),
      " | manifest IDs = ",
      paste(
        manifest_ids,
        collapse = ";"
      )
    )
  }

  z <- MANIFEST[
    idx
  ]

  t1 <- as.character(
    z[["trait1"]][1]
  )

  t2 <- as.character(
    z[["trait2"]][1]
  )

  target_trait <- as.character(
    trait_value
  )

  if (identical(target_trait, t1)) {
    return(t2)
  }

  if (identical(target_trait, t2)) {
    return(t1)
  }

  stop(
    "Flagged trait ",
    target_trait,
    " is not trait1/trait2 for test ",
    target_id,
    " | trait1=",
    t1,
    " | trait2=",
    t2
  )
}

# ============================================================
# 5. AUDIT ALL FLAGGED FITS
# ============================================================

SUMMARY_LIST <- list()
OVERLAP_LIST <- list()
VARIANT_LIST <- list()
PARTNER_LIST <- list()
SER_CS_LIST <- list()
QC_LIST <- list()

for (
  i in seq_len(
    nrow(FLAGGED)
  )
) {

  rr <- FLAGGED[i]

  test_id_i <- as.character(
    rr$test_id
  )

  trait_i <- as.character(
    rr$trait
  )

  cat(
    "\n====================================================\n",
    "[", i, "/", nrow(FLAGGED), "] ",
    test_id_i,
    " | ",
    trait_i,
    "\n====================================================\n",
    sep = ""
  )

  fit_file <- tryCatch(
    locate_fit(
      test_id_i,
      trait_i
    ),
    error = function(e) {
      e
    }
  )

  if (
    inherits(
      fit_file,
      "error"
    )
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      trait = trait_i,
      status = "FIT_NOT_FOUND",
      message = conditionMessage(
        fit_file
      )
    )

    next
  }

  fit <- tryCatch(
    readRDS(
      fit_file
    ),
    error = function(e) e
  )

  if (
    inherits(
      fit,
      "error"
    )
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      trait = trait_i,
      status = "FIT_READ_FAIL",
      message = conditionMessage(
        fit
      )
    )

    next
  }

  diag <- fit$R_finite_diagnostics

  if (
    is.null(diag)
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      trait = trait_i,
      status = "NO_R_FINITE_DIAGNOSTICS",
      message = ""
    )

    next
  }

  ser <- diag$ser_model

  if (
    is.null(ser)
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      trait = trait_i,
      status = "NO_SER_MODEL",
      message = paste(
        names(diag),
        collapse = ";"
      )
    )

    next
  }

  # ----------------------------------------------------------
  # Main EB model
  # ----------------------------------------------------------

  main_pip <- extract_pip_vector(
    fit
  )

  main_snps <- main_pip$SNP

  main_cs <- extract_main_cs(
    fit,
    main_snps
  )

  main_pip_ord <- copy(
    main_pip
  )

  setorder(
    main_pip_ord,
    -PIP
  )

  main_top_snp <- if (
    nrow(main_pip_ord)
  ) {
    main_pip_ord$SNP[1]
  } else {
    NA_character_
  }

  main_top_pip <- if (
    nrow(main_pip_ord)
  ) {
    main_pip_ord$PIP[1]
  } else {
    NA_real_
  }

  # ----------------------------------------------------------
  # SER diagnostic model
  # ----------------------------------------------------------

  ser_pip <- extract_pip_vector(
    ser,
    fallback_names = main_snps
  )

  # Force exact variant matching/order to the main EB fit.
  if (
    !setequal(
      ser_pip$SNP,
      main_snps
    )
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      trait = trait_i,
      status = "SER_VARIANT_SET_MISMATCH",
      message = paste0(
        "main_n=",
        length(main_snps),
        "; ser_n=",
        nrow(ser_pip)
      )
    )

    next
  }

  ser_pip <- ser_pip[
    match(
      main_snps,
      SNP
    )
  ]

  if (
    !identical(
      ser_pip$SNP,
      main_snps
    )
  ) {
    stop(
      "Internal SER/main SNP-order reconstruction failed for ",
      test_id_i,
      " | ",
      trait_i
    )
  }

  ser_pip_ord <- copy(
    ser_pip
  )

  setorder(
    ser_pip_ord,
    -PIP
  )

  ser_top_snp <- if (
    nrow(ser_pip_ord)
  ) {
    ser_pip_ord$SNP[1]
  } else {
    NA_character_
  }

  ser_top_pip <- if (
    nrow(ser_pip_ord)
  ) {
    ser_pip_ord$PIP[1]
  } else {
    NA_real_
  }

  ser_cs_info <- extract_ser_cs(
    ser,
    ser_pip,
    coverage = 0.95
  )

  ser_cs <- ser_cs_info$snps

  SER_CS_LIST[[
    length(SER_CS_LIST) + 1L
  ]] <- data.table(
    test_id = test_id_i,
    trait = trait_i,
    ser_cs_method =
      ser_cs_info$method,
    ser_cs_size =
      length(ser_cs),
    ser_cs_achieved_PIP =
      ser_cs_info$achieved_pip,
    ser_top_snp =
      ser_top_snp,
    ser_top_PIP =
      ser_top_pip,
    ser_cs_snps =
      paste(
        ser_cs,
        collapse = ";"
      )
  )

  # ----------------------------------------------------------
  # Variant-level comparison
  # ----------------------------------------------------------

  var_dt <- merge(
    main_pip[
      ,
      .(
        SNP,
        EB_PIP = PIP
      )
    ],
    ser_pip[
      ,
      .(
        SNP,
        SER_PIP = PIP
      )
    ],
    by = "SNP",
    all = TRUE
  )

  main_any_cs <- unique(
    unlist(
      main_cs,
      use.names = FALSE
    )
  )

  var_dt[
    ,
    `:=`(
      test_id = test_id_i,
      trait = trait_i,
      in_EB_any_CS =
        SNP %in%
        main_any_cs,
      in_SER_CS =
        SNP %in%
        ser_cs,
      delta_PIP =
        SER_PIP -
        EB_PIP
    )
  ]

  VARIANT_LIST[[
    length(VARIANT_LIST) + 1L
  ]] <- var_dt

  # ----------------------------------------------------------
  # EB CS versus SER CS overlap
  # ----------------------------------------------------------

  if (length(main_cs)) {

    fit_overlap_rows <- list()

    for (
      cs_name in names(
        main_cs
      )
    ) {

      eb_cs <- main_cs[[
        cs_name
      ]]

      mm <- set_overlap_metrics(
        eb_cs,
        ser_cs
      )

      mm[
        ,
        `:=`(
          test_id = test_id_i,
          trait = trait_i,
          EB_CS = cs_name,
          EB_CS_size =
            length(
              eb_cs
            ),
          SER_CS_size =
            length(
              ser_cs
            ),
          EB_lead_snp = {
            z <- main_pip[
              SNP %in%
                eb_cs
            ]
            if (
              nrow(z)
            ) {
              z$SNP[
                which.max(
                  z$PIP
                )
              ]
            } else {
              NA_character_
            }
          },
          SER_top_snp =
            ser_top_snp,
          EB_lead_in_SER_CS = {
            lead <- main_pip[
              SNP %in%
                eb_cs
            ]
            if (
              nrow(lead)
            ) {
              lead_snp <- lead$SNP[
                which.max(
                  lead$PIP
                )
              ]
              lead_snp %in%
                ser_cs
            } else {
              NA
            }
          },
          SER_top_in_EB_CS =
            ser_top_snp %in%
            eb_cs
        )
      ]

      fit_overlap_rows[[
        length(fit_overlap_rows) + 1L
      ]] <- mm
    }

    fit_overlap <- rbindlist(
      fit_overlap_rows,
      fill = TRUE
    )

    OVERLAP_LIST[[
      length(OVERLAP_LIST) + 1L
    ]] <- fit_overlap

    best_i <- which.max(
      fifelse(
        is.na(
          fit_overlap$overlap_coefficient
        ),
        -Inf,
        fit_overlap$overlap_coefficient
      )
    )

    best_overlap <- fit_overlap[
      best_i
    ]

  } else {

    best_overlap <- data.table(
      EB_CS = NA_character_,
      n_overlap = 0L,
      jaccard = NA_real_,
      overlap_coefficient =
        NA_real_,
      A_fraction_overlapped =
        NA_real_,
      B_fraction_overlapped =
        NA_real_
    )
  }

  # ----------------------------------------------------------
  # Partner-trait EB CS versus flagged SER CS
  # ----------------------------------------------------------

  partner_trait <- get_partner_trait(
    test_id_i,
    trait_i
  )

  partner_fit_file <- tryCatch(
    locate_fit(
      test_id_i,
      partner_trait
    ),
    error = function(e) e
  )

  partner_available <- !inherits(
    partner_fit_file,
    "error"
  )

  max_partner_overlap <- NA_real_
  max_partner_jaccard <- NA_real_
  best_partner_cs <- NA_character_
  ser_top_in_partner_any_cs <- NA

  if (partner_available) {

    partner_fit <- tryCatch(
      readRDS(
        partner_fit_file
      ),
      error = function(e) e
    )

    if (
      !inherits(
        partner_fit,
        "error"
      )
    ) {

      partner_pip <- extract_pip_vector(
        partner_fit
      )

      partner_cs <- extract_main_cs(
        partner_fit,
        partner_pip$SNP
      )

      if (length(partner_cs)) {

        partner_rows <- list()

        for (
          pcs in names(
            partner_cs
          )
        ) {

          pset <- partner_cs[[
            pcs
          ]]

          pm <- set_overlap_metrics(
            ser_cs,
            pset
          )

          pp <- partner_pip[
            SNP %in%
              pset
          ]

          partner_lead <- if (
            nrow(pp)
          ) {
            pp$SNP[
              which.max(
                pp$PIP
              )
            ]
          } else {
            NA_character_
          }

          pm[
            ,
            `:=`(
              test_id = test_id_i,
              flagged_trait =
                trait_i,
              partner_trait =
                partner_trait,
              partner_EB_CS =
                pcs,
              flagged_SER_CS_size =
                length(
                  ser_cs
                ),
              partner_EB_CS_size =
                length(
                  pset
                ),
              flagged_SER_top_snp =
                ser_top_snp,
              partner_EB_lead_snp =
                partner_lead,
              SER_top_in_partner_EB_CS =
                ser_top_snp %in%
                pset,
              partner_lead_in_SER_CS =
                partner_lead %in%
                ser_cs
            )
          ]

          partner_rows[[
            length(partner_rows) + 1L
          ]] <- pm
        }

        pdt <- rbindlist(
          partner_rows,
          fill = TRUE
        )

        PARTNER_LIST[[
          length(PARTNER_LIST) + 1L
        ]] <- pdt

        k <- which.max(
          fifelse(
            is.na(
              pdt$overlap_coefficient
            ),
            -Inf,
            pdt$overlap_coefficient
          )
        )

        max_partner_overlap <-
          pdt$overlap_coefficient[k]

        max_partner_jaccard <-
          pdt$jaccard[k]

        best_partner_cs <-
          pdt$partner_EB_CS[k]

        partner_any <- unique(
          unlist(
            partner_cs,
            use.names = FALSE
          )
        )

        ser_top_in_partner_any_cs <-
          ser_top_snp %in%
          partner_any
      }
    }
  }

  # ----------------------------------------------------------
  # Summary row
  # ----------------------------------------------------------

  spearman_pip <- suppressWarnings(
    cor(
      var_dt$EB_PIP,
      var_dt$SER_PIP,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )

  pearson_pip <- suppressWarnings(
    cor(
      var_dt$EB_PIP,
      var_dt$SER_PIP,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )

  SUMMARY_LIST[[
    length(SUMMARY_LIST) + 1L
  ]] <- data.table(
    test_id = test_id_i,
    trait = trait_i,
    partner_trait =
      partner_trait,

    V3_R_sensitivity_flag =
      get_diag_flag(
        fit,
        "R_sensitivity_flag"
      ),

    V3_R_reliability_flag =
      get_diag_flag(
        fit,
        "R_reliability_flag"
      ),

    V3_Q_art =
      get_diag_scalar(
        fit,
        "Q_art"
      ),

    V3_r_over_B =
      get_diag_scalar(
        fit,
        "r_over_B"
      ),

    V3_B_corrected =
      get_diag_scalar(
        fit,
        "B_corrected"
      ),

    V3_lambda_bias =
      get_diag_scalar(
        fit,
        "lambda_bias"
      ),

    EB_converged =
      isTRUE(
        fit$converged
      ),

    EB_n_CS =
      length(
        main_cs
      ),

    EB_top_snp =
      main_top_snp,

    EB_top_PIP =
      main_top_pip,

    SER_converged =
      if (
        !is.null(
          ser$converged
        )
      ) {
        isTRUE(
          ser$converged
        )
      } else {
        NA
      },

    SER_top_snp =
      ser_top_snp,

    SER_top_PIP =
      ser_top_pip,

    top_SNP_identical =
      identical(
        main_top_snp,
        ser_top_snp
      ),

    SER_CS_method =
      ser_cs_info$method,

    SER_CS_size =
      length(
        ser_cs
      ),

    SER_CS_achieved_PIP =
      ser_cs_info$achieved_pip,

    EB_top_in_SER_CS =
      main_top_snp %in%
      ser_cs,

    SER_top_in_EB_any_CS =
      ser_top_snp %in%
      main_any_cs,

    best_EB_CS =
      best_overlap$EB_CS[1],

    EB_SER_max_overlap_n =
      best_overlap$n_overlap[1],

    EB_SER_max_Jaccard =
      best_overlap$jaccard[1],

    EB_SER_max_overlap_coefficient =
      best_overlap$overlap_coefficient[1],

    EB_SER_max_EB_fraction =
      best_overlap$A_fraction_overlapped[1],

    EB_SER_max_SER_fraction =
      best_overlap$B_fraction_overlapped[1],

    PIP_Spearman =
      spearman_pip,

    PIP_Pearson =
      pearson_pip,

    partner_fit_available =
      partner_available,

    best_partner_EB_CS =
      best_partner_cs,

    SER_partner_max_Jaccard =
      max_partner_jaccard,

    SER_partner_max_overlap_coefficient =
      max_partner_overlap,

    SER_top_in_partner_any_EB_CS =
      ser_top_in_partner_any_cs
  )

  QC_LIST[[
    length(QC_LIST) + 1L
  ]] <- data.table(
    test_id = test_id_i,
    trait = trait_i,
    status = "PASS",
    message = paste0(
      "SER model found; EB_CS=",
      length(main_cs),
      "; SER_CS=",
      length(ser_cs),
      "; partner=",
      partner_trait
    )
  )
}

# ============================================================
# 6. COMBINE OUTPUTS
# ============================================================

SUMMARY <- if (
  length(SUMMARY_LIST)
) {
  rbindlist(
    SUMMARY_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

OVERLAP <- if (
  length(OVERLAP_LIST)
) {
  rbindlist(
    OVERLAP_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

VARIANT <- if (
  length(VARIANT_LIST)
) {
  rbindlist(
    VARIANT_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

PARTNER <- if (
  length(PARTNER_LIST)
) {
  rbindlist(
    PARTNER_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

SER_CS <- if (
  length(SER_CS_LIST)
) {
  rbindlist(
    SER_CS_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

QC <- if (
  length(QC_LIST)
) {
  rbindlist(
    QC_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

# ============================================================
# 7. ADD A TRANSPARENT DESCRIPTIVE AUDIT LABEL
#
# This is NOT an official susieR classification.
# It summarizes only EB-vs-SER signal-cluster agreement.
#
# Priority:
#   strong_cluster_agreement:
#      same top SNP OR overlap coefficient >= 0.80
#
#   partial_cluster_agreement:
#      nonzero overlap but <0.80
#
#   discordant_cluster:
#      no CS overlap
#
#   no_EB_CS:
#      EB fit itself reported no CS
# ============================================================

if (nrow(SUMMARY)) {

  SUMMARY[
    ,
    audit_label :=
      fifelse(
        EB_n_CS == 0,
        "no_EB_CS",
        fifelse(
          top_SNP_identical == TRUE |
          (
            is.finite(
              EB_SER_max_overlap_coefficient
            ) &
            EB_SER_max_overlap_coefficient >= 0.80
          ),
          "strong_cluster_agreement",
          fifelse(
            is.finite(
              EB_SER_max_overlap_n
            ) &
            EB_SER_max_overlap_n > 0,
            "partial_cluster_agreement",
            "discordant_cluster"
          )
        )
      )
  ]

  # Separate partner-overlap description.
  SUMMARY[
    ,
    partner_cluster_label :=
      fifelse(
        !partner_fit_available,
        "partner_fit_unavailable",
        fifelse(
          is.finite(
            SER_partner_max_overlap_coefficient
          ) &
          SER_partner_max_overlap_coefficient >= 0.80,
          "strong_partner_cluster_overlap",
          fifelse(
            is.finite(
              SER_partner_max_overlap_coefficient
            ) &
            SER_partner_max_overlap_coefficient > 0,
            "partial_partner_cluster_overlap",
            "no_partner_cluster_overlap"
          )
        )
      )
  ]
}

# ============================================================
# 8. SAVE OUTPUTS
# ============================================================

SUMMARY_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_SER_audit_summary.csv"
)

OVERLAP_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_EB_vs_SER_CS_overlap.csv"
)

SER_CS_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_SER_credible_sets.csv"
)

PARTNER_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_SER_vs_partner_EB_CS_overlap.csv"
)

VARIANT_FILE <- file.path(
  AUDIT_VARIANT_DIR,
  "STEP8C_V4_EB_vs_SER_variant_PIPs.csv.gz"
)

QC_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_processing_QC.csv"
)

fwrite(
  SUMMARY,
  SUMMARY_FILE
)

fwrite(
  OVERLAP,
  OVERLAP_FILE
)

fwrite(
  SER_CS,
  SER_CS_FILE
)

fwrite(
  PARTNER,
  PARTNER_FILE
)

fwrite(
  VARIANT,
  VARIANT_FILE
)

fwrite(
  QC,
  QC_FILE
)

# ============================================================
# 9. METHODS / PROVENANCE
# ============================================================

METHOD <- data.table(
  item = c(
    "analysis",
    "input_source",
    "susieR_version",
    "fits_refit",
    "region_changed",
    "SNP_filter_changed",
    "prior_changed",
    "SER_source",
    "SER_CS_priority",
    "descriptive_overlap_threshold_note"
  ),
  value = c(
    "STEP8C-V4.1 SER reliability audit",
    "Saved mismatch-aware STEP8C-V3 SuSiE fits",
    SUSIER_VERSION,
    "No",
    "No",
    "No",
    "No",
    "fit$R_finite_diagnostics$ser_model",
    "Use ser_model$sets$cs when available; otherwise construct cumulative 95% PIP set",
    "0.80 overlap coefficient is used only as a descriptive audit label, not an official susieR threshold"
  )
)

METHOD_FILE <- file.path(
  AUDIT_TABLE_DIR,
  "STEP8C_V4_method_provenance.csv"
)

fwrite(
  METHOD,
  METHOD_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    AUDIT_LOG_DIR,
    "STEP8C_V4_sessionInfo.txt"
  )
)

# ============================================================
# 10. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP8C-V4.1 COMPLETE — SER RELIABILITY AUDIT\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nProcessing QC:\n"
)

print(
  QC[
    ,
    .N,
    by = status
  ]
)

if (nrow(SUMMARY)) {

  cat(
    "\nFlagged-fit audit summary:\n"
  )

  print(
    SUMMARY[
      ,
      .(
        test_id,
        trait,
        partner_trait,
        EB_n_CS,
        EB_top_snp,
        EB_top_PIP,
        SER_top_snp,
        SER_top_PIP,
        SER_CS_size,
        top_SNP_identical,
        EB_SER_max_overlap_coefficient,
        PIP_Spearman,
        SER_partner_max_overlap_coefficient,
        audit_label,
        partner_cluster_label
      )
    ]
  )
}

cat(
  "\nUPLOAD THESE 5 FILES:\n",
  "1) ", SUMMARY_FILE, "\n",
  "2) ", OVERLAP_FILE, "\n",
  "3) ", SER_CS_FILE, "\n",
  "4) ", PARTNER_FILE, "\n",
  "5) ", QC_FILE, "\n",
  sep = ""
)

cat(
  "\nOptional only if we need SNP-level inspection:\n",
  VARIANT_FILE,
  "\n",
  sep = ""
)

cat(
  "\nInterpretation rule:\n",
  "- V4 does NOT create new colocalization evidence.\n",
  "- It audits whether reliability-flagged V3 fine-mapping ",
  "signals survive under the official one-effect SER diagnostic.\n",
  "- If EB and SER disagree materially, use the SER signal as ",
  "the safer fallback and downgrade fine-mapping specificity.\n",
  "====================================================\n",
  sep = ""
)
