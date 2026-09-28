# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7D_RUN_BIVARIATE_LOCAL_RG(1).R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7D — BIVARIATE LOCAL rg FOR THE 6 TRAIT PAIRS
# LAVA v0.1.5
#
# Uses:
#   - STEP7B2_V4_LAVA_input.RData
#   - STEP7C_V2_pairwise_bivariate_eligible_loci.csv
#
# Primary design:
#   - Only loci where BOTH phenotypes passed the STEP7C
#     Bonferroni local-h2 gate are tested.
#   - Pair-specific Bonferroni threshold = 0.05 / number of
#     eligible loci for that pair.
#   - BH-FDR across ALL valid bivariate tests is also reported
#     as a secondary/sensitivity correction.
#
# Outputs:
#   04_bivariate/STEP7D_local_rg_all.csv
#   04_bivariate/STEP7D_local_rg_Bonferroni_significant.csv
#   04_bivariate/STEP7D_local_rg_BH_FDR_significant.csv
#   07_tables/STEP7D_pair_summary.csv
#   07_tables/STEP7D_processing_QC.csv
#   07_tables/STEP7D_direction_summary.csv
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")

INPUT_RDATA <- file.path(
  LAVA_ROOT,
  "02_input",
  "STEP7B2_V4_LAVA_input.RData"
)

ELIG_FILE <- file.path(
  LAVA_ROOT,
  "07_tables",
  "STEP7C_V2_pairwise_bivariate_eligible_loci.csv"
)

PAIR_SUMMARY_FILE <- file.path(
  LAVA_ROOT,
  "07_tables",
  "STEP7C_V2_pairwise_eligibility_summary.csv"
)

OUT_DIR <- file.path(
  LAVA_ROOT,
  "04_bivariate"
)

TAB_DIR <- file.path(
  LAVA_ROOT,
  "07_tables"
)

LOG_DIR <- file.path(
  LAVA_ROOT,
  "08_logs"
)

SOFTWARE_DIR <- file.path(
  LAVA_ROOT,
  "00_software"
)

for (d in c(OUT_DIR, TAB_DIR, LOG_DIR, SOFTWARE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PHEN <- c("AF", "HFpEF", "BMI", "OSA")

# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)

if (!requireNamespace("LAVA", quietly = TRUE)) {
  stop("LAVA is not installed.")
}

library(LAVA)

lava_version <- as.character(
  packageVersion("LAVA")
)

if (lava_version != "0.1.5") {
  warning(
    "Expected LAVA v0.1.5, detected v",
    lava_version
  )
}

# ============================================================
# 2. LOAD INPUTS
# ============================================================

if (!file.exists(INPUT_RDATA)) {
  stop(
    "Missing processed LAVA input:\n",
    INPUT_RDATA
  )
}

if (!file.exists(ELIG_FILE)) {
  stop(
    "Missing STEP7C eligible loci file:\n",
    ELIG_FILE
  )
}

if (!file.exists(PAIR_SUMMARY_FILE)) {
  stop(
    "Missing STEP7C pair summary:\n",
    PAIR_SUMMARY_FILE
  )
}

load(INPUT_RDATA)

if (!exists("input") || !is.environment(input)) {
  stop(
    "Invalid LAVA input object in:\n",
    INPUT_RDATA
  )
}

ELIG <- fread(ELIG_FILE)
PAIR_SUMMARY <- fread(PAIR_SUMMARY_FILE)

required_elig <- c(
  "pair_key",
  "locus",
  "chr",
  "start",
  "stop",
  "trait1",
  "trait2"
)

if (!all(required_elig %in% names(ELIG))) {
  stop(
    "Eligible-loci file lacks columns: ",
    paste(
      setdiff(required_elig, names(ELIG)),
      collapse = ", "
    )
  )
}

if (!nrow(ELIG)) {
  stop("No loci are eligible for bivariate local rg.")
}

ELIG[, locus := as.character(locus)]

# ============================================================
# 3. BASIC PAIR QC
# ============================================================

expected_pairs <- c(
  "AF__HFpEF",
  "AF__BMI",
  "AF__OSA",
  "HFpEF__BMI",
  "HFpEF__OSA",
  "BMI__OSA"
)

if (!setequal(unique(ELIG$pair_key), expected_pairs)) {
  warning(
    "Eligible-loci file does not contain exactly all expected pair keys."
  )
}

# Pair-specific primary thresholds
PAIR_THRESH <- PAIR_SUMMARY[
  ,
  .(
    trait1,
    trait2,
    n_bivariate_eligible_loci,
    bivar_bonferroni_threshold
  )
]

PAIR_THRESH[
  ,
  pair_key := paste(
    trait1,
    trait2,
    sep = "__"
  )
]

# ============================================================
# 4. LOCUS DEFINITIONS DIRECTLY FROM STEP7C
#
# The eligible file already stores the exact standard LAVA
# locus coordinates, so no separate locus file is needed.
# ============================================================

make_loc_def <- function(row) {
  data.frame(
    LOC = as.character(row$locus),
    CHR = as.integer(row$chr),
    START = as.numeric(row$start),
    STOP = as.numeric(row$stop),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 5. SAFE WRAPPERS
# ============================================================

safe_process_locus <- function(loc_def, pair_phenos) {

  printed <- character(0)
  warnings_vec <- character(0)
  err <- NULL
  obj <- NULL

  obj <- tryCatch(
    withCallingHandlers(
      {
        printed <- capture.output(
          obj_inner <- process.locus(
            loc_def,
            input,
            phenos = pair_phenos
          ),
          type = "output"
        )
        obj_inner
      },
      warning = function(w) {
        warnings_vec <<- c(
          warnings_vec,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(
    value = obj,
    error = err,
    printed = unique(printed),
    warnings = unique(warnings_vec)
  )
}

safe_run_bivar <- function(locus, pair_phenos) {

  messages_vec <- character(0)
  warnings_vec <- character(0)
  err <- NULL
  obj <- NULL

  obj <- tryCatch(
    withCallingHandlers(
      {
        obj_inner <- withCallingHandlers(
          run.bivar(
            locus,
            phenos = pair_phenos,
            p.values = TRUE,
            CIs = TRUE,
            param.lim = 1.25,
            cap.estimates = TRUE
          ),
          message = function(m) {
            messages_vec <<- c(
              messages_vec,
              conditionMessage(m)
            )
            invokeRestart("muffleMessage")
          }
        )
        obj_inner
      },
      warning = function(w) {
        warnings_vec <<- c(
          warnings_vec,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(
    value = obj,
    error = err,
    messages = unique(messages_vec),
    warnings = unique(warnings_vec)
  )
}

compact_msg <- function(x, max_n = 6L) {

  x <- trimws(
    as.character(x)
  )

  x <- unique(
    x[nzchar(x)]
  )

  if (!length(x)) {
    return("")
  }

  if (length(x) > max_n) {
    x <- c(
      x[seq_len(max_n)],
      paste0(
        "... [",
        length(x) - max_n,
        " more]"
      )
    )
  }

  paste(
    x,
    collapse = " || "
  )
}

# ============================================================
# 6. RUN ALL ELIGIBLE PAIR-LOCUS TESTS
# ============================================================

RESULT_LIST <- list()
QC_LIST <- list()

analysis_start <- Sys.time()

n_tests <- nrow(ELIG)

cat(
  "\nTotal eligible pair-locus tests: ",
  n_tests,
  "\n",
  sep = ""
)

for (i in seq_len(n_tests)) {

  if (
    i == 1 ||
    i %% 25 == 0 ||
    i == n_tests
  ) {
    cat(
      "STEP7D progress: ",
      i,
      "/",
      n_tests,
      " (",
      round(
        100 * i / n_tests,
        1
      ),
      "%)\n",
      sep = ""
    )
  }

  row <- ELIG[i]

  pair_phenos <- c(
    as.character(row$trait1),
    as.character(row$trait2)
  )

  loc_def <- make_loc_def(row)

  proc <- safe_process_locus(
    loc_def,
    pair_phenos
  )

  proc_text <- compact_msg(
    c(
      proc$printed,
      proc$warnings,
      proc$error
    )
  )

  locus <- proc$value

  if (
    is.null(locus) ||
    !all(pair_phenos %in% locus$phenos)
  ) {

    QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
      pair_key = row$pair_key,
      trait1 = row$trait1,
      trait2 = row$trait2,
      locus = row$locus,
      chr = row$chr,
      start = row$start,
      stop = row$stop,
      status = "PROCESS_FAIL",
      message = proc_text
    )

    next
  }

  rb <- safe_run_bivar(
    locus,
    pair_phenos
  )

  b <- rb$value

  b_text <- compact_msg(
    c(
      rb$messages,
      rb$warnings,
      rb$error
    )
  )

  if (
    is.null(b) ||
    !is.data.frame(b) ||
    nrow(b) != 1
  ) {

    QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
      pair_key = row$pair_key,
      trait1 = row$trait1,
      trait2 = row$trait2,
      locus = row$locus,
      chr = row$chr,
      start = row$start,
      stop = row$stop,
      status = "BIVAR_FAIL",
      message = compact_msg(
        c(
          proc_text,
          b_text
        )
      )
    )

    next
  }

  required_biv <- c(
    "phen1",
    "phen2",
    "rho",
    "rho.lower",
    "rho.upper",
    "r2",
    "r2.lower",
    "r2.upper",
    "p"
  )

  if (!all(required_biv %in% names(b))) {

    QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
      pair_key = row$pair_key,
      trait1 = row$trait1,
      trait2 = row$trait2,
      locus = row$locus,
      chr = row$chr,
      start = row$start,
      stop = row$stop,
      status = "FORMAT_FAIL",
      message = paste0(
        "run.bivar columns: ",
        paste(
          names(b),
          collapse = ","
        )
      )
    )

    next
  }

  z <- as.data.table(b)

  z[
    ,
    `:=`(
      pair_key = as.character(row$pair_key),
      locus = as.character(row$locus),
      chr = as.integer(row$chr),
      start = as.numeric(row$start),
      stop = as.numeric(row$stop),
      n_snps = as.integer(locus$n.snps),
      n_pcs = as.integer(locus$K),
      local_h2_p_trait1 = as.numeric(row$p_trait1),
      local_h2_p_trait2 = as.numeric(row$p_trait2),
      local_h2_trait1 = as.numeric(row$h2_trait1),
      local_h2_trait2 = as.numeric(row$h2_trait2)
    )
  ]

  RESULT_LIST[[length(RESULT_LIST) + 1L]] <- z

  QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
    pair_key = row$pair_key,
    trait1 = row$trait1,
    trait2 = row$trait2,
    locus = row$locus,
    chr = row$chr,
    start = row$start,
    stop = row$stop,
    status = "PASS",
    message = compact_msg(
      c(
        proc_text,
        b_text
      )
    )
  )
}

analysis_end <- Sys.time()

# ============================================================
# 7. COMBINE RESULTS
# ============================================================

RESULTS <- if (length(RESULT_LIST)) {
  rbindlist(
    RESULT_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

QC <- if (length(QC_LIST)) {
  rbindlist(
    QC_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (!nrow(RESULTS)) {
  stop(
    "STEP7D produced no valid bivariate results."
  )
}

# ============================================================
# 8. ADD MULTIPLE-TESTING CORRECTION
# ============================================================

RESULTS <- merge(
  RESULTS,
  PAIR_THRESH[
    ,
    .(
      pair_key,
      n_bivariate_eligible_loci,
      bivar_bonferroni_threshold
    )
  ],
  by = "pair_key",
  all.x = TRUE
)

RESULTS[
  ,
  pass_pair_bonferroni :=
    is.finite(p) &
    p < bivar_bonferroni_threshold
]

# BH-FDR across all valid bivariate tests, secondary only
RESULTS[
  ,
  p_BH_all :=
    p.adjust(
      p,
      method = "BH"
    )
]

RESULTS[
  ,
  pass_BH_all :=
    is.finite(p_BH_all) &
    p_BH_all < 0.05
]

# Also BH within pair for descriptive sensitivity
RESULTS[
  ,
  p_BH_within_pair :=
    p.adjust(
      p,
      method = "BH"
    ),
  by = pair_key
]

RESULTS[
  ,
  pass_BH_within_pair :=
    is.finite(p_BH_within_pair) &
    p_BH_within_pair < 0.05
]

# Local rg direction
RESULTS[
  ,
  direction :=
    fifelse(
      is.na(rho),
      "NA",
      fifelse(
        rho > 0,
        "positive",
        fifelse(
          rho < 0,
          "negative",
          "zero"
        )
      )
    )
]

# ============================================================
# 9. SAVE FULL AND SIGNIFICANT RESULTS
# ============================================================

ALL_FILE <- file.path(
  OUT_DIR,
  "STEP7D_local_rg_all.csv"
)

BONF_FILE <- file.path(
  OUT_DIR,
  "STEP7D_local_rg_Bonferroni_significant.csv"
)

BH_FILE <- file.path(
  OUT_DIR,
  "STEP7D_local_rg_BH_FDR_significant.csv"
)

fwrite(
  RESULTS,
  ALL_FILE
)

fwrite(
  RESULTS[
    pass_pair_bonferroni == TRUE
  ],
  BONF_FILE
)

fwrite(
  RESULTS[
    pass_BH_all == TRUE
  ],
  BH_FILE
)

# ============================================================
# 10. PAIR SUMMARY
# ============================================================

PAIR_OUT <- RESULTS[
  ,
  .(
    eligible_loci = unique(
      n_bivariate_eligible_loci
    )[1],
    valid_rg_tests = .N,
    failed_or_filtered = unique(
      n_bivariate_eligible_loci
    )[1] - .N,
    bonferroni_threshold = unique(
      bivar_bonferroni_threshold
    )[1],
    n_bonferroni_sig = sum(
      pass_pair_bonferroni,
      na.rm = TRUE
    ),
    n_BH_all_sig = sum(
      pass_BH_all,
      na.rm = TRUE
    ),
    n_BH_within_pair_sig = sum(
      pass_BH_within_pair,
      na.rm = TRUE
    ),
    n_positive = sum(
      direction == "positive",
      na.rm = TRUE
    ),
    n_negative = sum(
      direction == "negative",
      na.rm = TRUE
    ),
    median_rho = median(
      rho,
      na.rm = TRUE
    ),
    min_p = min(
      p,
      na.rm = TRUE
    )
  ),
  by = pair_key
]

PAIR_OUT[
  ,
  c("trait1", "trait2") :=
    tstrsplit(
      pair_key,
      "__",
      fixed = TRUE
    )
]

setcolorder(
  PAIR_OUT,
  c(
    "pair_key",
    "trait1",
    "trait2",
    setdiff(
      names(PAIR_OUT),
      c(
        "pair_key",
        "trait1",
        "trait2"
      )
    )
  )
)

PAIR_OUT_FILE <- file.path(
  TAB_DIR,
  "STEP7D_pair_summary.csv"
)

fwrite(
  PAIR_OUT,
  PAIR_OUT_FILE
)

# ============================================================
# 11. DIRECTION SUMMARY FOR PRIMARY-SIGNIFICANT RESULTS
# ============================================================

DIRECTION_SUMMARY <- RESULTS[
  pass_pair_bonferroni == TRUE,
  .(
    n = .N,
    median_rho = median(
      rho,
      na.rm = TRUE
    ),
    min_rho = min(
      rho,
      na.rm = TRUE
    ),
    max_rho = max(
      rho,
      na.rm = TRUE
    )
  ),
  by = .(
    pair_key,
    direction
  )
]

DIRECTION_FILE <- file.path(
  TAB_DIR,
  "STEP7D_direction_summary.csv"
)

fwrite(
  DIRECTION_SUMMARY,
  DIRECTION_FILE
)

# ============================================================
# 12. PROCESSING QC
# ============================================================

QC_SUMMARY <- QC[
  ,
  .N,
  by = .(
    pair_key,
    status
  )
]

QC_SUMMARY[
  ,
  runtime_minutes :=
    as.numeric(
      difftime(
        analysis_end,
        analysis_start,
        units = "mins"
      )
    )
]

QC_FILE <- file.path(
  TAB_DIR,
  "STEP7D_processing_QC.csv"
)

fwrite(
  QC_SUMMARY,
  QC_FILE
)

QC_DETAIL_FILE <- file.path(
  OUT_DIR,
  "STEP7D_processing_QC_detail.csv"
)

fwrite(
  QC,
  QC_DETAIL_FILE
)

# ============================================================
# 13. PRIMARY TARGET-REGION STATUS
#
# These are the previously frozen recurrent/important FUMA
# regions. This table does NOT force a bivariate test.
# It only records whether they qualified under the primary
# STEP7C local-h2 screen.
# ============================================================

TARGET_REGIONS <- data.table(
  region = c(
    "chr11_recurrent",
    "chr3_BMI_centered",
    "chr2_AF_HFpEF"
  ),
  chr = c(
    11L,
    3L,
    2L
  ),
  start = c(
    16410000,
    185490000,
    200700000
  ),
  stop = c(
    16790000,
    185850000,
    201120000
  )
)

TARGET_STATUS_LIST <- list()

for (r in seq_len(nrow(TARGET_REGIONS))) {

  rr <- TARGET_REGIONS[r]

  x <- ELIG[
    chr == rr$chr &
    stop >= rr$start &
    start <= rr$stop
  ]

  TARGET_STATUS_LIST[[r]] <- data.table(
    region = rr$region,
    chr = rr$chr,
    start = rr$start,
    stop = rr$stop,
    n_primary_eligible_pair_tests = nrow(x),
    eligible_pairs = if (nrow(x)) {
      paste(
        unique(x$pair_key),
        collapse = ";"
      )
    } else {
      ""
    },
    primary_local_rg_status = if (nrow(x)) {
      "Eligible in STEP7D"
    } else {
      "Not eligible: both traits did not pass strict local-h2 gate in same LAVA block"
    }
  )
}

TARGET_STATUS <- rbindlist(
  TARGET_STATUS_LIST
)

TARGET_STATUS_FILE <- file.path(
  TAB_DIR,
  "STEP7D_key_FUMA_region_primary_LAVA_status.csv"
)

fwrite(
  TARGET_STATUS,
  TARGET_STATUS_FILE
)

# ============================================================
# 14. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP7D_sessionInfo.txt"
  )
)

# ============================================================
# 15. FINAL REPORT
# ============================================================

cat("\n====================================================\n")
cat("STEP7D COMPLETE — BIVARIATE LOCAL rg\n")
cat("====================================================\n")

cat(
  "\nValid bivariate tests: ",
  nrow(RESULTS),
  " / ",
  n_tests,
  "\n",
  sep = ""
)

cat("\nPair summary:\n")
print(
  PAIR_OUT
)

cat("\nPrimary-significant direction summary:\n")
print(
  DIRECTION_SUMMARY
)

cat("\nKey recurrent FUMA region status:\n")
print(
  TARGET_STATUS
)

cat("\nProcessing QC:\n")
print(
  QC_SUMMARY
)

cat("\nUPLOAD THESE 6 FILES:\n")
cat("1) ", PAIR_OUT_FILE, "\n", sep = "")
cat("2) ", BONF_FILE, "\n", sep = "")
cat("3) ", ALL_FILE, "\n", sep = "")
cat("4) ", QC_FILE, "\n", sep = "")
cat("5) ", DIRECTION_FILE, "\n", sep = "")
cat("6) ", TARGET_STATUS_FILE, "\n", sep = "")

cat("\nDo NOT start multivariate LAVA until these results are checked.\n")
cat("====================================================\n")
