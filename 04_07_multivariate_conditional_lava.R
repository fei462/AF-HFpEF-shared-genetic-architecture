# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7E_RUN_MULTIVARIATE_CONDITIONAL_LAVA.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7E — MULTIVARIATE / CONDITIONAL LAVA
# Automatically identify triads in which ALL THREE pairwise
# local genetic correlations passed STEP7D pair-specific
# Bonferroni correction, then run:
#
#   1) Partial correlations:
#      A-B | C
#      A-C | B
#      B-C | A
#
#   2) Multiple regression:
#      A ~ B + C
#      B ~ A + C
#      C ~ A + B
#
# Primary scientific focus:
#   AF / BMI / OSA triad if selected by STEP7D.
#
# LAVA v0.1.5
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

STEP7D_ALL <- file.path(
  LAVA_ROOT,
  "04_bivariate",
  "STEP7D_local_rg_all.csv"
)

OUT_DIR <- file.path(
  LAVA_ROOT,
  "05_multivariate"
)

TAB_DIR <- file.path(
  LAVA_ROOT,
  "07_tables"
)

LOG_DIR <- file.path(
  LAVA_ROOT,
  "08_logs"
)

for (d in c(OUT_DIR, TAB_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

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

lava_version <- as.character(packageVersion("LAVA"))

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
  stop("Missing: ", INPUT_RDATA)
}

if (!file.exists(STEP7D_ALL)) {
  stop("Missing: ", STEP7D_ALL)
}

load(INPUT_RDATA)

if (!exists("input") || !is.environment(input)) {
  stop("Invalid LAVA input object.")
}

BIV <- fread(STEP7D_ALL)

required_biv <- c(
  "pair_key",
  "locus",
  "chr",
  "start",
  "stop",
  "phen1",
  "phen2",
  "rho",
  "p",
  "pass_pair_bonferroni"
)

if (!all(required_biv %in% names(BIV))) {
  stop(
    "STEP7D file lacks: ",
    paste(
      setdiff(required_biv, names(BIV)),
      collapse = ", "
    )
  )
}

BIV[, locus := as.character(locus)]

# ============================================================
# 3. IDENTIFY COMPLETE SIGNIFICANT TRIADS
#
# We do not manually cherry-pick locus 2135.
# A triad is eligible only when all 3 unique pairwise edges
# at the same LAVA locus passed pair-specific Bonferroni.
# ============================================================

ALL_PHEN <- c("AF", "HFpEF", "BMI", "OSA")

TRIADS <- combn(
  ALL_PHEN,
  3,
  simplify = FALSE
)

canonical_pair <- function(a, b) {
  ord <- match(c(a, b), ALL_PHEN)
  paste(
    c(a, b)[order(ord)],
    collapse = "__"
  )
}

candidate_list <- list()

for (tri in TRIADS) {

  required_pairs <- c(
    canonical_pair(tri[1], tri[2]),
    canonical_pair(tri[1], tri[3]),
    canonical_pair(tri[2], tri[3])
  )

  sig_edges <- BIV[
    pass_pair_bonferroni == TRUE &
    pair_key %in% required_pairs
  ]

  if (!nrow(sig_edges)) next

  by_locus <- sig_edges[
    ,
    .(
      n_required_edges = uniqueN(pair_key),
      pair_keys = paste(
        sort(unique(pair_key)),
        collapse = ";"
      )
    ),
    by = .(
      locus,
      chr,
      start,
      stop
    )
  ]

  ok <- by_locus[
    n_required_edges == 3
  ]

  if (nrow(ok)) {
    ok[
      ,
      `:=`(
        triad = paste(tri, collapse = "__"),
        phen1 = tri[1],
        phen2 = tri[2],
        phen3 = tri[3]
      )
    ]
    candidate_list[[length(candidate_list) + 1L]] <- ok
  }
}

CAND <- if (length(candidate_list)) {
  rbindlist(candidate_list, fill = TRUE)
} else {
  data.table()
}

if (!nrow(CAND)) {
  stop(
    "No locus has all three pairwise local rg values ",
    "significant after STEP7D Bonferroni correction."
  )
}

CAND_FILE <- file.path(
  TAB_DIR,
  "STEP7E_multivariate_candidate_triads.csv"
)

fwrite(CAND, CAND_FILE)

cat("\n================ STEP7E CANDIDATE TRIADS ================\n")
print(CAND)

# ============================================================
# 4. SAFE HELPERS
# ============================================================

safe_call <- function(expr) {

  warnings_vec <- character(0)
  messages_vec <- character(0)
  err <- NULL
  value <- NULL

  value <- tryCatch(
    withCallingHandlers(
      withCallingHandlers(
        eval.parent(substitute(expr)),
        message = function(m) {
          messages_vec <<- c(
            messages_vec,
            conditionMessage(m)
          )
          invokeRestart("muffleMessage")
        }
      ),
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
    ok = is.null(err),
    value = value,
    error = err,
    warnings = unique(warnings_vec),
    messages = unique(messages_vec)
  )
}

compact_msg <- function(x, max_n = 8L) {

  x <- trimws(as.character(x))
  x <- unique(x[nzchar(x)])

  if (!length(x)) return("")

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

  paste(x, collapse = " || ")
}

flatten_multireg <- function(x) {

  if (is.null(x)) return(data.table())

  rows <- list()

  walk <- function(obj) {
    if (is.data.frame(obj)) {
      rows[[length(rows) + 1L]] <<- as.data.table(obj)
    } else if (is.list(obj)) {
      for (j in seq_along(obj)) {
        walk(obj[[j]])
      }
    }
  }

  walk(x)

  if (!length(rows)) {
    return(data.table())
  }

  rbindlist(rows, fill = TRUE)
}

# ============================================================
# 5. RUN CONDITIONAL ANALYSES
# ============================================================

PCOR_LIST <- list()
MULTIREG_LIST <- list()
QC_LIST <- list()
BIV_RECHECK_LIST <- list()
UNIV_RECHECK_LIST <- list()

for (i in seq_len(nrow(CAND))) {

  row <- CAND[i]

  phenos3 <- c(
    row$phen1,
    row$phen2,
    row$phen3
  )

  loc_def <- data.frame(
    LOC = as.character(row$locus),
    CHR = as.integer(row$chr),
    START = as.numeric(row$start),
    STOP = as.numeric(row$stop),
    stringsAsFactors = FALSE
  )

  cat(
    "\nProcessing multivariate locus ",
    row$locus,
    " : ",
    paste(phenos3, collapse = " / "),
    "\n",
    sep = ""
  )

  proc <- safe_call(
    process.locus(
      loc_def,
      input,
      phenos = phenos3
    )
  )

  if (
    !proc$ok ||
    is.null(proc$value) ||
    !all(phenos3 %in% proc$value$phenos)
  ) {

    QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
      triad = row$triad,
      locus = row$locus,
      status = "PROCESS_FAIL",
      message = compact_msg(
        c(
          proc$error,
          proc$warnings,
          proc$messages
        )
      )
    )

    next
  }

  locus <- proc$value

  # ----------------------------------------------------------
  # Re-check univariate and bivariate results at same locus
  # ----------------------------------------------------------

  uu <- safe_call(
    run.univ(
      locus,
      phenos = phenos3
    )
  )

  if (uu$ok && is.data.frame(uu$value)) {
    u <- as.data.table(uu$value)
    u[
      ,
      `:=`(
        triad = as.character(row$triad),
        locus = as.character(row$locus),
        chr = as.integer(row$chr),
        start = as.numeric(row$start),
        stop = as.numeric(row$stop)
      )
    ]
    UNIV_RECHECK_LIST[[length(UNIV_RECHECK_LIST) + 1L]] <- u
  }

  bb <- safe_call(
    run.bivar(
      locus,
      phenos = phenos3,
      p.values = TRUE,
      CIs = TRUE,
      param.lim = 1.25,
      cap.estimates = TRUE
    )
  )

  if (bb$ok && is.data.frame(bb$value)) {
    b <- as.data.table(bb$value)
    b[
      ,
      `:=`(
        triad = as.character(row$triad),
        locus = as.character(row$locus),
        chr = as.integer(row$chr),
        start = as.numeric(row$start),
        stop = as.numeric(row$stop)
      )
    ]
    BIV_RECHECK_LIST[[length(BIV_RECHECK_LIST) + 1L]] <- b
  }

  # ----------------------------------------------------------
  # Partial correlations: each pair conditioned on third
  # ----------------------------------------------------------

  targets <- list(
    c(phenos3[1], phenos3[2], phenos3[3]),
    c(phenos3[1], phenos3[3], phenos3[2]),
    c(phenos3[2], phenos3[3], phenos3[1])
  )

  pcor_ok <- 0L

  for (tt in targets) {

    target_pair <- tt[1:2]
    conditioner <- tt[3]

    rr <- safe_call(
      run.pcor(
        locus,
        target = target_pair,
        phenos = conditioner,
        p.values = TRUE,
        CIs = TRUE,
        max.r2 = 0.95,
        param.lim = 1.25
      )
    )

    if (
      rr$ok &&
      is.data.frame(rr$value) &&
      nrow(rr$value) == 1
    ) {

      pcor_ok <- pcor_ok + 1L

      z <- as.data.table(rr$value)

      z[
        ,
        `:=`(
          triad = as.character(row$triad),
          locus = as.character(row$locus),
          chr = as.integer(row$chr),
          start = as.numeric(row$start),
          stop = as.numeric(row$stop),
          analysis_label = paste0(
            target_pair[1],
            "__",
            target_pair[2],
            "_given_",
            conditioner
          )
        )
      ]

      PCOR_LIST[[length(PCOR_LIST) + 1L]] <- z

    } else {

      QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
        triad = row$triad,
        locus = row$locus,
        status = "PCOR_FAIL",
        message = paste0(
          paste(target_pair, collapse = "-"),
          "|",
          conditioner,
          ": ",
          compact_msg(
            c(
              rr$error,
              rr$warnings,
              rr$messages
            )
          )
        )
      )
    }
  }

  # ----------------------------------------------------------
  # Multiple regression: each phenotype as outcome
  # ----------------------------------------------------------

  multireg_ok <- 0L

  for (target in phenos3) {

    predictors <- setdiff(
      phenos3,
      target
    )

    rr <- safe_call(
      run.multireg(
        locus,
        target = target,
        phenos = predictors,
        only.full.model = TRUE,
        p.values = TRUE,
        CIs = TRUE,
        param.lim = 1.5,
        suppress.message = TRUE
      )
    )

    m <- flatten_multireg(
      rr$value
    )

    if (
      rr$ok &&
      nrow(m) > 0
    ) {

      multireg_ok <- multireg_ok + 1L

      m[
        ,
        `:=`(
          triad = as.character(row$triad),
          locus = as.character(row$locus),
          chr = as.integer(row$chr),
          start = as.numeric(row$start),
          stop = as.numeric(row$stop),
          model_label = paste0(
            target,
            "~",
            paste(
              predictors,
              collapse = "+"
            )
          )
        )
      ]

      MULTIREG_LIST[[length(MULTIREG_LIST) + 1L]] <- m

    } else {

      QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
        triad = row$triad,
        locus = row$locus,
        status = "MULTIREG_FAIL",
        message = paste0(
          target,
          "~",
          paste(predictors, collapse = "+"),
          ": ",
          compact_msg(
            c(
              rr$error,
              rr$warnings,
              rr$messages
            )
          )
        )
      )
    }
  }

  QC_LIST[[length(QC_LIST) + 1L]] <- data.table(
    triad = row$triad,
    locus = row$locus,
    status = "PASS",
    message = paste0(
      "pcor_success=",
      pcor_ok,
      "/3; multireg_success=",
      multireg_ok,
      "/3"
    )
  )
}

# ============================================================
# 6. COMBINE
# ============================================================

PCOR <- if (length(PCOR_LIST)) {
  rbindlist(
    PCOR_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

MULTIREG <- if (length(MULTIREG_LIST)) {
  rbindlist(
    MULTIREG_LIST,
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

UNIV_RECHECK <- if (length(UNIV_RECHECK_LIST)) {
  rbindlist(
    UNIV_RECHECK_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

BIV_RECHECK <- if (length(BIV_RECHECK_LIST)) {
  rbindlist(
    BIV_RECHECK_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (!nrow(PCOR)) {
  stop("No valid partial-correlation results produced.")
}

if (!nrow(MULTIREG)) {
  stop("No valid multiple-regression results produced.")
}

# ============================================================
# 7. MULTIPLE-TESTING CORRECTION
# ============================================================

# Partial correlations: all conditional pair tests generated here
PCOR[
  ,
  p_Bonferroni :=
    p.adjust(
      p,
      method = "bonferroni"
    )
]

PCOR[
  ,
  p_BH :=
    p.adjust(
      p,
      method = "BH"
    )
]

PCOR[
  ,
  pass_Bonferroni :=
    is.finite(p_Bonferroni) &
    p_Bonferroni < 0.05
]

PCOR[
  ,
  pass_BH :=
    is.finite(p_BH) &
    p_BH < 0.05
]

# Multiple regression: coefficient-level tests
MULTIREG[
  ,
  p_Bonferroni :=
    p.adjust(
      p,
      method = "bonferroni"
    )
]

MULTIREG[
  ,
  p_BH :=
    p.adjust(
      p,
      method = "BH"
    )
]

MULTIREG[
  ,
  pass_Bonferroni :=
    is.finite(p_Bonferroni) &
    p_Bonferroni < 0.05
]

MULTIREG[
  ,
  pass_BH :=
    is.finite(p_BH) &
    p_BH < 0.05
]

# ============================================================
# 8. SAVE RESULTS
# ============================================================

PCOR_FILE <- file.path(
  OUT_DIR,
  "STEP7E_partial_correlations.csv"
)

MULTIREG_FILE <- file.path(
  OUT_DIR,
  "STEP7E_multiple_regression.csv"
)

UNIV_FILE <- file.path(
  OUT_DIR,
  "STEP7E_univ_recheck.csv"
)

BIV_FILE <- file.path(
  OUT_DIR,
  "STEP7E_bivar_recheck.csv"
)

QC_FILE <- file.path(
  TAB_DIR,
  "STEP7E_processing_QC.csv"
)

fwrite(PCOR, PCOR_FILE)
fwrite(MULTIREG, MULTIREG_FILE)
fwrite(UNIV_RECHECK, UNIV_FILE)
fwrite(BIV_RECHECK, BIV_FILE)
fwrite(QC, QC_FILE)

# ============================================================
# 9. PRIMARY AF-FOCUSED SUMMARY
# ============================================================

AF_PCOR <- PCOR[
  phen1 == "AF" |
  phen2 == "AF"
]

AF_MREG <- MULTIREG[
  outcome == "AF"
]

AF_SUMMARY_FILE <- file.path(
  TAB_DIR,
  "STEP7E_AF_focused_summary.csv"
)

# Store both classes in long form
AF_SUMMARY <- rbindlist(
  list(
    if (nrow(AF_PCOR)) {
      data.table(
        analysis_type = "partial_correlation",
        label = AF_PCOR$analysis_label,
        estimate = AF_PCOR$pcor,
        ci_lower = AF_PCOR$ci.lower,
        ci_upper = AF_PCOR$ci.upper,
        p = AF_PCOR$p,
        p_Bonferroni = AF_PCOR$p_Bonferroni,
        pass_Bonferroni = AF_PCOR$pass_Bonferroni
      )
    } else {
      data.table()
    },
    if (nrow(AF_MREG)) {
      data.table(
        analysis_type = "multiple_regression",
        label = paste0(
          AF_MREG$outcome,
          "~",
          AF_MREG$predictors
        ),
        estimate = AF_MREG$gamma,
        ci_lower = AF_MREG$gamma.lower,
        ci_upper = AF_MREG$gamma.upper,
        p = AF_MREG$p,
        p_Bonferroni = AF_MREG$p_Bonferroni,
        pass_Bonferroni = AF_MREG$pass_Bonferroni
      )
    } else {
      data.table()
    }
  ),
  fill = TRUE
)

fwrite(
  AF_SUMMARY,
  AF_SUMMARY_FILE
)

# ============================================================
# 10. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP7E_sessionInfo.txt"
  )
)

# ============================================================
# 11. FINAL REPORT
# ============================================================

cat("\n====================================================\n")
cat("STEP7E COMPLETE — CONDITIONAL / MULTIVARIATE LAVA\n")
cat("====================================================\n")

cat("\nCandidate triads:\n")
print(CAND)

cat("\nPartial correlations:\n")
print(PCOR)

cat("\nMultiple regression:\n")
print(MULTIREG)

cat("\nAF-focused conditional summary:\n")
print(AF_SUMMARY)

cat("\nUPLOAD THESE 5 FILES:\n")
cat("1) ", CAND_FILE, "\n", sep = "")
cat("2) ", PCOR_FILE, "\n", sep = "")
cat("3) ", MULTIREG_FILE, "\n", sep = "")
cat("4) ", AF_SUMMARY_FILE, "\n", sep = "")
cat("5) ", QC_FILE, "\n", sep = "")

cat("\nAfter STEP7E is checked, proceed to locus-level colocalization.\n")
cat("====================================================\n")
