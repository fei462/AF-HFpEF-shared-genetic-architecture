# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7C_V2_RESUME_LOCAL_H2_SCREENING.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7C_V2 — ROBUST / RESUMABLE LAVA LOCAL h2 SCREENING
# AF / HFpEF / BMI / OSA
#
# FIXES vs STEP7C:
#   1) Handles unexpected run.univ() phenotype-label formats safely.
#   2) If "phen" is absent, tries:
#        - phenotype/trait columns
#        - row names
#        - locus$phenos when row counts match
#   3) If format still cannot be recovered, records the locus as
#      UNIV_FORMAT_FAIL and CONTINUES instead of stopping.
#   4) Captures LAVA print/warning messages to QC rather than flooding console.
#   5) Reuses the EXISTING STEP7C checkpoint, so completed loci are not lost.
#
# Primary univariate threshold:
#   Bonferroni = 0.05 / number of loci
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. USER SETTINGS
# ============================================================

RESUME_IF_POSSIBLE <- TRUE
CHECKPOINT_EVERY <- 25L

# ============================================================
# 1. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")

INPUT_RDATA <- file.path(
  LAVA_ROOT,
  "02_input",
  "STEP7B2_V4_LAVA_input.RData"
)

OUT_DIR <- file.path(LAVA_ROOT, "03_univariate")
TAB_DIR <- file.path(LAVA_ROOT, "07_tables")
LOG_DIR <- file.path(LAVA_ROOT, "08_logs")
SOFTWARE_DIR <- file.path(LAVA_ROOT, "00_software")

for (d in c(OUT_DIR, TAB_DIR, LOG_DIR, SOFTWARE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# IMPORTANT: same checkpoint file as STEP7C V1.
# This allows resume from the most recent completed checkpoint.
STATE_FILE <- file.path(
  OUT_DIR,
  "STEP7C_checkpoint_state.rds"
)

PHEN <- c("AF", "HFpEF", "BMI", "OSA")

# ============================================================
# 2. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
library(data.table)

if (!requireNamespace("LAVA", quietly = TRUE)) {
  stop("LAVA is not installed.")
}
library(LAVA)

lava_version <- as.character(packageVersion("LAVA"))

if (lava_version != "0.1.5") {
  warning(
    "Expected LAVA v0.1.5; detected v",
    lava_version
  )
}

# ============================================================
# 3. LOAD PROCESSED INPUT
# ============================================================

if (!file.exists(INPUT_RDATA)) {
  stop(
    "Missing processed LAVA input:\n",
    INPUT_RDATA
  )
}

load(INPUT_RDATA)

if (!exists("input") || !is.environment(input)) {
  stop("Invalid STEP7B2_V4 LAVA input object.")
}

# ============================================================
# 4. FIND OFFICIAL LOCUS FILE
# ============================================================

LOCUS_BASENAME <- "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"

search_dirs <- unique(c(
  SOFTWARE_DIR,
  LAVA_ROOT,
  dirname(system.file(package = "LAVA"))
))

search_dirs <- search_dirs[dir.exists(search_dirs)]

locus_hits <- character(0)

for (d in search_dirs) {
  hit <- list.files(
    d,
    pattern = paste0(
      "^",
      gsub("\\.", "\\\\.", LOCUS_BASENAME),
      "$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  locus_hits <- unique(c(locus_hits, hit))
}

if (!length(locus_hits)) {
  stop(
    "Cannot find official LAVA locus file:\n",
    LOCUS_BASENAME
  )
}

preferred <- locus_hits[
  grepl(
    normalizePath(LAVA_ROOT, winslash = "/", mustWork = FALSE),
    normalizePath(locus_hits, winslash = "/", mustWork = FALSE),
    fixed = TRUE
  )
]

LOC_FILE <- if (length(preferred)) preferred[1] else locus_hits[1]

cat("\nUsing locus file:\n", LOC_FILE, "\n")

# ============================================================
# 5. READ LOCI
# ============================================================

loci <- read.loci(LOC_FILE)

required_locus_cols <- c("LOC", "CHR", "START", "STOP")

if (
  is.null(loci) ||
  !all(required_locus_cols %in% names(loci))
) {
  stop("Invalid LAVA locus file.")
}

loci <- loci[
  loci$CHR %in% 1:22,
  ,
  drop = FALSE
]

n_loci <- nrow(loci)

if (n_loci < 2000) {
  stop("Unexpectedly few autosomal loci: ", n_loci)
}

UNIV_P_THRESHOLD <- 0.05 / n_loci

cat(
  "\nAutosomal loci: ", n_loci,
  "\nLocal h2 Bonferroni threshold: ",
  format(UNIV_P_THRESHOLD, scientific = TRUE, digits = 6),
  "\n",
  sep = ""
)

# ============================================================
# 6. HELPERS
# ============================================================

safe_eval <- function(expr) {

  printed <- character(0)
  warnings_vec <- character(0)
  value <- NULL
  err <- NULL

  value <- tryCatch(
    withCallingHandlers(
      {
        printed <- capture.output(
          value_inner <- eval.parent(substitute(expr)),
          type = "output"
        )
        value_inner
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
    ok = is.null(err),
    value = value,
    error = err,
    printed = printed,
    warnings = warnings_vec
  )
}

compact_messages <- function(x, max_n = 8L) {

  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  x <- unique(x)

  if (!length(x)) return("")

  if (length(x) > max_n) {
    x <- c(
      x[seq_len(max_n)],
      paste0("... [", length(x) - max_n, " more messages]")
    )
  }

  paste(x, collapse = " || ")
}

standardise_univ <- function(univ_obj, locus_phenos) {

  if (is.null(univ_obj)) {
    return(list(
      ok = FALSE,
      data = NULL,
      method = NA_character_,
      message = "run.univ returned NULL"
    ))
  }

  u <- as.data.table(univ_obj)

  if (nrow(u) == 0L) {
    return(list(
      ok = FALSE,
      data = NULL,
      method = NA_character_,
      message = paste0(
        "run.univ returned 0 rows; columns=",
        paste(names(u), collapse = ",")
      )
    ))
  }

  if ("phen" %in% names(u)) {
    u[, phenotype := as.character(phen)]
    return(list(
      ok = TRUE,
      data = u,
      method = "phen_column",
      message = ""
    ))
  }

  alias <- c(
    "phenotype",
    "trait",
    "pheno",
    "Phenotype",
    "Trait"
  )

  alt <- alias[alias %in% names(u)]

  if (length(alt)) {
    u[, phenotype := as.character(get(alt[1]))]
    return(list(
      ok = TRUE,
      data = u,
      method = paste0("alias_column:", alt[1]),
      message = ""
    ))
  }

  rn <- rownames(univ_obj)

  default_rn <- identical(
    rn,
    as.character(seq_len(nrow(u)))
  )

  if (
    length(rn) == nrow(u) &&
    !default_rn &&
    all(nzchar(rn))
  ) {
    u[, phenotype := as.character(rn)]
    return(list(
      ok = TRUE,
      data = u,
      method = "rownames",
      message = ""
    ))
  }

  if (
    length(locus_phenos) == nrow(u) &&
    nrow(u) > 0
  ) {
    u[, phenotype := as.character(locus_phenos)]
    return(list(
      ok = TRUE,
      data = u,
      method = "locus_phenos_by_row_order",
      message = paste0(
        "phen column absent; inferred from locus$phenos; original columns=",
        paste(names(univ_obj), collapse = ",")
      )
    ))
  }

  list(
    ok = FALSE,
    data = NULL,
    method = NA_character_,
    message = paste0(
      "Cannot recover phenotype labels. nrow=",
      nrow(u),
      "; length(locus$phenos)=",
      length(locus_phenos),
      "; columns=",
      paste(names(u), collapse = ",")
    )
  )
}

# ============================================================
# 7. CHECKPOINT / RESUME
# ============================================================

start_i <- 1L
UNIV_LIST <- list()
LOCUS_QC_LIST <- list()

if (
  RESUME_IF_POSSIBLE &&
  file.exists(STATE_FILE)
) {

  state <- readRDS(STATE_FILE)

  required_state <- c(
    "next_i",
    "univ_list",
    "locus_qc_list",
    "n_loci"
  )

  if (
    all(required_state %in% names(state)) &&
    identical(
      as.integer(state$n_loci),
      as.integer(n_loci)
    )
  ) {

    start_i <- as.integer(state$next_i)
    UNIV_LIST <- state$univ_list
    LOCUS_QC_LIST <- state$locus_qc_list

    cat(
      "\nRESUME checkpoint found.\n",
      "Already checkpointed through locus index ",
      start_i - 1L,
      ".\n",
      "Continuing at locus index ",
      start_i,
      " of ",
      n_loci,
      ".\n",
      sep = ""
    )
  } else {
    warning(
      "Existing checkpoint is incompatible; starting from locus 1."
    )
  }
}

# ============================================================
# 8. RUN LOCAL h2
# ============================================================

analysis_start <- Sys.time()

progress_marks <- unique(
  pmax(
    1L,
    ceiling(
      seq(0.05, 1, by = 0.05) * n_loci
    )
  )
)

if (start_i <= n_loci) {

  for (i in seq.int(start_i, n_loci)) {

    if (i %in% progress_marks) {
      cat(
        "\nProgress: ",
        i, "/", n_loci,
        " (",
        round(100 * i / n_loci, 1),
        "%)\n",
        sep = ""
      )
    }

    loc_def <- loci[i, , drop = FALSE]
    loc_id <- as.character(loc_def$LOC)

    proc <- safe_eval(
      process.locus(
        loc_def,
        input
      )
    )

    proc_msg <- compact_messages(
      c(
        proc$printed,
        proc$warnings,
        proc$error
      )
    )

    locus <- proc$value

    if (!proc$ok || is.null(locus)) {

      LOCUS_QC_LIST[[length(LOCUS_QC_LIST) + 1L]] <- data.table(
        locus = loc_id,
        chr = as.integer(loc_def$CHR),
        start = as.numeric(loc_def$START),
        stop = as.numeric(loc_def$STOP),
        process_status = "PROCESS_FAIL",
        n_snps = NA_integer_,
        n_pcs = NA_integer_,
        n_phenotypes_processed = 0L,
        phen_label_method = NA_character_,
        message = proc_msg
      )

    } else {

      ru <- safe_eval(
        run.univ(
          locus,
          phenos = locus$phenos
        )
      )

      ru_msg <- compact_messages(
        c(
          ru$printed,
          ru$warnings,
          ru$error
        )
      )

      if (!ru$ok || is.null(ru$value)) {

        LOCUS_QC_LIST[[length(LOCUS_QC_LIST) + 1L]] <- data.table(
          locus = loc_id,
          chr = as.integer(locus$chr),
          start = as.numeric(locus$start),
          stop = as.numeric(locus$stop),
          process_status = "UNIV_FAIL",
          n_snps = as.integer(locus$n.snps),
          n_pcs = as.integer(locus$K),
          n_phenotypes_processed = length(locus$phenos),
          phen_label_method = NA_character_,
          message = compact_messages(c(proc_msg, ru_msg))
        )

      } else {

        su <- standardise_univ(
          ru$value,
          locus$phenos
        )

        if (!su$ok) {

          LOCUS_QC_LIST[[length(LOCUS_QC_LIST) + 1L]] <- data.table(
            locus = loc_id,
            chr = as.integer(locus$chr),
            start = as.numeric(locus$start),
            stop = as.numeric(locus$stop),
            process_status = "UNIV_FORMAT_FAIL",
            n_snps = as.integer(locus$n.snps),
            n_pcs = as.integer(locus$K),
            n_phenotypes_processed = length(locus$phenos),
            phen_label_method = NA_character_,
            message = compact_messages(
              c(
                proc_msg,
                ru_msg,
                su$message
              )
            )
          )

        } else {

          u <- su$data

          required_univ_cols <- c(
            "phenotype",
            "h2.obs",
            "p"
          )

          if (!all(required_univ_cols %in% names(u))) {

            LOCUS_QC_LIST[[length(LOCUS_QC_LIST) + 1L]] <- data.table(
              locus = loc_id,
              chr = as.integer(locus$chr),
              start = as.numeric(locus$start),
              stop = as.numeric(locus$stop),
              process_status = "UNIV_FORMAT_FAIL",
              n_snps = as.integer(locus$n.snps),
              n_pcs = as.integer(locus$K),
              n_phenotypes_processed = length(locus$phenos),
              phen_label_method = su$method,
              message = compact_messages(
                c(
                  proc_msg,
                  ru_msg,
                  su$message,
                  paste0(
                    "Required cols missing; columns=",
                    paste(names(u), collapse = ",")
                  )
                )
              )
            )

          } else {

            u[
              ,
              `:=`(
                locus = loc_id,
                chr = as.integer(locus$chr),
                start = as.numeric(locus$start),
                stop = as.numeric(locus$stop),
                n_snps = as.integer(locus$n.snps),
                n_pcs = as.integer(locus$K),
                phen_label_method = su$method
              )
            ]

            if ("phen" %in% names(u)) {
              u[, phen := NULL]
            }

            leading_cols <- c(
              "locus",
              "chr",
              "start",
              "stop",
              "n_snps",
              "n_pcs",
              "phenotype",
              "phen_label_method"
            )

            setcolorder(
              u,
              c(
                leading_cols,
                setdiff(names(u), leading_cols)
              )
            )

            UNIV_LIST[[length(UNIV_LIST) + 1L]] <- u

            LOCUS_QC_LIST[[length(LOCUS_QC_LIST) + 1L]] <- data.table(
              locus = loc_id,
              chr = as.integer(locus$chr),
              start = as.numeric(locus$start),
              stop = as.numeric(locus$stop),
              process_status = "PASS",
              n_snps = as.integer(locus$n.snps),
              n_pcs = as.integer(locus$K),
              n_phenotypes_processed = length(locus$phenos),
              phen_label_method = su$method,
              message = compact_messages(
                c(
                  proc_msg,
                  ru_msg,
                  su$message
                )
              )
            )
          }
        }
      }
    }

    if (
      i %% CHECKPOINT_EVERY == 0L ||
      i == n_loci
    ) {

      saveRDS(
        list(
          next_i = i + 1L,
          univ_list = UNIV_LIST,
          locus_qc_list = LOCUS_QC_LIST,
          n_loci = n_loci,
          loc_file = LOC_FILE,
          univ_p_threshold = UNIV_P_THRESHOLD,
          LAVA_version = lava_version,
          script_version = "STEP7C_V2"
        ),
        STATE_FILE,
        compress = FALSE
      )
    }
  }
}

analysis_end <- Sys.time()

# ============================================================
# 9. COMBINE RESULTS
# ============================================================

UNIV <- if (length(UNIV_LIST)) {
  rbindlist(
    UNIV_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

LOCUS_QC <- if (length(LOCUS_QC_LIST)) {
  rbindlist(
    LOCUS_QC_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (!nrow(UNIV)) {
  stop("No usable local h2 results were produced.")
}

required_final_cols <- c(
  "locus",
  "phenotype",
  "h2.obs",
  "p"
)

if (!all(required_final_cols %in% names(UNIV))) {
  stop(
    "Combined local h2 output lacks expected columns: ",
    paste(
      setdiff(required_final_cols, names(UNIV)),
      collapse = ", "
    )
  )
}

UNIV[
  ,
  pass_bonferroni :=
    is.finite(p) &
    p < UNIV_P_THRESHOLD
]

UNIV[
  ,
  pass_nominal :=
    is.finite(p) &
    p < 0.05
]

# ============================================================
# 10. SAVE FULL RESULTS + QC
# ============================================================

FULL_FILE <- file.path(
  OUT_DIR,
  "STEP7C_V2_local_h2_all_loci.csv"
)

fwrite(UNIV, FULL_FILE)

LOCUS_QC_OUT <- file.path(
  OUT_DIR,
  "STEP7C_V2_locus_processing_QC.csv"
)

fwrite(LOCUS_QC, LOCUS_QC_OUT)

# ============================================================
# 11. PHENOTYPE SUMMARY
# ============================================================

PHENO_SUMMARY <- UNIV[
  ,
  .(
    loci_tested = uniqueN(locus),
    loci_bonferroni = sum(
      pass_bonferroni,
      na.rm = TRUE
    ),
    loci_nominal = sum(
      pass_nominal,
      na.rm = TRUE
    ),
    min_p = if (any(is.finite(p))) {
      min(p, na.rm = TRUE)
    } else {
      NA_real_
    },
    median_h2_obs = if (any(is.finite(h2.obs))) {
      median(h2.obs, na.rm = TRUE)
    } else {
      NA_real_
    },
    max_h2_obs = if (any(is.finite(h2.obs))) {
      max(h2.obs, na.rm = TRUE)
    } else {
      NA_real_
    }
  ),
  by = phenotype
]

PHENO_SUMMARY[
  ,
  univ_bonferroni_threshold := UNIV_P_THRESHOLD
]

PHENO_SUMMARY_FILE <- file.path(
  TAB_DIR,
  "STEP7C_V2_local_h2_phenotype_summary.csv"
)

fwrite(
  PHENO_SUMMARY,
  PHENO_SUMMARY_FILE
)

# ============================================================
# 12. PAIRWISE ELIGIBLE LOCI
# ============================================================

PAIRS <- data.table(
  trait1 = c(
    "AF", "AF", "AF",
    "HFpEF", "HFpEF",
    "BMI"
  ),
  trait2 = c(
    "HFpEF", "BMI", "OSA",
    "BMI", "OSA",
    "OSA"
  )
)

PAIR_LIST <- list()

for (i in seq_len(nrow(PAIRS))) {

  a <- PAIRS$trait1[i]
  b <- PAIRS$trait2[i]

  aa <- UNIV[
    phenotype == a &
    pass_bonferroni == TRUE,
    .(
      locus,
      chr,
      start,
      stop,
      p_trait1 = p,
      h2_trait1 = h2.obs
    )
  ]

  bb <- UNIV[
    phenotype == b &
    pass_bonferroni == TRUE,
    .(
      locus,
      p_trait2 = p,
      h2_trait2 = h2.obs
    )
  ]

  x <- merge(
    aa,
    bb,
    by = "locus",
    all = FALSE
  )

  if (nrow(x)) {
    x[
      ,
      `:=`(
        trait1 = a,
        trait2 = b
      )
    ]
  }

  PAIR_LIST[[paste(a, b, sep = "__")]] <- x
}

PAIR_ELIGIBLE <- rbindlist(
  PAIR_LIST,
  fill = TRUE,
  idcol = "pair_key"
)

PAIR_ELIGIBLE_FILE <- file.path(
  TAB_DIR,
  "STEP7C_V2_pairwise_bivariate_eligible_loci.csv"
)

fwrite(
  PAIR_ELIGIBLE,
  PAIR_ELIGIBLE_FILE
)

PAIR_SUMMARY <- rbindlist(
  lapply(
    names(PAIR_LIST),
    function(k) {

      parts <- strsplit(
        k,
        "__",
        fixed = TRUE
      )[[1]]

      n_eligible <- nrow(
        PAIR_LIST[[k]]
      )

      data.table(
        trait1 = parts[1],
        trait2 = parts[2],
        n_bivariate_eligible_loci = n_eligible,
        bivar_bonferroni_threshold = if (
          n_eligible > 0
        ) {
          0.05 / n_eligible
        } else {
          NA_real_
        }
      )
    }
  )
)

PAIR_SUMMARY_FILE <- file.path(
  TAB_DIR,
  "STEP7C_V2_pairwise_eligibility_summary.csv"
)

fwrite(
  PAIR_SUMMARY,
  PAIR_SUMMARY_FILE
)

# ============================================================
# 13. PROCESSING SUMMARY
# ============================================================

PROCESS_SUMMARY <- LOCUS_QC[
  ,
  .N,
  by = process_status
]

PROCESS_SUMMARY[
  ,
  total_input_loci := n_loci
]

PROCESS_SUMMARY[
  ,
  runtime_minutes := as.numeric(
    difftime(
      analysis_end,
      analysis_start,
      units = "mins"
    )
  )
]

PROCESS_SUMMARY_FILE <- file.path(
  TAB_DIR,
  "STEP7C_V2_processing_summary.csv"
)

fwrite(
  PROCESS_SUMMARY,
  PROCESS_SUMMARY_FILE
)

# ============================================================
# 14. PHENOTYPE-LABEL RECOVERY AUDIT
# ============================================================

LABEL_AUDIT <- LOCUS_QC[
  process_status == "PASS",
  .N,
  by = phen_label_method
][
  order(-N)
]

LABEL_AUDIT_FILE <- file.path(
  TAB_DIR,
  "STEP7C_V2_phenotype_label_recovery_audit.csv"
)

fwrite(
  LABEL_AUDIT,
  LABEL_AUDIT_FILE
)

# ============================================================
# 15. SIGNIFICANT LOCAL h2
# ============================================================

SIG_FILE <- file.path(
  OUT_DIR,
  "STEP7C_V2_local_h2_Bonferroni_significant.csv"
)

fwrite(
  UNIV[
    pass_bonferroni == TRUE
  ],
  SIG_FILE
)

# ============================================================
# 16. SESSION INFO + REMOVE CHECKPOINT ONLY AFTER SUCCESS
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP7C_V2_sessionInfo.txt"
  )
)

if (file.exists(STATE_FILE)) {
  file.remove(STATE_FILE)
}

# ============================================================
# 17. FINAL REPORT
# ============================================================

cat("\n====================================================\n")
cat("STEP7C_V2 COMPLETE — GENOME-WIDE LOCAL h2\n")
cat("====================================================\n")

cat(
  "\nBonferroni threshold: ",
  format(
    UNIV_P_THRESHOLD,
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)

cat("\nPhenotype summary:\n")
print(PHENO_SUMMARY)

cat("\nPairwise loci eligible for STEP7D local rg:\n")
print(PAIR_SUMMARY)

cat("\nLocus processing summary:\n")
print(PROCESS_SUMMARY)

cat("\nPhenotype-label recovery audit:\n")
print(LABEL_AUDIT)

cat("\nUPLOAD THESE 6 FILES:\n")
cat("1) ", PHENO_SUMMARY_FILE, "\n", sep = "")
cat("2) ", PAIR_SUMMARY_FILE, "\n", sep = "")
cat("3) ", PAIR_ELIGIBLE_FILE, "\n", sep = "")
cat("4) ", PROCESS_SUMMARY_FILE, "\n", sep = "")
cat("5) ", LABEL_AUDIT_FILE, "\n", sep = "")
cat("6) ", SIG_FILE, "\n", sep = "")

cat("\nDo NOT start STEP7D until these results are checked.\n")
cat("====================================================\n")
