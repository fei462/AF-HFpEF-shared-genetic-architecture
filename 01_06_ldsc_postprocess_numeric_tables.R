# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP1D_POSTPROCESS_V5_V2_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP 1D POSTPROCESS V5 — V2 ROBUST FIXED-STRING PARSER
# AF–HFpEF–Obesity/OSA project
#
# This script ONLY reads existing raw LDSC logs.
# It does NOT rerun LDSC and does NOT overwrite raw logs.
#
# Main fix:
# - No regex parsing of the key numeric lines.
# - Uses fixed-string splitting for:
#     Total Observed scale h2:
#     Intercept:
#     Genetic Correlation:
#     Z-score:
#     P:
# - Writes the exact raw lines used for parsing into a debug file.
# ============================================================

DATA_DIR <- "D:/A/data"
ROOT <- file.path(DATA_DIR, "STEP1D_LDSC")

MUNGE_DIR <- file.path(ROOT, "01_munged")
H2_DIR    <- file.path(ROOT, "02_h2")
RG_DIR    <- file.path(ROOT, "03_rg")
MASTER_LOG <- file.path(ROOT, "STEP1D_master_console.log")

OUT_DIR <- file.path(ROOT, "04_reports_FINAL_V2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_lines_safe <- function(path) {
  if (!file.exists(path)) return(character(0))
  readLines(path, warn = FALSE)
}

last_fixed_line <- function(lines, label) {
  x <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (length(x) == 0) return(NA_character_)
  tail(x, 1)
}

after_fixed <- function(line, label) {
  if (is.na(line) || !nzchar(line)) return(NA_character_)
  parts <- strsplit(line, label, fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  trimws(tail(parts, 1))
}

parse_value_se_fixed <- function(line, label) {
  rhs <- after_fixed(line, label)

  if (is.na(rhs) || !nzchar(rhs)) {
    return(c(value = NA_real_, se = NA_real_))
  }

  # Example rhs: 0.0479 (0.005)
  p1 <- strsplit(rhs, "(", fixed = TRUE)[[1]]

  value_txt <- trimws(p1[1])

  if (length(p1) >= 2) {
    se_txt <- strsplit(p1[2], ")", fixed = TRUE)[[1]][1]
    se_txt <- trimws(se_txt)
  } else {
    se_txt <- NA_character_
  }

  c(
    value = suppressWarnings(as.numeric(value_txt)),
    se = suppressWarnings(as.numeric(se_txt))
  )
}

parse_scalar_fixed <- function(line, label) {
  rhs <- after_fixed(line, label)
  if (is.na(rhs) || !nzchar(rhs)) return(NA_real_)

  # Keep only the first whitespace-delimited token
  tok <- strsplit(trimws(rhs), "[[:space:]]+")[[1]][1]
  suppressWarnings(as.numeric(tok))
}

first_integer_in_line <- function(line) {
  if (is.na(line) || !nzchar(line)) return(NA_real_)
  m <- regexpr("[0-9]+", line)
  if (m[1] < 0) return(NA_real_)
  suppressWarnings(as.numeric(regmatches(line, m)))
}

# ------------------------------------------------------------
# Locate trait logs robustly
# ------------------------------------------------------------
find_h2_log <- function(trait) {
  exact <- file.path(H2_DIR, paste0(trait, ".log"))
  if (file.exists(exact)) return(exact)

  all_logs <- list.files(H2_DIR, pattern = "\\.log$", full.names = TRUE)
  hit <- all_logs[grepl(trait, basename(all_logs), fixed = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_rg_log <- function(t1, t2) {
  exact <- file.path(RG_DIR, paste0(t1, "__", t2, ".log"))
  if (file.exists(exact)) return(exact)

  all_logs <- list.files(RG_DIR, pattern = "\\.log$", full.names = TRUE)
  hit <- all_logs[
    grepl(t1, basename(all_logs), fixed = TRUE) &
    grepl(t2, basename(all_logs), fixed = TRUE)
  ]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ------------------------------------------------------------
# 1) Rebuild h2 numeric table
# ------------------------------------------------------------
traits <- c("AF", "HFpEF", "BMI", "OSA")
h2_rows <- list()
h2_debug <- c()

for (tr in traits) {
  path <- find_h2_log(tr)

  if (is.na(path) || !file.exists(path)) {
    h2_rows[[tr]] <- data.frame(
      trait = tr, log_file = NA,
      h2_obs = NA, h2_se = NA, h2_z = NA, h2_p = NA,
      mean_chi2 = NA, lambda_gc = NA,
      intercept = NA, intercept_se = NA,
      ratio = NA, ratio_se = NA,
      n_regression_snps = NA,
      status = "LOG_MISSING"
    )
    h2_debug <- c(h2_debug, paste0("[", tr, "] LOG MISSING"))
    next
  }

  lines <- read_lines_safe(path)

  hline <- last_fixed_line(lines, "Total Observed scale h2:")
  hs <- parse_value_se_fixed(hline, "Total Observed scale h2:")

  mean_line <- last_fixed_line(lines, "Mean Chi^2:")
  mean_chi2 <- parse_scalar_fixed(mean_line, "Mean Chi^2:")

  lambda_line <- last_fixed_line(lines, "Lambda GC:")
  lambda_gc <- parse_scalar_fixed(lambda_line, "Lambda GC:")

  int_line <- last_fixed_line(lines, "Intercept:")
  ints <- parse_value_se_fixed(int_line, "Intercept:")

  ratio_line <- last_fixed_line(lines, "Ratio:")
  ratio_vs <- parse_value_se_fixed(ratio_line, "Ratio:")

  merge_line <- last_fixed_line(lines, "After merging with regression SNP LD,")
  n_reg <- if (!is.na(merge_line)) {
    rhs <- after_fixed(merge_line, "After merging with regression SNP LD,")
    first_integer_in_line(rhs)
  } else {
    NA_real_
  }

  h2_z <- unname(hs["value"]) / unname(hs["se"])
  h2_p <- if (is.finite(h2_z)) 2 * pnorm(-abs(h2_z)) else NA_real_

  ok <- is.finite(hs["value"]) && is.finite(hs["se"]) && hs["se"] > 0

  h2_rows[[tr]] <- data.frame(
    trait = tr,
    log_file = path,
    h2_obs = unname(hs["value"]),
    h2_se = unname(hs["se"]),
    h2_z = h2_z,
    h2_p = h2_p,
    mean_chi2 = mean_chi2,
    lambda_gc = lambda_gc,
    intercept = unname(ints["value"]),
    intercept_se = unname(ints["se"]),
    ratio = unname(ratio_vs["value"]),
    ratio_se = unname(ratio_vs["se"]),
    n_regression_snps = n_reg,
    status = ifelse(ok, "OK", "PARSE_FAIL")
  )

  h2_debug <- c(
    h2_debug,
    paste0("[", tr, "] ", path),
    paste0("H2_RAW: ", hline),
    paste0("MEAN_CHI2_RAW: ", mean_line),
    paste0("LAMBDA_RAW: ", lambda_line),
    paste0("INTERCEPT_RAW: ", int_line),
    paste0("RATIO_RAW: ", ratio_line),
    paste0("PARSED h2=", hs["value"], " se=", hs["se"],
           " z=", h2_z, " intercept=", ints["value"]),
    ""
  )
}

h2 <- do.call(rbind, h2_rows)
rownames(h2) <- NULL

write.csv(
  h2,
  file.path(OUT_DIR, "STEP1D_FINAL_h2_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  h2_debug,
  file.path(OUT_DIR, "DEBUG_h2_raw_lines_and_parsed_values.txt")
)

cat("\n================ H2 FINAL TABLE ================\n")
print(h2[, c("trait", "h2_obs", "h2_se", "h2_z",
             "mean_chi2", "intercept", "status")])
cat("================================================\n\n")

if (sum(h2$status == "OK") != 4) {
  stop(
    "STOP: only ", sum(h2$status == "OK"), "/4 h2 traits parsed.\n",
    "Please upload these two files:\n",
    file.path(OUT_DIR, "STEP1D_FINAL_h2_results.csv"), "\n",
    file.path(OUT_DIR, "DEBUG_h2_raw_lines_and_parsed_values.txt")
  )
}

# ------------------------------------------------------------
# 2) Rebuild rg numeric table
# ------------------------------------------------------------
pair_map <- data.frame(
  trait1 = c("AF", "AF", "AF", "HFpEF", "HFpEF", "BMI"),
  trait2 = c("HFpEF", "BMI", "OSA", "BMI", "OSA", "OSA"),
  stringsAsFactors = FALSE
)

rg_rows <- list()
rg_debug <- c()

for (i in seq_len(nrow(pair_map))) {
  t1 <- pair_map$trait1[i]
  t2 <- pair_map$trait2[i]
  pair <- paste0(t1, "–", t2)

  path <- find_rg_log(t1, t2)

  if (is.na(path) || !file.exists(path)) {
    rg_rows[[pair]] <- data.frame(
      trait1 = t1, trait2 = t2, pair = pair,
      log_file = NA,
      rg = NA, se = NA, z = NA, p = NA,
      gcov_intercept = NA, gcov_intercept_se = NA,
      n_valid_allele_snps = NA,
      status = "LOG_MISSING"
    )
    rg_debug <- c(rg_debug, paste0("[", pair, "] LOG MISSING"))
    next
  }

  lines <- read_lines_safe(path)

  rg_line <- last_fixed_line(lines, "Genetic Correlation:")
  rgs <- parse_value_se_fixed(rg_line, "Genetic Correlation:")

  z_line <- last_fixed_line(lines, "Z-score:")
  zval <- parse_scalar_fixed(z_line, "Z-score:")

  # Use the LAST line starting with P:
  p_idx <- grep("^[[:space:]]*P:[[:space:]]*", lines)
  p_line <- if (length(p_idx) > 0) lines[tail(p_idx, 1)] else NA_character_
  pval <- if (!is.na(p_line)) {
    # remove leading spaces first
    p_clean <- trimws(p_line)
    parse_scalar_fixed(p_clean, "P:")
  } else NA_real_

  valid_line <- last_fixed_line(lines, "SNPs with valid alleles.")
  n_valid <- if (!is.na(valid_line)) first_integer_in_line(valid_line) else NA_real_

  # Parse genetic covariance intercept ONLY inside Genetic Covariance section
  gcov_int <- NA_real_
  gcov_int_se <- NA_real_

  gcov_start <- grep("Genetic Covariance", lines, fixed = TRUE)
  gcorr_start <- grep("Genetic Correlation", lines, fixed = TRUE)

  if (length(gcov_start) > 0) {
    s <- tail(gcov_start, 1)
    later_corr <- gcorr_start[gcorr_start > s]
    e <- if (length(later_corr) > 0) min(later_corr) - 1 else min(length(lines), s + 50)
    block <- lines[s:e]

    cov_int_line <- last_fixed_line(block, "Intercept:")
    cov_ints <- parse_value_se_fixed(cov_int_line, "Intercept:")

    gcov_int <- unname(cov_ints["value"])
    gcov_int_se <- unname(cov_ints["se"])
  }

  has_error <- any(grepl("ERROR computing rg", lines, fixed = TRUE))

  ok <- !has_error &&
        is.finite(rgs["value"]) &&
        is.finite(rgs["se"]) &&
        rgs["se"] > 0 &&
        is.finite(pval)

  rg_rows[[pair]] <- data.frame(
    trait1 = t1,
    trait2 = t2,
    pair = pair,
    log_file = path,
    rg = unname(rgs["value"]),
    se = unname(rgs["se"]),
    z = zval,
    p = pval,
    gcov_intercept = gcov_int,
    gcov_intercept_se = gcov_int_se,
    n_valid_allele_snps = n_valid,
    status = ifelse(ok, "OK", "PARSE_OR_ANALYSIS_FAIL")
  )

  rg_debug <- c(
    rg_debug,
    paste0("[", pair, "] ", path),
    paste0("RG_RAW: ", rg_line),
    paste0("Z_RAW: ", z_line),
    paste0("P_RAW: ", p_line),
    paste0("VALID_ALLELES_RAW: ", valid_line),
    paste0("PARSED rg=", rgs["value"],
           " se=", rgs["se"],
           " z=", zval,
           " p=", pval,
           " status=", ifelse(ok, "OK", "FAIL")),
    ""
  )
}

rg <- do.call(rbind, rg_rows)
rownames(rg) <- NULL

# Multiple testing
rg$BH_FDR_q <- NA_real_
ok_idx <- which(rg$status == "OK" & is.finite(rg$p))

if (length(ok_idx) > 0) {
  rg$BH_FDR_q[ok_idx] <- p.adjust(rg$p[ok_idx], method = "BH")
}

rg$bonferroni_threshold <- 0.05 / 6
rg$bonferroni_significant <- with(
  rg,
  status == "OK" & is.finite(p) & p < bonferroni_threshold
)
rg$nominal_significant <- with(
  rg,
  status == "OK" & is.finite(p) & p < 0.05
)

write.csv(
  rg,
  file.path(OUT_DIR, "STEP1D_FINAL_rg_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  rg_debug,
  file.path(OUT_DIR, "DEBUG_rg_raw_lines_and_parsed_values.txt")
)

cat("\n================ RG FINAL TABLE ================\n")
print(rg[, c("pair", "rg", "se", "z", "p", "BH_FDR_q", "status")])
cat("===============================================\n\n")

# ------------------------------------------------------------
# 3) Power gate
# ------------------------------------------------------------
gate <- c(
  "STEP 1D FINAL POWER GATE",
  "========================",
  ""
)

for (i in seq_len(nrow(h2))) {
  z <- h2$h2_z[i]
  decision <- if (z >= 4) {
    "ROBUST"
  } else if (z >= 2) {
    "BORDERLINE"
  } else {
    "LOW POWER"
  }

  gate <- c(
    gate,
    sprintf(
      "%s: h2=%.6f, SE=%.6f, Z=%.3f -> %s",
      h2$trait[i],
      h2$h2_obs[i],
      h2$h2_se[i],
      z,
      decision
    )
  )
}

gate <- c(
  gate,
  "",
  paste0("Valid rg pairs: ", sum(rg$status == "OK"), "/6")
)

writeLines(
  gate,
  file.path(OUT_DIR, "STEP1D_FINAL_power_gate.txt")
)

cat("\n====================================================\n")
cat("POSTPROCESS V5 V2 COMPLETE\n")
cat("Output folder:\n", OUT_DIR, "\n")
cat("Valid h2 traits: ", sum(h2$status == "OK"), "/4\n", sep = "")
cat("Valid rg pairs:  ", sum(rg$status == "OK"), "/6\n", sep = "")
cat("====================================================\n")
