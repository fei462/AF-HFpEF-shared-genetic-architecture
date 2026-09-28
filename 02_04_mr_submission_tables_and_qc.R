# CODE RELEASE v1.0
# Curated final script. Original working filename: FINAL_MR_VISUALIZATION_AND_SUBMISSION_PACKAGE_V3_CORRECTED.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# FINAL MR VISUALIZATION & SUBMISSION PACKAGE — V3 FINAL
# Project: AF – HFpEF – BMI – OSA
#
# IMPORTANT CORRECTION
#   This V3 script uses ONLY the corrected STEP3B V2 results:
#       D:/A/data/STEP3_MR_V2
#   It NEVER recursively picks the older STEP3_MR results.
#
# FINAL PRESENTATION STRATEGY
#   Main Figure 5 follows the reference article:
#       Figure 5A = primary UVMR with HFpEF as outcome
#       Figure 5B = primary UVMR with AF as outcome
#
#   MVMR is placed in Supplementary Figures because OSA conditional F < 10.
#   qhet-MVMR is presented as weak-instrument robust sensitivity analysis.
#
#   No new MR is performed. This script only visualizes and packages completed
#   analyses.
#
# OUTPUT
#   D:/A/data/MR_FINAL_SUBMISSION_V3/
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

DATA_DIR <- "D:/A/data"

# ---- Exact analysis directories: DO NOT change unless your folders differ ----
STEP3B_DIR <- file.path(DATA_DIR, "STEP3_MR_V2")
STEP3C_DIR <- file.path(DATA_DIR, "STEP3_MR_STEP3C")
STEP3D_DIR <- file.path(DATA_DIR, "STEP3_MR_STEP3D")

OUT_ROOT <- file.path(DATA_DIR, "MR_FINAL_SUBMISSION_V3")

FIG_DPI <- 600
BASE_FAMILY <- "Arial"

# Reference-paper-like significance colours:
COL_BONF <- "#C65353"   # red/pink: Bonferroni significant
COL_NOM  <- "#3F78A8"   # blue: nominal
COL_NS   <- "#8B9098"   # grey: non-significant
COL_QHET <- "#7A5C43"
COL_LINE <- "#252525"
COL_REF  <- "#B8B8B8"


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (p in c("data.table", "ggplot2")) {
  install_if_missing(p)
}

library(data.table)
library(ggplot2)


# ==============================================================================
# 2. OUTPUT FOLDERS
# ==============================================================================

DIR_MAIN    <- file.path(OUT_ROOT, "01_Main_Figures")
DIR_SUPPFIG <- file.path(OUT_ROOT, "02_Supplementary_Figures")
DIR_SUPPTAB <- file.path(OUT_ROOT, "03_Supplementary_Tables")
DIR_SOURCE  <- file.path(OUT_ROOT, "04_Source_Data")
DIR_LEGEND  <- file.path(OUT_ROOT, "05_Figure_Legends")
DIR_RAW     <- file.path(OUT_ROOT, "06_Key_Analysis_Outputs")
DIR_QC      <- file.path(OUT_ROOT, "07_QC_and_Logs")

for (d in c(
  OUT_ROOT, DIR_MAIN, DIR_SUPPFIG, DIR_SUPPTAB,
  DIR_SOURCE, DIR_LEGEND, DIR_RAW, DIR_QC
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIR_QC, "MR_FINAL_V3_console.log")
if (file.exists(LOG_FILE)) unlink(LOG_FILE)

log_msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x, "\n"),
    file = LOG_FILE, append = TRUE
  )
}

log_msg("============================================================")
log_msg("MR FINAL V3 START")
log_msg("============================================================")


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

must_exist <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Required file missing [", label, "]:\n", path)
  }
  path
}

optional_file <- function(path) {
  if (file.exists(path)) path else NA_character_
}

copy_safe <- function(src, dest_dir, new_name = NULL) {
  if (length(src) != 1 || is.na(src) || !file.exists(src)) return(FALSE)
  if (is.null(new_name)) new_name <- basename(src)
  dest <- file.path(dest_dir, new_name)

  src_norm <- normalizePath(src, winslash = "/", mustWork = TRUE)
  dest_norm <- normalizePath(dest, winslash = "/", mustWork = FALSE)

  if (tolower(src_norm) == tolower(dest_norm)) return(TRUE)

  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!isTRUE(ok)) warning("Copy failed: ", src, " -> ", dest)
  ok
}

save_plot <- function(p, stub, out_dir, width = 6.6, height = 4.2) {
  ggsave(
    file.path(out_dir, paste0(stub, ".tiff")),
    p, width = width, height = height, units = "in",
    dpi = FIG_DPI, compression = "lzw"
  )
  ggsave(
    file.path(out_dir, paste0(stub, ".pdf")),
    p, width = width, height = height, units = "in",
    device = cairo_pdf
  )
}

write_legend <- function(stub, text) {
  writeLines(text, file.path(DIR_LEGEND, paste0(stub, "_legend.txt")))
}

forest_theme <- function() {
  theme_classic(base_size = 11, base_family = BASE_FAMILY) +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(8, 14, 8, 8)
    )
}

significance_colour <- function(x) {
  fifelse(
    grepl("Bonferroni", x, ignore.case = TRUE),
    COL_BONF,
    fifelse(
      grepl("Nominal", x, ignore.case = TRUE),
      COL_NOM,
      COL_NS
    )
  )
}

add_or <- function(d) {
  d <- copy(d)
  d[, OR := exp(b)]
  d[, OR_low := exp(ci_low)]
  d[, OR_high := exp(ci_high)]
  d
}

pair_label <- function(exposure, outcome) {
  paste0(exposure, " \u2192 ", outcome)
}

base_sens_tag <- function(filename) {
  x <- basename(filename)
  x <- sub("__heterogeneity\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("__Egger_intercept\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("__MRPRESSO_QC\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("__leave_one_out\\.csv$", "", x, ignore.case = TRUE)
  x
}

tag_to_traits <- function(tag) {
  core <- sub("__(primary_GWS|exploratory_P1e_6|exploratory_P5e_6)$", "", tag)
  sp <- strsplit(core, "__to__", fixed = TRUE)[[1]]
  if (length(sp) != 2) return(c(NA_character_, NA_character_))
  sp
}

tag_to_tier <- function(tag) {
  if (grepl("__primary_GWS$", tag)) return("Primary GWS")
  if (grepl("__exploratory_P1e_6$", tag)) return("Exploratory P<1e-6")
  if (grepl("__exploratory_P5e_6$", tag)) return("Exploratory P<5e-6")
  NA_character_
}


# ==============================================================================
# 4. EXACT INPUT FILES — CORRECTED V2 ONLY
# ==============================================================================

# STEP3B V2
F_UVMR_PRIMARY <- must_exist(
  file.path(STEP3B_DIR, "03_uvmr", "STEP3B_UVMR_PRIMARY_final.csv"),
  "corrected V2 UVMR primary"
)

F_UVMR_EXPLORE <- must_exist(
  file.path(STEP3B_DIR, "03_uvmr", "STEP3B_UVMR_EXPLORATORY_RELAXED_final.csv"),
  "corrected V2 UVMR exploratory"
)

F_UVMR_ALLMETHODS <- optional_file(
  file.path(STEP3B_DIR, "03_uvmr", "STEP3B_UVMR_all_methods_all_tiers.csv")
)

F_INST_QC <- must_exist(
  file.path(STEP3B_DIR, "07_reports", "STEP3B_instrument_selection_QC.csv")
)

F_LD_QC <- must_exist(
  file.path(STEP3B_DIR, "07_reports", "STEP3B_V2_LD_reference_QC.csv")
)

SENS_DIR <- file.path(STEP3B_DIR, "05_sensitivity")
if (!dir.exists(SENS_DIR)) {
  stop("Sensitivity directory missing:\n", SENS_DIR)
}

# STEP3C
F_MVMR_SELECTED <- must_exist(
  file.path(STEP3C_DIR, "03_mvmr_results", "STEP3C_SELECTED_MVMR_effects.csv")
)

F_MVMR_GRID <- must_exist(
  file.path(STEP3C_DIR, "05_reports", "STEP3C_MVMR_grid_strength_QC.csv")
)

F_STEP3C_DECISION <- optional_file(
  file.path(STEP3C_DIR, "05_reports", "STEP3C_FINAL_DECISION.txt")
)

# STEP3D
F_QHET_GRID <- must_exist(
  file.path(STEP3D_DIR, "01_results", "STEP3D_qhet_rho_grid_results.csv")
)

F_QHET_STAB <- must_exist(
  file.path(STEP3D_DIR, "01_results", "STEP3D_qhet_stability_summary.csv")
)

F_QHET_COMP <- must_exist(
  file.path(STEP3D_DIR, "01_results", "STEP3D_IVW_vs_qhet_comparison.csv")
)

F_STEP3D_DECISION <- optional_file(
  file.path(STEP3D_DIR, "03_reports", "STEP3D_FINAL_DECISION.txt")
)

F_STEP3D_LOG <- optional_file(
  file.path(STEP3D_DIR, "03_reports", "STEP3D_master_console.log")
)

log_msg("Exact corrected input paths verified.")


# ==============================================================================
# 5. HARD SAFETY CHECK — REJECT THE OLD CHR10-BIASED RESULT
# ==============================================================================

uvmr <- fread(F_UVMR_PRIMARY)
uvmr_exp <- fread(F_UVMR_EXPLORE)

expected_min_iv <- list(
  AF = 300,
  BMI = 20,
  OSA = 10
)

af_n <- uvmr[exposure == "AF" & outcome == "HFpEF", nsnp][1]
bmi_n <- min(uvmr[exposure == "BMI", nsnp], na.rm = TRUE)
osa_n <- min(uvmr[exposure == "OSA", nsnp], na.rm = TRUE)

if (
  is.na(af_n) || af_n < expected_min_iv$AF ||
  is.na(bmi_n) || bmi_n < expected_min_iv$BMI ||
  is.na(osa_n) || osa_n < expected_min_iv$OSA
) {
  stop(
    "SAFETY STOP: the loaded UVMR file does not look like the corrected STEP3B V2 result.\n",
    "Observed IV counts: AF=", af_n, ", BMI=", bmi_n, ", OSA=", osa_n, "\n",
    "Expected approximately: AF hundreds, BMI dozens, OSA tens.\n",
    "Do NOT use the older chr10-biased STEP3_MR result."
  )
}

log_msg(
  "Corrected V2 IV-count safety check PASS: AF=", af_n,
  "; BMI=", bmi_n, "; OSA=", osa_n
)

uvmr <- add_or(uvmr)
uvmr_exp <- add_or(uvmr_exp)

mvmr <- fread(F_MVMR_SELECTED)
mvmr_grid <- fread(F_MVMR_GRID)
qhet_grid <- fread(F_QHET_GRID)
qhet_stab <- fread(F_QHET_STAB)
qhet_comp <- fread(F_QHET_COMP)

if (!"OR" %in% names(mvmr)) mvmr <- add_or(mvmr)


# ==============================================================================
# 6. SENSITIVITY SUMMARY — CORRECT TAG PARSING
# ==============================================================================

het_files <- list.files(
  SENS_DIR, pattern = "__heterogeneity\\.csv$", full.names = TRUE
)
egger_files <- list.files(
  SENS_DIR, pattern = "__Egger_intercept\\.csv$", full.names = TRUE
)
presso_files <- list.files(
  SENS_DIR, pattern = "__MRPRESSO_QC\\.csv$", full.names = TRUE
)
loo_files <- list.files(
  SENS_DIR, pattern = "__leave_one_out\\.csv$", full.names = TRUE
)

# Build master tags from all sensitivity files
sens_tags <- unique(c(
  vapply(het_files, base_sens_tag, character(1)),
  vapply(egger_files, base_sens_tag, character(1)),
  vapply(presso_files, base_sens_tag, character(1)),
  vapply(loo_files, base_sens_tag, character(1))
))

sens <- data.table(tag = sens_tags)

traits <- t(vapply(
  sens$tag,
  tag_to_traits,
  FUN.VALUE = c("", "")
))

sens[, exposure := traits[,1]]
sens[, outcome := traits[,2]]
sens[, tier := vapply(tag, tag_to_tier, character(1))]
sens[, pair := pair_label(exposure, outcome)]

# Heterogeneity
if (length(het_files) > 0) {
  hd <- rbindlist(lapply(het_files, function(f) {
    x <- fread(f)
    x[, tag := base_sens_tag(f)]
    x
  }), fill = TRUE)

  # use IVW heterogeneity preferentially, otherwise first available
  hd[, is_ivw := grepl("Inverse variance weighted", method, ignore.case = TRUE)]
  hsum <- hd[
    order(-is_ivw),
    .(
      heterogeneity_Q = Q[1],
      heterogeneity_df = Q_df[1],
      heterogeneity_p = Q_pval[1]
    ),
    by = tag
  ]
  sens <- merge(sens, hsum, by = "tag", all.x = TRUE)
}

# Egger intercept
if (length(egger_files) > 0) {
  ed <- rbindlist(lapply(egger_files, function(f) {
    x <- fread(f)
    x[, tag := base_sens_tag(f)]
    x
  }), fill = TRUE)

  esum <- ed[, .(
    egger_intercept = egger_intercept[1],
    egger_intercept_se = se[1],
    egger_intercept_p = pval[1]
  ), by = tag]

  sens <- merge(sens, esum, by = "tag", all.x = TRUE)
}

# MR-PRESSO global
if (length(presso_files) > 0) {
  pd <- rbindlist(lapply(presso_files, function(f) {
    x <- fread(f)
    # Trust the filename as authoritative tag
    x[, tag := base_sens_tag(f)]
    x
  }), fill = TRUE)

  psum <- pd[, .(
    MRPRESSO_ran = MRPRESSO_ran[1],
    MRPRESSO_global_p = MRPRESSO_global_p[1],
    MRPRESSO_note = MRPRESSO_note[1]
  ), by = tag]

  sens <- merge(sens, psum, by = "tag", all.x = TRUE)
}

# Leave-one-out
loo_all <- data.table()

if (length(loo_files) > 0) {
  loo_all <- rbindlist(lapply(loo_files, function(f) {
    x <- fread(f)
    x[, tag := base_sens_tag(f)]
    x[, OR := exp(b)]
    x
  }), fill = TRUE)

  lsum <- loo_all[, .(
    loo_n = .N,
    loo_OR_min = min(OR, na.rm = TRUE),
    loo_OR_max = max(OR, na.rm = TRUE)
  ), by = tag]

  sens <- merge(sens, lsum, by = "tag", all.x = TRUE)
}

setcolorder(
  sens,
  c(
    "tag","exposure","outcome","pair","tier",
    setdiff(names(sens), c("tag","exposure","outcome","pair","tier"))
  )
)

fwrite(sens, file.path(DIR_SUPPTAB, "Table_S3_MR_sensitivity_diagnostics.csv"))

log_msg(
  "Sensitivity parser PASS: ",
  sum(!is.na(sens$egger_intercept_p)), " Egger tests; ",
  sum(!is.na(sens$MRPRESSO_global_p)), " MR-PRESSO tests; ",
  sum(!is.na(sens$heterogeneity_p)), " heterogeneity tests."
)


# ==============================================================================
# 7. MAIN FIGURE 5A — PRIMARY UVMR: HFpEF OUTCOME
# ==============================================================================

A <- copy(uvmr[outcome == "HFpEF"])
A <- A[match(c("AF","BMI","OSA"), exposure)]

A[, display := paste0(exposure, " (", nsnp, " IVs)")]
A[, display := factor(display, levels = rev(display))]
A[, point_col := significance_colour(significance_class)]

pA <- ggplot(A, aes(x = OR, y = display)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    colour = COL_REF
  ) +
  geom_segment(
    aes(x = OR_low, xend = OR_high, yend = display),
    linewidth = 0.75,
    colour = COL_LINE
  ) +
  geom_point(
    aes(fill = point_col),
    shape = 22,
    size = 3.5,
    stroke = 0.45,
    colour = COL_LINE
  ) +
  scale_fill_identity() +
  scale_x_log10() +
  labs(
    x = "Odds ratio for HFpEF (95% CI)",
    y = NULL
  ) +
  forest_theme()

save_plot(pA, "Figure5A_UVMR_HFpEF_outcome", DIR_MAIN, 6.5, 4.1)
fwrite(A, file.path(DIR_SOURCE, "Figure5A_UVMR_HFpEF_outcome_source_data.csv"))

write_legend(
  "Figure5A_UVMR_HFpEF_outcome",
  c(
    "Figure 5A | Effects of genetically predicted AF, BMI, and OSA on HFpEF.",
    "Odds ratios (squares) with 95% confidence intervals (horizontal lines) are shown from the primary MR model.",
    "For analyses with at least two instruments, inverse-variance weighted multiplicative random-effects estimates are shown; the single-instrument analysis, where applicable, uses the Wald ratio.",
    "Red/pink markers denote Bonferroni-significant associations (P < 0.05/6), blue markers denote nominal associations (P < 0.05 but not Bonferroni significant), and grey markers denote non-significant associations."
  )
)


# ==============================================================================
# 8. MAIN FIGURE 5B — PRIMARY UVMR: AF OUTCOME
# ==============================================================================

B <- copy(uvmr[outcome == "AF"])
B <- B[match(c("HFpEF","BMI","OSA"), exposure)]

B[, display := fifelse(
  exposure == "HFpEF",
  paste0("HFpEF (", nsnp, " IV; Wald ratio)"),
  paste0(exposure, " (", nsnp, " IVs)")
)]
B[, display := factor(display, levels = rev(display))]
B[, point_col := significance_colour(significance_class)]

pB <- ggplot(B, aes(x = OR, y = display)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    colour = COL_REF
  ) +
  geom_segment(
    aes(x = OR_low, xend = OR_high, yend = display),
    linewidth = 0.75,
    colour = COL_LINE
  ) +
  geom_point(
    aes(fill = point_col),
    shape = 22,
    size = 3.5,
    stroke = 0.45,
    colour = COL_LINE
  ) +
  scale_fill_identity() +
  scale_x_log10() +
  labs(
    x = "Odds ratio for AF (95% CI)",
    y = NULL
  ) +
  forest_theme()

save_plot(pB, "Figure5B_UVMR_AF_outcome", DIR_MAIN, 6.5, 4.1)
fwrite(B, file.path(DIR_SOURCE, "Figure5B_UVMR_AF_outcome_source_data.csv"))

write_legend(
  "Figure5B_UVMR_AF_outcome",
  c(
    "Figure 5B | Effects of genetically predicted HFpEF, BMI, and OSA on AF.",
    "Odds ratios (squares) with 95% confidence intervals (horizontal lines) are shown from the primary MR model.",
    "HFpEF had one genome-wide-significant instrument and therefore the primary reverse-direction estimate was calculated using the Wald ratio.",
    "Red/pink markers denote Bonferroni-significant associations (P < 0.05/6), blue markers denote nominal associations, and grey markers denote non-significant associations."
  )
)


# ==============================================================================
# 9. SUPPLEMENTARY FIGURE S1 — REVERSE HFpEF -> AF, PRIMARY + RELAXED
# ==============================================================================

S1 <- rbindlist(
  list(
    uvmr[exposure == "HFpEF" & outcome == "AF"][
      , .(
        label = "Primary GWS",
        nsnp, OR, OR_low, OR_high, pval
      )
    ],
    uvmr_exp[exposure == "HFpEF" & outcome == "AF"][
      , .(
        label = "Exploratory P<5e-6",
        nsnp, OR, OR_low, OR_high, pval
      )
    ]
  ),
  fill = TRUE
)

if (nrow(S1) > 0) {
  S1[, display := factor(
    paste0(label, " (", nsnp, " IV", ifelse(nsnp == 1, "", "s"), ")"),
    levels = rev(paste0(label, " (", nsnp, " IV", ifelse(nsnp == 1, "", "s"), ")"))
  )]

  S1[, point_col := fifelse(label == "Primary GWS", COL_NS, COL_NOM)]

  p <- ggplot(S1, aes(x = OR, y = display)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = COL_REF) +
    geom_segment(
      aes(x = OR_low, xend = OR_high, yend = display),
      linewidth = 0.75, colour = COL_LINE
    ) +
    geom_point(
      aes(fill = point_col),
      shape = 22, size = 3.4, stroke = 0.45, colour = COL_LINE
    ) +
    scale_fill_identity() +
    scale_x_log10() +
    labs(x = "Odds ratio for AF (95% CI)", y = NULL) +
    forest_theme()

  save_plot(p, "FigureS1_HFpEF_to_AF_reverse_MR", DIR_SUPPFIG, 6.3, 3.5)
  fwrite(S1, file.path(DIR_SOURCE, "FigureS1_HFpEF_to_AF_reverse_MR_source_data.csv"))

  write_legend(
    "FigureS1_HFpEF_to_AF_reverse_MR",
    c(
      "Figure S1 | Reverse-direction MR analyses for HFpEF on AF.",
      "The primary GWS analysis contained one IV and used the Wald ratio.",
      "The exploratory relaxed-threshold analysis used P < 5×10^-6 and remained non-significant."
    )
  )
}


# ==============================================================================
# 10. SUPPLEMENTARY FIGURES S2-S4 — DIAGNOSTIC P-VALUE SUMMARIES
# ==============================================================================

plot_diag <- function(
  d,
  p_col,
  fill_col,
  xlabel,
  stub
) {
  d <- copy(d[!is.na(get(p_col))])
  if (nrow(d) == 0) return(FALSE)

  d[, plot_p := pmax(get(p_col), .Machine$double.xmin)]
  d[, neglog10p := -log10(plot_p)]
  d[, display := paste0(pair, " [", tier, "]")]
  d[, display := factor(display, levels = rev(display))]

  p <- ggplot(d, aes(x = neglog10p, y = display)) +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.45,
      colour = COL_REF
    ) +
    geom_segment(
      aes(x = 0, xend = neglog10p, yend = display),
      linewidth = 0.65,
      colour = COL_LINE
    ) +
    geom_point(
      shape = 21,
      size = 2.7,
      stroke = 0.4,
      fill = fill_col,
      colour = COL_LINE
    ) +
    labs(x = xlabel, y = NULL) +
    forest_theme()

  save_plot(p, stub, DIR_SUPPFIG, 7.2, max(4.0, 0.50*nrow(d)+1.6))
  fwrite(d, file.path(DIR_SOURCE, paste0(stub, "_source_data.csv")))
  TRUE
}

plot_diag(
  sens,
  "heterogeneity_p",
  COL_BONF,
  expression(-log[10](P)~"for Cochran Q"),
  "FigureS2_CochranQ_heterogeneity"
)

plot_diag(
  sens,
  "egger_intercept_p",
  COL_NOM,
  expression(-log[10](P)~"for MR-Egger intercept"),
  "FigureS3_MR_Egger_intercept"
)

plot_diag(
  sens,
  "MRPRESSO_global_p",
  COL_QHET,
  expression(-log[10](P)~"for MR-PRESSO global test"),
  "FigureS4_MR_PRESSO_global_test"
)

write_legend(
  "FigureS2_CochranQ_heterogeneity",
  c(
    "Figure S2 | Cochran Q heterogeneity tests.",
    "The vertical dashed line corresponds to P = 0.05."
  )
)

write_legend(
  "FigureS3_MR_Egger_intercept",
  c(
    "Figure S3 | MR-Egger intercept tests for directional horizontal pleiotropy.",
    "The vertical dashed line corresponds to P = 0.05."
  )
)

write_legend(
  "FigureS4_MR_PRESSO_global_test",
  c(
    "Figure S4 | MR-PRESSO global tests for horizontal pleiotropy/outliers.",
    "The vertical dashed line corresponds to P = 0.05."
  )
)


# ==============================================================================
# 11. SUPPLEMENTARY FIGURES S5-S10 — LEAVE-ONE-OUT
# ==============================================================================

loo_order <- c(
  "AF__to__HFpEF__primary_GWS",
  "BMI__to__HFpEF__primary_GWS",
  "BMI__to__AF__primary_GWS",
  "OSA__to__HFpEF__primary_GWS",
  "OSA__to__AF__primary_GWS",
  "HFpEF__to__AF__exploratory_P5e_6"
)

if (nrow(loo_all) > 0) {

  for (ii in seq_along(loo_order)) {
    tg <- loo_order[ii]
    d <- copy(loo_all[tag == tg])
    if (nrow(d) == 0) next

    tr <- tag_to_traits(tg)
    tier <- tag_to_tier(tg)

    setorder(d, OR)
    d[, SNP_display := factor(SNP, levels = SNP)]

    p <- ggplot(d, aes(x = OR, y = SNP_display)) +
      geom_vline(
        xintercept = 1,
        linetype = "dashed",
        linewidth = 0.4,
        colour = COL_REF
      ) +
      geom_point(
        shape = 21,
        size = 1.7,
        fill = COL_NOM,
        colour = COL_LINE,
        stroke = 0.25
      ) +
      labs(
        x = paste0("Leave-one-out OR for ", pair_label(tr[1], tr[2])),
        y = "Omitted SNP"
      ) +
      theme_classic(base_size = 8.5, base_family = BASE_FAMILY) +
      theme(legend.position = "none")

    stub <- paste0(
      "FigureS", ii + 4, "_LOO_",
      tr[1], "_to_", tr[2],
      ifelse(grepl("Exploratory", tier), "_exploratory", "")
    )

    # large AF LOO plot needs more height
    hh <- min(16, max(5.0, 0.032*nrow(d) + 3.0))
    save_plot(p, stub, DIR_SUPPFIG, 7.0, hh)
    fwrite(d, file.path(DIR_SOURCE, paste0(stub, "_source_data.csv")))

    write_legend(
      stub,
      c(
        paste0("Figure S", ii + 4, " | Leave-one-out analysis for ", pair_label(tr[1], tr[2]), "."),
        paste0("Analysis tier: ", tier, "."),
        "Each point shows the causal estimate after removing one instrument."
      )
    )
  }
}


# ==============================================================================
# 12. SUPPLEMENTARY FIGURE S11 — CONDITIONAL F GRID
# ==============================================================================

# Conditional F is identical for AF and HFpEF because it evaluates exposure
# instrument strength; therefore show only one non-redundant figure.

Fgrid <- unique(
  mvmr_grid[
    outcome == "HFpEF",
    .(
      combo_id, BMI_tier, OSA_tier, relax_score,
      F_BMI, F_OSA, min_conditional_F
    )
  ]
)

setorder(Fgrid, relax_score, BMI_tier, OSA_tier)
Fgrid[, combo_label := paste0("BMI ", BMI_tier, "\nOSA ", OSA_tier)]
Fgrid[, combo_label := factor(combo_label, levels = combo_label)]

Flong <- melt(
  Fgrid,
  id.vars = c("combo_id","combo_label","relax_score"),
  measure.vars = c("F_BMI","F_OSA"),
  variable.name = "exposure",
  value.name = "conditional_F"
)

Flong[, exposure := fifelse(exposure == "F_BMI", "BMI", "OSA")]

pF <- ggplot(
  Flong,
  aes(x = combo_label, y = conditional_F, group = exposure, shape = exposure)
) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.5, colour = COL_REF) +
  geom_line(linewidth = 0.65, colour = COL_LINE) +
  geom_point(size = 2.5, fill = "white", colour = COL_LINE) +
  labs(
    x = NULL,
    y = "Conditional F-statistic"
  ) +
  theme_classic(base_size = 10, base_family = BASE_FAMILY) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

save_plot(pF, "FigureS11_MVMR_conditional_F_grid", DIR_SUPPFIG, 8.2, 4.5)
fwrite(Flong, file.path(DIR_SOURCE, "FigureS11_MVMR_conditional_F_grid_source_data.csv"))

write_legend(
  "FigureS11_MVMR_conditional_F_grid",
  c(
    "Figure S11 | Conditional instrument strength across the pre-specified BMI and OSA instrument-selection grid.",
    "The horizontal dashed line denotes the conventional conditional F = 10 threshold.",
    "No tested combination achieved conditional F >= 10 for both exposures simultaneously."
  )
)


# ==============================================================================
# 13. SUPPLEMENTARY FIGURES S12-S13 — SELECTED CONVENTIONAL MVMR
# ==============================================================================

make_mvmr_forest <- function(outcome_name, fig_num) {
  d <- copy(mvmr[outcome == outcome_name])
  d <- d[match(c("BMI","OSA"), exposure)]

  d[, display := paste0(exposure, " (", n_final, " SNPs)")]
  d[, display := factor(display, levels = rev(display))]
  d[, point_col := significance_colour(significance_class)]

  p <- ggplot(d, aes(x = OR, y = display)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = COL_REF) +
    geom_segment(
      aes(x = OR_low, xend = OR_high, yend = display),
      linewidth = 0.75, colour = COL_LINE
    ) +
    geom_point(
      aes(fill = point_col),
      shape = 22, size = 3.5, stroke = 0.45, colour = COL_LINE
    ) +
    scale_fill_identity() +
    scale_x_log10() +
    labs(
      x = paste0("Direct-effect OR for ", outcome_name, " (95% CI)"),
      y = NULL
    ) +
    forest_theme()

  stub <- paste0("FigureS", fig_num, "_Conventional_MVMR_", outcome_name)
  save_plot(p, stub, DIR_SUPPFIG, 6.3, 3.7)
  fwrite(d, file.path(DIR_SOURCE, paste0(stub, "_source_data.csv")))

  write_legend(
    stub,
    c(
      paste0("Figure S", fig_num, " | Conventional multivariable MR with ", outcome_name, " as the outcome."),
      "BMI and OSA were modeled simultaneously using the STEP3C-selected instrument set.",
      "Because the OSA conditional F-statistic remained <10, these estimates are considered exploratory and are not used as primary causal evidence."
    )
  )
}

make_mvmr_forest("HFpEF", 12)
make_mvmr_forest("AF", 13)


# ==============================================================================
# 14. SUPPLEMENTARY FIGURES S14-S17 — qhet RHO SENSITIVITY
# ==============================================================================

qsel <- qhet_grid[analysis_set == "STEP3C_SELECTED" & status == "OK"]

make_qhet_plot <- function(ex, oy, fig_num) {
  d <- copy(qsel[exposure == ex & outcome == oy])
  if (nrow(d) == 0) return(FALSE)

  ref <- qhet_comp[exposure == ex & outcome == oy]
  ref_or <- if (nrow(ref) > 0) ref$IVW_OR[1] else NA_real_

  p <- ggplot(d, aes(x = rho, y = OR)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.45, colour = COL_REF) +
    geom_line(linewidth = 0.75, colour = COL_QHET) +
    geom_point(shape = 21, size = 1.9, fill = COL_QHET, colour = COL_LINE, stroke = 0.3)

  if (is.finite(ref_or)) {
    p <- p +
      geom_hline(
        yintercept = ref_or,
        linetype = "dotted",
        linewidth = 0.5,
        colour = COL_NOM
      )
  }

  p <- p +
    labs(
      x = "Assumed BMI–OSA phenotypic correlation (rho)",
      y = paste0("qhet OR for ", ex, " \u2192 ", oy)
    ) +
    theme_classic(base_size = 11, base_family = BASE_FAMILY) +
    theme(legend.position = "none")

  stub <- paste0("FigureS", fig_num, "_qhet_", ex, "_to_", oy)
  save_plot(p, stub, DIR_SUPPFIG, 6.4, 4.0)
  fwrite(d, file.path(DIR_SOURCE, paste0(stub, "_source_data.csv")))

  write_legend(
    stub,
    c(
      paste0("Figure S", fig_num, " | Weak-instrument robust qhet-MVMR sensitivity analysis for ", ex, " -> ", oy, "."),
      "The curve shows Q-statistic-minimization estimates across the pre-specified BMI–OSA phenotypic correlation grid.",
      "The dashed line indicates OR = 1; the dotted line indicates the corresponding conventional STEP3C IVW-MVMR estimate.",
      "No single rho was selected as primary because no externally justified phenotypic correlation was supplied."
    )
  )
  TRUE
}

make_qhet_plot("BMI","HFpEF",14)
make_qhet_plot("OSA","HFpEF",15)
make_qhet_plot("BMI","AF",16)
make_qhet_plot("OSA","AF",17)


# ==============================================================================
# 15. SUPPLEMENTARY FIGURE S18 — qhet STABILITY RANGE
# ==============================================================================

S18 <- copy(qhet_stab[analysis_set == "STEP3C_SELECTED"])
S18[, pair := pair_label(exposure, outcome)]

ord <- c("BMI \u2192 HFpEF","OSA \u2192 HFpEF","BMI \u2192 AF","OSA \u2192 AF")
S18[, pair := factor(pair, levels = rev(ord))]

p18 <- ggplot(S18, aes(y = pair)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = COL_REF) +
  geom_segment(
    aes(x = OR_min, xend = OR_max, yend = pair),
    linewidth = 0.8, colour = COL_LINE
  ) +
  geom_point(
    aes(x = OR_median),
    shape = 22, size = 3.4,
    fill = COL_QHET, colour = COL_LINE, stroke = 0.4
  ) +
  scale_x_log10() +
  labs(x = "qhet OR range across rho grid", y = NULL) +
  forest_theme()

save_plot(p18, "FigureS18_qhet_stability_summary", DIR_SUPPFIG, 6.4, 4.0)
fwrite(S18, file.path(DIR_SOURCE, "FigureS18_qhet_stability_summary_source_data.csv"))

write_legend(
  "FigureS18_qhet_stability_summary",
  c(
    "Figure S18 | Stability of weak-instrument robust qhet-MVMR estimates across the phenotypic-correlation sensitivity grid.",
    "Horizontal lines show the minimum-to-maximum OR across the rho grid and squares show the median OR.",
    "All four estimates remained directionally positive across the evaluated rho range."
  )
)


# ==============================================================================
# 16. SUPPLEMENTARY TABLES
# ==============================================================================

# S1: primary UVMR
T1 <- uvmr[
  ,
  .(
    exposure, outcome, method, nsnp,
    b, se, pval, OR, OR_low, OR_high,
    mean_F, min_F, BH_FDR_q,
    significance_class
  )
]
fwrite(T1, file.path(DIR_SUPPTAB, "Table_S1_Primary_UVMR.csv"))

# S2: all MR estimators if available
if (!is.na(F_UVMR_ALLMETHODS) && file.exists(F_UVMR_ALLMETHODS)) {
  allm <- fread(F_UVMR_ALLMETHODS)

  # add ORs where estimable
  allm[, OR := exp(b)]
  allm[, OR_low := exp(b - 1.96*se)]
  allm[, OR_high := exp(b + 1.96*se)]

  T2 <- allm[
    analysis_class == "PRIMARY" & status == "OK",
    .(
      exposure_trait, outcome_trait, method, nsnp,
      b, se, pval, OR, OR_low, OR_high,
      tier, analysis_class, mean_F, min_F
    )
  ]

  fwrite(T2, file.path(DIR_SUPPTAB, "Table_S2_Primary_UVMR_all_estimators.csv"))

} else {
  writeLines(
    c(
      "STEP3B_UVMR_all_methods_all_tiers.csv was not found.",
      "Primary IVW/Wald results are available in Table S1.",
      "If needed, rerun only the final table-export portion of STEP3B V2; do not rerun GWAS preprocessing."
    ),
    file.path(DIR_SUPPTAB, "Table_S2_NOTICE.txt")
  )
}

# S3 already written: diagnostics

# S4 exploratory reverse
T4 <- uvmr_exp[
  ,
  .(
    exposure, outcome, method, nsnp,
    b, se, pval, OR, OR_low, OR_high,
    mean_F, min_F, tier, analysis_class
  )
]
fwrite(T4, file.path(DIR_SUPPTAB, "Table_S4_Exploratory_reverse_MR.csv"))

# S5 selected conventional MVMR
T5 <- mvmr[
  ,
  .(
    exposure, outcome, combo_id,
    BMI_tier, OSA_tier,
    BMI_p, OSA_p,
    n_final,
    b, se, pval,
    OR, OR_low, OR_high,
    F_BMI, F_OSA, min_conditional_F,
    Q, Q_p,
    significance_class,
    analysis_class
  )
]
fwrite(T5, file.path(DIR_SUPPTAB, "Table_S5_Conventional_MVMR_selected.csv"))

# S6 full conditional F grid
fwrite(mvmr_grid, file.path(DIR_SUPPTAB, "Table_S6_MVMR_conditional_F_grid.csv"))

# S7 qhet stability
fwrite(qhet_stab, file.path(DIR_SUPPTAB, "Table_S7_qhet_stability_summary.csv"))

# S8 full qhet grid
fwrite(qhet_grid, file.path(DIR_SUPPTAB, "Table_S8_qhet_full_rho_grid.csv"))

# S9 IVW versus qhet
fwrite(qhet_comp, file.path(DIR_SUPPTAB, "Table_S9_IVW_vs_qhet_comparison.csv"))

# S10 instrument selection QC
fwrite(fread(F_INST_QC), file.path(DIR_SUPPTAB, "Table_S10_Instrument_selection_QC.csv"))


# ==============================================================================
# 17. COPY KEY RAW OUTPUTS AND SENSITIVITY FILES
# ==============================================================================

key_files <- c(
  F_UVMR_PRIMARY, F_UVMR_EXPLORE, F_UVMR_ALLMETHODS,
  F_INST_QC, F_LD_QC,
  F_MVMR_SELECTED, F_MVMR_GRID, F_STEP3C_DECISION,
  F_QHET_GRID, F_QHET_STAB, F_QHET_COMP, F_STEP3D_DECISION, F_STEP3D_LOG
)

for (f in key_files) {
  if (!is.na(f) && file.exists(f)) copy_safe(f, DIR_RAW)
}

raw_sens_dir <- file.path(DIR_RAW, "Sensitivity_outputs")
dir.create(raw_sens_dir, recursive = TRUE, showWarnings = FALSE)

for (f in list.files(SENS_DIR, full.names = TRUE)) {
  if (file.exists(f)) copy_safe(f, raw_sens_dir)
}


# ==============================================================================
# 18. MASTER LEGEND FILE
# ==============================================================================

legend_files <- list.files(
  DIR_LEGEND,
  pattern = "_legend\\.txt$",
  full.names = TRUE
)

legend_text <- c(
  "MR FIGURE LEGENDS",
  ""
)

for (f in sort(legend_files)) {
  legend_text <- c(
    legend_text,
    readLines(f, warn = FALSE),
    "",
    "------------------------------------------------------------",
    ""
  )
}

writeLines(
  legend_text,
  file.path(DIR_LEGEND, "ALL_MR_FIGURE_LEGENDS.txt")
)


# ==============================================================================
# 19. SUBMISSION INDEX
# ==============================================================================

index <- data.table(
  item = c(
    "Figure 5A",
    "Figure 5B",
    paste0("Figure S", 1:18),
    paste0("Table S", 1:10)
  ),
  recommended_location = c(
    rep("Main manuscript", 2),
    rep("Supplementary material", 18),
    rep("Supplementary material", 10)
  )
)

fwrite(index, file.path(OUT_ROOT, "SUBMISSION_INDEX.csv"))

readme <- c(
  "MR FINAL SUBMISSION PACKAGE — V3",
  "",
  "IMPORTANT:",
  "This package was built only from the corrected STEP3B V2 results.",
  "The older chr10-biased STEP3_MR results are not used.",
  "",
  "FINAL MANUSCRIPT STRUCTURE",
  "Main Figure 5:",
  "  Figure 5A — Primary UVMR, HFpEF outcome",
  "  Figure 5B — Primary UVMR, AF outcome",
  "",
  "MVMR placement:",
  "  Conventional MVMR and qhet-MVMR are supplementary because OSA conditional F remained <10.",
  "",
  "Supplementary Figures:",
  "  S1 reverse HFpEF -> AF",
  "  S2 Cochran Q",
  "  S3 MR-Egger intercept",
  "  S4 MR-PRESSO",
  "  S5-S10 leave-one-out",
  "  S11 conditional-F grid",
  "  S12-S13 conventional MVMR",
  "  S14-S17 qhet rho sensitivity",
  "  S18 qhet stability summary",
  "",
  "Supplementary Tables:",
  "  S1 primary UVMR",
  "  S2 all primary estimators (IVW / weighted median / MR-Egger, when source file available)",
  "  S3 sensitivity diagnostics",
  "  S4 exploratory reverse MR",
  "  S5 conventional MVMR",
  "  S6 conditional-F grid",
  "  S7 qhet stability",
  "  S8 full qhet rho grid",
  "  S9 IVW vs qhet",
  "  S10 instrument-selection QC",
  "",
  "No additional MR analysis is required for the planned manuscript framework."
)

writeLines(readme, file.path(OUT_ROOT, "README_FINAL_MR.txt"))


# ==============================================================================
# 20. ZIP PACKAGE
# ==============================================================================

zip_path <- paste0(OUT_ROOT, ".zip")
if (file.exists(zip_path)) unlink(zip_path)

# Prefer zip package if installed; otherwise utils::zip
if (requireNamespace("zip", quietly = TRUE)) {
  zip::zipr(
    zipfile = zip_path,
    files = list.files(
      OUT_ROOT,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    ),
    root = OUT_ROOT
  )
} else {
  oldwd <- getwd()
  setwd(OUT_ROOT)
  utils::zip(
    zipfile = zip_path,
    files = list.files(".", recursive = TRUE, all.files = FALSE)
  )
  setwd(oldwd)
}


# ==============================================================================
# 21. FINAL QC
# ==============================================================================

# Confirm correct UVMR values made it into source data
qc_A <- fread(file.path(DIR_SOURCE, "Figure5A_UVMR_HFpEF_outcome_source_data.csv"))
qc_B <- fread(file.path(DIR_SOURCE, "Figure5B_UVMR_AF_outcome_source_data.csv"))

if (
  qc_A[exposure == "AF", nsnp][1] < 300 ||
  qc_A[exposure == "BMI", nsnp][1] < 20 ||
  qc_A[exposure == "OSA", nsnp][1] < 10
) {
  stop("FINAL SAFETY STOP: Figure 5A source data failed V2 IV-count QC.")
}

if (
  qc_B[exposure == "BMI", nsnp][1] < 20 ||
  qc_B[exposure == "OSA", nsnp][1] < 10
) {
  stop("FINAL SAFETY STOP: Figure 5B source data failed V2 IV-count QC.")
}

# Ensure key supplementary diagnostics exist
expected_supp <- c(
  "FigureS2_CochranQ_heterogeneity.tiff",
  "FigureS3_MR_Egger_intercept.tiff",
  "FigureS4_MR_PRESSO_global_test.tiff",
  "FigureS11_MVMR_conditional_F_grid.tiff",
  "FigureS18_qhet_stability_summary.tiff"
)

missing_supp <- expected_supp[
  !file.exists(file.path(DIR_SUPPFIG, expected_supp))
]

if (length(missing_supp) > 0) {
  warning(
    "Some expected supplementary figures were not generated: ",
    paste(missing_supp, collapse = ", ")
  )
}

log_msg("============================================================")
log_msg("MR FINAL V3 COMPLETE")
log_msg("Correct V2 data confirmed.")
log_msg("Output folder: ", OUT_ROOT)
log_msg("Zip: ", zip_path)
log_msg("============================================================")

cat("\nFINAL MR PACKAGE COMPLETE\n")
cat("Folder: ", OUT_ROOT, "\n", sep = "")
cat("Zip:    ", zip_path, "\n", sep = "")
