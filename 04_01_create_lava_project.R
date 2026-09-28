# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7A_CREATE_LAVA_FOLDERS.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7A — CREATE LAVA PROJECT FOLDERS
# ============================================================

rm(list = ls())

ROOT <- "D:/A/data/STEP7_LAVA"

dirs <- c(
  ROOT,
  file.path(ROOT, "00_reference"),
  file.path(ROOT, "00_reference", "UKB_v1.1"),
  file.path(ROOT, "01_sumstats"),
  file.path(ROOT, "02_input"),
  file.path(ROOT, "03_univariate"),
  file.path(ROOT, "04_bivariate"),
  file.path(ROOT, "05_multivariate"),
  file.path(ROOT, "06_figures"),
  file.path(ROOT, "07_tables"),
  file.path(ROOT, "08_logs")
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\n====================================================\n")
cat("STEP7A COMPLETE — LAVA folders created\n")
cat("====================================================\n")

for (d in dirs) {
  cat(d, "\n")
}

cat("\nIMPORTANT:\n")
cat("Extract ALL 7 UKB v1.1 reference ZIP files into:\n")
cat(file.path(ROOT, "00_reference", "UKB_v1.1"), "\n")
cat("\nDo NOT rename the extracted reference files.\n")
cat("Do NOT leave the reference files inside the ZIP archives.\n")
cat("After extraction, send me a screenshot of the folder contents\n")
cat("or run the reference QC script I will provide next.\n")
cat("====================================================\n")
