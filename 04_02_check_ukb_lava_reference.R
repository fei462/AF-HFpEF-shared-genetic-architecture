# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7A2_CHECK_UKB_REFERENCE.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7A2 — QC UKB v1.1 LAVA REFERENCE AFTER EXTRACTION
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

REF_DIR <- "D:/A/data/STEP7_LAVA/00_reference/UKB_v1.1"

if (!dir.exists(REF_DIR)) {
  stop(
    "Reference directory does not exist:\n",
    REF_DIR,
    "\nRun STEP7A_CREATE_LAVA_FOLDERS.R first."
  )
}

files <- list.files(
  REF_DIR,
  recursive = TRUE,
  full.names = TRUE
)

if (!length(files)) {
  stop(
    "UKB reference folder is empty.\n",
    "Extract all 7 UKB v1.1 archives into:\n",
    REF_DIR
  )
}

info <- data.frame(
  file = basename(files),
  path = files,
  size_MB = round(file.info(files)$size / 1024^2, 3),
  stringsAsFactors = FALSE
)

# Detect chromosome labels from filenames
detect_chr <- function(x) {
  out <- NA_integer_
  for (k in 1:22) {
    pat <- paste0("chr", k, "([^0-9]|$)")
    if (grepl(pat, x, ignore.case = TRUE)) {
      out <- k
      break
    }
  }
  out
}

info$chr <- vapply(
  info$file,
  detect_chr,
  integer(1)
)

chr_present <- sort(unique(na.omit(info$chr)))
chr_missing <- setdiff(1:22, chr_present)

cat("\n====================================================\n")
cat("UKB v1.1 LAVA REFERENCE QC\n")
cat("====================================================\n")
cat("Reference directory:\n", REF_DIR, "\n\n")
cat("Total files found:", nrow(info), "\n")
cat(
  "Total size (GB):",
  round(sum(file.info(files)$size, na.rm = TRUE) / 1024^3, 2),
  "\n"
)

cat("\nChromosomes detected from filenames:\n")
cat(paste(chr_present, collapse = ", "), "\n")

if (length(chr_missing)) {
  cat("\nWARNING — chromosomes not detected:\n")
  cat(paste(chr_missing, collapse = ", "), "\n")
} else {
  cat("\nPASS: chr1–chr22 are all represented.\n")
}

# Detect likely common prefix by checking official naming convention
official_like <- grepl(
  "lava[-_]?ukb[-_]?v1\\.1",
  info$file,
  ignore.case = TRUE
)

cat(
  "\nFiles matching 'lava-ukb-v1.1' naming pattern:",
  sum(official_like),
  "\n"
)

# Save QC table
OUT <- "D:/A/data/STEP7_LAVA/07_tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write.csv(
  info,
  file.path(OUT, "STEP7A2_UKB_reference_file_QC.csv"),
  row.names = FALSE
)

cat("\nExpected LAVA reference prefix for multi-chromosome analysis:\n")
cat(
  "D:/A/data/STEP7_LAVA/00_reference/UKB_v1.1/lava-ukb-v1.1\n"
)

cat("\nNOTE:\n")
cat(
  "The exact prefix will be confirmed from the extracted filenames.\n",
  "Do not rename the original reference files.\n",
  sep = ""
)

cat("\nQC table saved to:\n")
cat(
  file.path(OUT, "STEP7A2_UKB_reference_file_QC.csv"),
  "\n"
)

cat("====================================================\n")
