# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP1B0_HERMES_nested_audit_V2_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP 1B-0 — HERMES2 nested ZIP audit (V2 FIXED)
# AF–HFpEF–Obesity/OSA project
#
# Purpose:
#   1) Open the outer HERMES2_GWAS_HF_EUR.zip
#   2) Extract ONLY the nested inner ZIP (not all GWAS files)
#   3) List the inner ZIP contents
#   4) Automatically identify the ni-HFpEF / HFpEF file
#   5) Extract only the candidate HFpEF file(s) and preview headers
#
# This script does NOT modify your original HERMES ZIP.
# ============================================================

# =========================
# [USER EDIT: ONLY IF NEEDED]
# =========================
DATA_DIR <- "D:/A/data"

OUT_DIR <- file.path(DATA_DIR, "STEP1B0_HERMES_AUDIT")
TMP_DIR <- file.path(OUT_DIR, "tmp_inner_zip")
CAND_DIR <- file.path(OUT_DIR, "HFpEF_candidates")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CAND_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helper
# -------------------------
find_hermes_outer <- function() {
  fs <- list.files(DATA_DIR, full.names = TRUE, recursive = FALSE)
  hit <- fs[grepl("^HERMES2_GWAS_HF_EUR.*\\.zip$", basename(fs), ignore.case = TRUE)]
  if (length(hit) == 0) {
    hit <- fs[grepl("HERMES2_GWAS_HF_EUR", basename(fs), ignore.case = TRUE)]
  }
  if (length(hit) == 0) stop("Cannot find HERMES2_GWAS_HF_EUR.zip in ", DATA_DIR)
  hit[1]
}

read_first_lines <- function(path, n = 8) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
  } else {
    con <- file(path, "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

# -------------------------
# 1) Locate outer ZIP
# -------------------------
outer_zip <- find_hermes_outer()
cat("Outer HERMES ZIP:\n", outer_zip, "\n\n")

outer_list <- unzip(outer_zip, list = TRUE)
write.csv(
  outer_list,
  file.path(OUT_DIR, "STEP1B0_01_outer_zip_contents.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------
# 2) Find nested ZIP
# -------------------------
inner_names <- outer_list$Name[grepl("\\.zip$", outer_list$Name, ignore.case = TRUE)]

if (length(inner_names) == 0) {
  stop("No nested ZIP found inside the HERMES outer ZIP.")
}

cat("Nested ZIP(s) found:\n")
print(inner_names)

# If multiple ZIPs exist, prioritize the one containing HERMES2_GWAS_HF_EUR
priority <- inner_names[grepl("HERMES2_GWAS_HF_EUR", inner_names, ignore.case = TRUE)]
inner_name <- if (length(priority) > 0) priority[1] else inner_names[1]

# Extract only inner ZIP
unzip(outer_zip, files = inner_name, exdir = TMP_DIR, overwrite = TRUE)

inner_zip <- file.path(TMP_DIR, inner_name)

if (!file.exists(inner_zip)) {
  # fallback: find recursively because inner_name can contain folders
  candidates <- list.files(TMP_DIR, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) == 0) stop("Nested ZIP extraction failed.")
  inner_zip <- candidates[1]
}

cat("\nInner HERMES ZIP extracted to:\n", inner_zip, "\n\n")

# -------------------------
# 3) List inner ZIP
# -------------------------
inner_list <- unzip(inner_zip, list = TRUE)

write.csv(
  inner_list,
  file.path(OUT_DIR, "STEP1B0_02_inner_zip_contents.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Inner ZIP contains ", nrow(inner_list), " entries.\n", sep = "")

# -------------------------
# 4) Identify HFpEF candidates
# -------------------------
# IMPORTANT:
# Do NOT use character-class patterns such as [_-.] here.
# In R's default TRE regex engine, an unescaped hyphen can be interpreted
# as an invalid character range. Instead, normalize punctuation first.
#
# Examples:
#   "ni-HFpEF"  -> "nihfpef"
#   "ni_HFpEF"  -> "nihfpef"
#   "HF_pEF"    -> "hfpef"
#
normalized_names <- tolower(
  gsub("[^A-Za-z0-9]", "", inner_list$Name)
)

candidate_flag <- (
  grepl("hfpef", normalized_names, fixed = TRUE) |
  grepl("preservedejectionfraction", normalized_names, fixed = TRUE) |
  grepl("heartfailurepreserved", normalized_names, fixed = TRUE)
)

# A broader fallback: some archives may abbreviate the phenotype as PEF.
# Only use this fallback if the specific search above found nothing.
if (!any(candidate_flag)) {
  candidate_flag <- grepl("pef", normalized_names, fixed = TRUE)
}

cand <- inner_list[
  candidate_flag,
  ,
  drop = FALSE
]

write.csv(
  cand,
  file.path(OUT_DIR, "STEP1B0_03_HFpEF_candidate_files.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nHFpEF candidate entries:\n")
print(cand)

# -------------------------
# 5) Extract only likely data files
# -------------------------
if (nrow(cand) > 0) {
  data_candidates <- cand$Name[
    grepl("\\.(txt|tsv|csv|gz|bgz)$", cand$Name, ignore.case = TRUE)
  ]

  if (length(data_candidates) == 0) {
    data_candidates <- cand$Name
  }

  for (nm in data_candidates) {
    cat("\nExtracting candidate: ", nm, "\n", sep = "")
    try(
      unzip(inner_zip, files = nm, exdir = CAND_DIR, overwrite = TRUE),
      silent = TRUE
    )
  }
}

# -------------------------
# 6) Preview candidate file headers
# -------------------------
extracted_candidates <- list.files(
  CAND_DIR,
  recursive = TRUE,
  full.names = TRUE
)

preview_summary <- list()

if (length(extracted_candidates) > 0) {
  for (f in extracted_candidates) {
    if (file.info(f)$isdir) next

    lines <- tryCatch(
      read_first_lines(f, 8),
      error = function(e) paste0("READ_ERROR: ", conditionMessage(e))
    )

    preview_file <- file.path(
      OUT_DIR,
      paste0("preview_", gsub("[^A-Za-z0-9_.-]", "_", basename(f)), ".txt")
    )
    writeLines(lines, preview_file, useBytes = TRUE)

    preview_summary[[length(preview_summary) + 1]] <- data.frame(
      file = f,
      size_GB = round(file.info(f)$size / 1024^3, 4),
      first_line = ifelse(length(lines) > 0, lines[1], ""),
      stringsAsFactors = FALSE
    )
  }
}

if (length(preview_summary) > 0) {
  preview_summary <- do.call(rbind, preview_summary)
} else {
  preview_summary <- data.frame(
    file = character(0),
    size_GB = numeric(0),
    first_line = character(0)
  )
}

write.csv(
  preview_summary,
  file.path(OUT_DIR, "STEP1B0_04_HFpEF_candidate_header_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------
# 7) Also report LDSC reference archives
#    (audit only; do NOT extract huge archives here)
# -------------------------
all_top <- list.files(DATA_DIR, full.names = TRUE, recursive = FALSE)

ref_pat <- paste(
  c(
    "1000G_Phase3_baselineLD_v2.2_ldscores",
    "1000G_Phase3_frq",
    "1000G_Phase3_ldscores",
    "1000G_Phase3_plinkfiles",
    "1000G_Phase3_weights_hm3_no_MHC",
    "baselineLD_v2.2"
  ),
  collapse = "|"
)

refs <- all_top[grepl(ref_pat, basename(all_top), ignore.case = TRUE)]

ref_report <- data.frame(
  path = refs,
  basename = basename(refs),
  is_directory = file.info(refs)$isdir,
  size_GB = round(file.info(refs)$size / 1024^3, 4),
  stringsAsFactors = FALSE
)

write.csv(
  ref_report,
  file.path(OUT_DIR, "STEP1B0_05_LDSC_reference_archive_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------
# 8) README
# -------------------------
readme <- c(
  paste0("STEP 1B-0 completed: ", Sys.time()),
  "",
  "Please send back:",
  "1) STEP1B0_02_inner_zip_contents.csv",
  "2) STEP1B0_03_HFpEF_candidate_files.csv",
  "3) STEP1B0_04_HFpEF_candidate_header_summary.csv",
  "4) all preview_*.txt files",
  "5) STEP1B0_05_LDSC_reference_archive_audit.csv",
  "",
  "Do NOT delete the original HERMES2_GWAS_HF_EUR.zip.",
  "Do NOT start LDSC yet.",
  "After the exact ni-HFpEF header is confirmed, Step 1B will harmonize AF/HFpEF/BMI/OSA GWAS into analysis-ready files."
)

writeLines(
  readme,
  file.path(OUT_DIR, "STEP1B0_06_README.txt"),
  useBytes = TRUE
)

cat("\n====================================================\n")
cat("STEP 1B-0 COMPLETE\n")
cat("Output folder: ", OUT_DIR, "\n", sep = "")
cat("Please send the STEP1B0_HERMES_AUDIT folder back to ChatGPT.\n")
cat("====================================================\n")
