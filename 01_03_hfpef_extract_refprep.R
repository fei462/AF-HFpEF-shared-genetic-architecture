# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP1B_HFpEF_extract_and_refprep.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP 1B — HERMES ni-HFpEF extraction + header audit
# AF–HFpEF–Obesity/OSA project
#
# We now know from the official HERMES documentation:
#   Pheno1 = HF
#   Pheno2 = non-ischaemic HF
#   Pheno3 = non-ischaemic HFrEF
#   Pheno4 = non-ischaemic HFpEF (LVEF >= 50%)
#
# Therefore the exact primary HFpEF file is:
#   FORMAT-METAL_Pheno4_EUR.tsv.gz
#
# This script:
#   1) finds the already-extracted inner HERMES ZIP if present;
#   2) otherwise extracts only the inner ZIP from the outer ZIP;
#   3) extracts ONLY FORMAT-METAL_Pheno4_EUR.tsv.gz;
#   4) reads the first 10 lines and reports exact columns;
#   5) extracts LDSC reference .tgz archives into folders.
#
# It does NOT modify your original raw files.
# ============================================================

# =========================
# [USER EDIT ONLY IF PATH CHANGED]
# =========================
DATA_DIR <- "D:/A/data"

OUT_DIR <- file.path(DATA_DIR, "STEP1B_HFpEF_PREP")
HERMES_OUT <- file.path(OUT_DIR, "HERMES_HFpEF")
REF_OUT <- file.path(OUT_DIR, "LDSC_refs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(HERMES_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_OUT, recursive = TRUE, showWarnings = FALSE)

cat("STEP 1B started\n")
cat("DATA_DIR =", DATA_DIR, "\n\n")

# ------------------------------------------------------------
# helper functions
# ------------------------------------------------------------
find_first <- function(pattern, paths) {
  hit <- paths[grepl(pattern, basename(paths), ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_first_lines <- function(path, n = 10) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

# ------------------------------------------------------------
# 1) Locate outer HERMES ZIP
# ------------------------------------------------------------
top_files <- list.files(DATA_DIR, full.names = TRUE, recursive = FALSE)

outer_zip <- find_first("^HERMES2_GWAS_HF_EUR.*\\.zip$", top_files)
if (is.na(outer_zip)) {
  outer_zip <- find_first("HERMES2_GWAS_HF_EUR", top_files)
}
if (is.na(outer_zip) || !file.exists(outer_zip)) {
  stop("Cannot find the outer HERMES2_GWAS_HF_EUR ZIP in: ", DATA_DIR)
}

cat("Outer ZIP:", outer_zip, "\n")

# ------------------------------------------------------------
# 2) Locate the inner HERMES ZIP
#    First search the previous STEP1B0 audit folder.
# ------------------------------------------------------------
inner_candidates <- list.files(
  DATA_DIR,
  pattern = "HERMES2_GWAS_HF_EUR\\.zip$",
  recursive = TRUE,
  full.names = TRUE
)

# remove the top-level outer zip itself
inner_candidates <- inner_candidates[
  normalizePath(inner_candidates, winslash = "/", mustWork = FALSE) !=
    normalizePath(outer_zip, winslash = "/", mustWork = FALSE)
]

inner_zip <- NA_character_

if (length(inner_candidates) > 0) {
  # Prefer a copy inside STEP1B0_HERMES_AUDIT/tmp_inner_zip
  pr <- inner_candidates[
    grepl("STEP1B0_HERMES_AUDIT", inner_candidates, ignore.case = TRUE)
  ]
  if (length(pr) > 0) {
    inner_zip <- pr[1]
  } else {
    inner_zip <- inner_candidates[1]
  }
}

# If the inner ZIP is not already available, extract only that nested ZIP.
if (is.na(inner_zip) || !file.exists(inner_zip)) {
  outer_list <- unzip(outer_zip, list = TRUE)
  nested_names <- outer_list$Name[
    grepl("HERMES2_GWAS_HF_EUR\\.zip$", outer_list$Name, ignore.case = TRUE)
  ]
  if (length(nested_names) == 0) {
    stop("Could not find nested HERMES ZIP inside outer ZIP.")
  }

  nested_tmp <- file.path(OUT_DIR, "nested_zip")
  dir.create(nested_tmp, recursive = TRUE, showWarnings = FALSE)
  unzip(outer_zip, files = nested_names[1], exdir = nested_tmp, overwrite = TRUE)

  z <- list.files(
    nested_tmp,
    pattern = "HERMES2_GWAS_HF_EUR\\.zip$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(z) == 0) stop("Nested ZIP extraction failed.")
  inner_zip <- z[1]
}

cat("Inner ZIP:", inner_zip, "\n")

# ------------------------------------------------------------
# 3) Verify Pheno4 exists
# ------------------------------------------------------------
inner_list <- unzip(inner_zip, list = TRUE)
write.csv(
  inner_list,
  file.path(OUT_DIR, "STEP1B_01_inner_zip_contents.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

target_name <- "FORMAT-METAL_Pheno4_EUR.tsv.gz"

if (!target_name %in% inner_list$Name) {
  stop(
    "Expected ni-HFpEF file not found: ", target_name,
    "\nPlease send STEP1B_01_inner_zip_contents.csv."
  )
}

# ------------------------------------------------------------
# 4) Extract ONLY Pheno4 = ni-HFpEF
# ------------------------------------------------------------
target_path <- file.path(HERMES_OUT, target_name)

if (!file.exists(target_path)) {
  cat("Extracting only ni-HFpEF file (~312 MB compressed)...\n")
  unzip(
    inner_zip,
    files = target_name,
    exdir = HERMES_OUT,
    overwrite = TRUE
  )
} else {
  cat("ni-HFpEF file already extracted; reusing existing file.\n")
}

if (!file.exists(target_path)) {
  # fallback recursive search
  z <- list.files(
    HERMES_OUT,
    pattern = "^FORMAT-METAL_Pheno4_EUR\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(z) == 0) stop("Could not locate extracted Pheno4 file.")
  target_path <- z[1]
}

cat("HFpEF file:", target_path, "\n")

# ------------------------------------------------------------
# 5) Header / preview audit
# ------------------------------------------------------------
lines <- read_first_lines(target_path, 10)

writeLines(
  lines,
  file.path(OUT_DIR, "STEP1B_02_HFpEF_first10lines.txt"),
  useBytes = TRUE
)

header <- if (length(lines) > 0) lines[1] else ""
cols <- if (grepl("\t", header, fixed = TRUE)) {
  strsplit(header, "\t", fixed = TRUE)[[1]]
} else {
  strsplit(trimws(header), "\\s+")[[1]]
}

header_report <- data.frame(
  phenotype_id = "Pheno4",
  phenotype = "Non-ischaemic HFpEF",
  ancestry = "EUR",
  file = basename(target_path),
  size_GB = round(file.info(target_path)$size / 1024^3, 4),
  n_columns = length(cols),
  columns = paste(cols, collapse = " | "),
  stringsAsFactors = FALSE
)

write.csv(
  header_report,
  file.path(OUT_DIR, "STEP1B_03_HFpEF_header_report.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 6) Extract LDSC reference archives
# ------------------------------------------------------------
ref_archives <- c(
  "1000G_Phase3_baselineLD_v2.2_ldscores.tgz",
  "1000G_Phase3_frq.tgz",
  "1000G_Phase3_ldscores.tgz",
  "1000G_Phase3_plinkfiles.tgz",
  "1000G_Phase3_weights_hm3_no_MHC.tgz",
  "baselineLD_v2.2_bedfiles.tgz"
)

ref_results <- list()

for (nm in ref_archives) {
  src <- file.path(DATA_DIR, nm)

  if (!file.exists(src)) {
    ref_results[[length(ref_results) + 1]] <- data.frame(
      archive = nm,
      status = "MISSING",
      extracted_to = NA_character_,
      n_files = 0,
      stringsAsFactors = FALSE
    )
    next
  }

  dest <- file.path(REF_OUT, sub("\\.tgz$", "", nm))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  # Avoid re-extracting if folder already contains files.
  existing <- list.files(dest, recursive = TRUE, all.files = FALSE)

  if (length(existing) == 0) {
    cat("Extracting reference:", nm, "\n")
    ok <- tryCatch({
      utils::untar(src, exdir = dest)
      TRUE
    }, error = function(e) {
      message("Reference extraction error for ", nm, ": ", conditionMessage(e))
      FALSE
    })
  } else {
    ok <- TRUE
  }

  nfiles <- length(list.files(dest, recursive = TRUE, all.files = FALSE))

  ref_results[[length(ref_results) + 1]] <- data.frame(
    archive = nm,
    status = ifelse(ok && nfiles > 0, "OK", "FAILED"),
    extracted_to = dest,
    n_files = nfiles,
    stringsAsFactors = FALSE
  )
}

ref_report <- do.call(rbind, ref_results)

write.csv(
  ref_report,
  file.path(OUT_DIR, "STEP1B_04_LDSC_reference_extraction.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 7) Summary README
# ------------------------------------------------------------
readme <- c(
  paste0("STEP 1B completed: ", Sys.time()),
  "",
  "HERMES official phenotype mapping:",
  "Pheno1 = Heart Failure",
  "Pheno2 = Non-ischaemic HF",
  "Pheno3 = Non-ischaemic HFrEF",
  "Pheno4 = Non-ischaemic HFpEF (LVEF >= 50%)",
  "",
  paste0("Primary HFpEF file = ", target_path),
  "",
  "Please upload ONLY these small files:",
  "1) STEP1B_02_HFpEF_first10lines.txt",
  "2) STEP1B_03_HFpEF_header_report.csv",
  "3) STEP1B_04_LDSC_reference_extraction.csv",
  "",
  "Do NOT upload FORMAT-METAL_Pheno4_EUR.tsv.gz.",
  "Do NOT upload the extracted LDSC reference directories."
)

writeLines(
  readme,
  file.path(OUT_DIR, "STEP1B_05_README.txt"),
  useBytes = TRUE
)

cat("\n====================================================\n")
cat("STEP 1B COMPLETE\n")
cat("Output folder:", OUT_DIR, "\n")
cat("Primary HERMES HFpEF = Pheno4 = ni-HFpEF\n")
cat("Upload only the 3 small result files listed in STEP1B_05_README.txt\n")
cat("====================================================\n")
