# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP1A_GWAS_inventory_audit.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP 1A — AF–HFpEF–Obesity/OSA project
# Raw data integrity + file/header audit
# Purpose:
#   1) Verify that all core GWAS / single-cell / eQTL / LDSC-reference files exist
#   2) Inspect GWAS headers WITHOUT altering raw data
#   3) List HERMES2 ZIP contents to identify the exact HFpEF phenotype file
#   4) Produce reports for the next step (exact harmonization + LDSC munging)
#
# IMPORTANT:
#   - Only edit DATA_DIR below if your folder is different.
#   - This script DOES NOT overwrite or decompress your raw data.
# ============================================================

# =========================
# [USER EDIT: ONLY HERE]
# =========================
DATA_DIR <- "D:/A/data"

# =========================
# Output folder
# =========================
OUT_DIR <- file.path(DATA_DIR, "STEP1A_QC")
PREVIEW_DIR <- file.path(OUT_DIR, "header_previews")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PREVIEW_DIR, showWarnings = FALSE, recursive = TRUE)

cat("STEP 1A started\n")
cat("DATA_DIR =", DATA_DIR, "\n\n")

if (!dir.exists(DATA_DIR)) {
  stop("DATA_DIR does not exist: ", DATA_DIR)
}

# =========================
# Helper functions
# =========================
all_files <- list.files(
  DATA_DIR,
  full.names = TRUE,
  recursive = FALSE,
  all.files = FALSE
)

basename_files <- basename(all_files)

find_one <- function(keyword) {
  idx <- grep(keyword, basename_files, fixed = TRUE, ignore.case = TRUE)
  if (length(idx) == 0) return(NA_character_)
  all_files[idx[1]]
}

status_of <- function(x) {
  if (is.na(x) || !file.exists(x)) "MISSING" else "OK"
}

size_gb <- function(x) {
  if (is.na(x) || !file.exists(x)) return(NA_real_)
  round(file.info(x)$size / 1024^3, 3)
}

read_first_lines <- function(path, n = 6) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(readLines(con, n = n, warn = FALSE))
  } else {
    con <- file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(readLines(con, n = n, warn = FALSE))
  }
}

detect_columns <- function(header_line) {
  if (length(header_line) == 0 || is.na(header_line) || header_line == "") {
    return(character(0))
  }
  if (grepl("\t", header_line, fixed = TRUE)) {
    return(strsplit(header_line, "\t", fixed = TRUE)[[1]])
  }
  # fallback: whitespace-delimited
  return(strsplit(trimws(header_line), "\\s+")[[1]])
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

# =========================
# 1. Core file inventory
# =========================
targets <- data.frame(
  role = c(
    "AF primary European GWAS",
    "AF European excluding UKB (sensitivity)",
    "AF UKB-only (sensitivity)",
    "AF FinnGen replication",
    "HFpEF cross-ancestry GWAS (secondary)",
    "HFpEF BBJ GWAS (ancestry-specific secondary)",
    "HFpEF MTAG (secondary; not primary LDSC)",
    "HFpEF/HF HERMES2 European package",
    "OSA FinnGen GWAS",
    "BMI GIANT Locke 2015 GWAS",
    "AF LAA sn-multiome raw/processed TAR",
    "AF LAA sn-multiome metadata",
    "AF snRNA expression matrix",
    "AF snRNA metadata",
    "AF snRNA barcodes",
    "AF snRNA genes",
    "HFpEF snRNA H5AD",
    "GTEx v8 cis-eQTL",
    "SMR cis-eQTL resource"
  ),
  keyword = c(
    "GCST90624412",
    "GCST90624413",
    "GCST90624414",
    "finngen_R12_I9_AF",
    "GCST90654629",
    "GCST90668011",
    "GCST90668015",
    "HERMES2_GWAS_HF_EUR",
    "finngen_R9_G6_SLEEPAPNO",
    "SNP_gwas_mc_merge_nogc.tbl.uniq",
    "GSE238242_RAW",
    "GSE238242_snAF.metadata",
    "GSE255612_AF_snRNA_Matrix",
    "GSE255612_AF_snRNA_MetaData",
    "GSE255612_AF_snRNA_Processed_Expression_Matrix_barcodes",
    "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes",
    "HFpEF_snRNAseq_single_cell_portal",
    "GTEx_V8_cis_eqtl_summary_lite",
    "cis-eQTL-SMR_20191212"
  ),
  stringsAsFactors = FALSE
)

targets$path <- vapply(targets$keyword, find_one, character(1))
targets$status <- vapply(targets$path, status_of, character(1))
targets$size_GB <- vapply(targets$path, size_gb, numeric(1))

write.csv(
  targets,
  file.path(OUT_DIR, "STEP1A_01_core_file_inventory.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Core inventory written.\n")
print(targets[, c("role", "status", "size_GB")], row.names = FALSE)

# =========================
# 2. LDSC reference inventory
# =========================
ref_dirs <- c(
  "1000G_Phase3_baselineLD_v2.2_ldscores",
  "1000G_Phase3_frq",
  "1000G_Phase3_ldscores",
  "1000G_Phase3_plinkfiles",
  "1000G_Phase3_weights_hm3_no_MHC",
  "baselineLD_v2.2_bedfiles"
)

ref_report <- lapply(ref_dirs, function(d) {
  p <- file.path(DATA_DIR, d)
  n <- if (dir.exists(p)) length(list.files(p, recursive = TRUE)) else 0
  data.frame(
    resource = d,
    exists = dir.exists(p),
    n_files = n,
    stringsAsFactors = FALSE
  )
})
ref_report <- do.call(rbind, ref_report)

extra_ref <- data.frame(
  resource = c("hm3_no_MHC.list"),
  exists = c(file.exists(file.path(DATA_DIR, "hm3_no_MHC.list"))),
  n_files = c(ifelse(file.exists(file.path(DATA_DIR, "hm3_no_MHC.list")), 1, 0))
)

ref_report <- rbind(ref_report, extra_ref)

write.csv(
  ref_report,
  file.path(OUT_DIR, "STEP1A_02_LDSC_reference_inventory.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# =========================
# 3. GWAS header audit
#    We do NOT read whole files yet.
# =========================
gwas_keywords <- c(
  "GCST90624412",
  "GCST90624413",
  "GCST90624414",
  "finngen_R12_I9_AF",
  "GCST90654629",
  "GCST90668011",
  "GCST90668015",
  "finngen_R9_G6_SLEEPAPNO",
  "SNP_gwas_mc_merge_nogc.tbl.uniq"
)

header_rows <- list()

for (kw in gwas_keywords) {
  p <- find_one(kw)
  if (is.na(p) || !file.exists(p)) {
    header_rows[[length(header_rows) + 1]] <- data.frame(
      dataset = kw,
      file = NA_character_,
      status = "MISSING",
      n_detected_columns = NA_integer_,
      columns = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  lines <- tryCatch(
    read_first_lines(p, n = 6),
    error = function(e) paste0("ERROR: ", conditionMessage(e))
  )

  preview_file <- file.path(PREVIEW_DIR, paste0(safe_name(kw), "_first6lines.txt"))
  writeLines(lines, preview_file, useBytes = TRUE)

  if (length(lines) > 0 && !startsWith(lines[1], "ERROR:")) {
    cols <- detect_columns(lines[1])
    header_rows[[length(header_rows) + 1]] <- data.frame(
      dataset = kw,
      file = basename(p),
      status = "OK",
      n_detected_columns = length(cols),
      columns = paste(cols, collapse = " | "),
      stringsAsFactors = FALSE
    )
  } else {
    header_rows[[length(header_rows) + 1]] <- data.frame(
      dataset = kw,
      file = basename(p),
      status = "READ_ERROR",
      n_detected_columns = NA_integer_,
      columns = paste(lines, collapse = " "),
      stringsAsFactors = FALSE
    )
  }
}

header_report <- do.call(rbind, header_rows)

write.csv(
  header_report,
  file.path(OUT_DIR, "STEP1A_03_GWAS_header_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# =========================
# 4. HERMES2 ZIP content audit
# =========================
hermes_zip <- find_one("HERMES2_GWAS_HF_EUR")

if (!is.na(hermes_zip) && file.exists(hermes_zip) &&
    grepl("\\.zip$", hermes_zip, ignore.case = TRUE)) {

  z <- tryCatch(
    unzip(hermes_zip, list = TRUE),
    error = function(e) NULL
  )

  if (!is.null(z)) {
    write.csv(
      z,
      file.path(OUT_DIR, "STEP1A_04_HERMES2_zip_contents.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    candidate <- z[grepl(
      "HFpEF|HFPEF|preserved|Pheno",
      z$Name,
      ignore.case = TRUE
    ), , drop = FALSE]

    write.csv(
      candidate,
      file.path(OUT_DIR, "STEP1A_05_HERMES2_candidate_HFpEF_files.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
}

# =========================
# 5. GSE238242 TAR content audit
# =========================
gse_tar <- find_one("GSE238242_RAW")

if (!is.na(gse_tar) && file.exists(gse_tar) &&
    grepl("\\.tar$", gse_tar, ignore.case = TRUE)) {
  tar_names <- tryCatch(
    utils::untar(gse_tar, list = TRUE),
    error = function(e) character(0)
  )
  writeLines(
    tar_names,
    file.path(OUT_DIR, "STEP1A_06_GSE238242_RAW_tar_contents.txt"),
    useBytes = TRUE
  )
}

# =========================
# 6. Basic project readiness summary
# =========================
critical_roles <- c(
  "AF primary European GWAS",
  "HFpEF/HF HERMES2 European package",
  "OSA FinnGen GWAS",
  "BMI GIANT Locke 2015 GWAS",
  "AF LAA sn-multiome raw/processed TAR",
  "AF LAA sn-multiome metadata",
  "AF snRNA expression matrix",
  "AF snRNA metadata",
  "HFpEF snRNA H5AD",
  "GTEx v8 cis-eQTL",
  "SMR cis-eQTL resource"
)

critical <- targets[targets$role %in% critical_roles, ]
ready <- all(critical$status == "OK")

summary_txt <- c(
  paste0("STEP 1A completed: ", Sys.time()),
  paste0("Core critical files present: ", sum(critical$status == "OK"), "/", nrow(critical)),
  paste0("Ready to proceed to exact GWAS harmonization: ", ifelse(ready, "YES", "CHECK MISSING ITEMS")),
  "",
  "IMPORTANT:",
  "1) For European LDSC, use a EUROPEAN HFpEF GWAS from HERMES2 as primary.",
  "2) GCST90654629 is cross-ancestry and GCST90668011 is BBJ; keep these for secondary/replication analyses.",
  "3) GCST90668015 is MTAG; do not use it as the sole primary HFpEF phenotype for cross-trait LDSC.",
  "4) OSA FinnGen R9 is usable for pipeline testing; a newer FinnGen release can be substituted before final analysis.",
  "",
  "Please send back these files:",
  "STEP1A_03_GWAS_header_audit.csv",
  "STEP1A_04_HERMES2_zip_contents.csv",
  "STEP1A_05_HERMES2_candidate_HFpEF_files.csv",
  "and, if convenient, screenshots of the STEP1A_QC folder."
)

writeLines(
  summary_txt,
  file.path(OUT_DIR, "STEP1A_07_README_results.txt"),
  useBytes = TRUE
)

cat("\n========================================\n")
cat("STEP 1A COMPLETE\n")
cat("Output folder:", OUT_DIR, "\n")
cat("Ready status:", ifelse(ready, "YES", "CHECK MISSING ITEMS"), "\n")
cat("========================================\n")
