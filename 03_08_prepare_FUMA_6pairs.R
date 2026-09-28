# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP6A_PREPARE_FUMA_6PAIRS_STANDARDIZED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP6A — PREPARE STANDARDIZED FUMA INPUTS FOR ALL 6 PAIRS
#
# Purpose:
#   1) Read official pleioFDR *_zscore_conjfdr_0.05_all.csv files
#   2) Keep conjFDR < 0.05 SNPs
#   3) Export one standardized FUMA SNP2GENE input per pair
#   4) Use conjFDR as the FUMA "P" column for cross-pair locus definition
#
# IMPORTANT:
#   - This script DOES NOT run FUMA automatically.
#   - Upload each generated *.tsv.gz manually to FUMA SNP2GENE.
#   - For the standardized six-pair FUMA analysis:
#       * DO NOT upload pre-defined lead SNPs
#       * Maximum lead SNP P-value = 0.05
#       * Maximum GWAS P-value = 0.05
#       * 1000G Phase 3 EUR
#       * r2 independent significant SNP = 0.6
#       * r2 lead SNP = 0.1
#       * merge distance = 250 kb
#       * MAF = 0.01
#       * GRCh37 / hg19
#       * MAGMA OFF for these truncated conjFDR subsets
#
# Project root expected:
# D:/A/data/STEP5_PLEIOFDR
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# =========================
# 0. USER PATH
# =========================
ROOT <- "D:/A/data/STEP5_PLEIOFDR"

RUN_ROOT <- file.path(ROOT, "07_official_runs")
OUT_ROOT <- file.path(ROOT, "09_FUMA_STANDARDIZED_6PAIRS")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1. PACKAGES
# =========================
pkgs <- c("data.table")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(data.table)

# =========================
# 2. PAIRS
# =========================
pairs <- c(
  "AF_HFpEF",
  "AF_BMI",
  "AF_OSA",
  "HFpEF_BMI",
  "BMI_OSA",
  "HFpEF_OSA"
)

# =========================
# 3. HELPERS
# =========================

find_zscore_all <- function(pair_dir) {
  hits <- list.files(
    pair_dir,
    pattern = "_zscore_conjfdr_0\\.05_all\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) != 1) {
    stop(
      "Expected exactly one zscore all CSV in:\n",
      pair_dir,
      "\nFound: ", length(hits),
      "\n", paste(hits, collapse = "\n")
    )
  }
  hits[1]
}

pick_col <- function(nms, patterns, required = TRUE) {
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  if (required) {
    stop(
      "Could not identify required column. Patterns: ",
      paste(patterns, collapse = ", ")
    )
  }
  NA_character_
}

# =========================
# 4. PROCESS EACH PAIR
# =========================

manifest <- list()

for (pair in pairs) {

  cat("\n====================================================\n")
  cat("PAIR:", pair, "\n")
  cat("====================================================\n")

  pair_dir <- file.path(RUN_ROOT, pair)
  if (!dir.exists(pair_dir)) {
    stop("Missing pair directory: ", pair_dir)
  }

  infile <- find_zscore_all(pair_dir)
  cat("Input:", infile, "\n")

  dt <- fread(infile)

  nms <- names(dt)

  snp_col <- pick_col(
    nms,
    c("^snpid$", "^rsid$", "^snp$")
  )

  chr_col <- pick_col(
    nms,
    c("^chrnum$", "^chr$", "chromosome")
  )

  bp_col <- pick_col(
    nms,
    c("^chrpos$", "^bp$", "^pos$", "position")
  )

  a1_col <- pick_col(
    nms,
    c("^A1$", "effect_allele"),
    required = FALSE
  )

  a2_col <- pick_col(
    nms,
    c("^A2$", "other_allele"),
    required = FALSE
  )

  conj_col <- pick_col(
    nms,
    c("^conjfdr_", "^min_conjfdr$")
  )

  # Prefer explicit conjFDR column over min_conjfdr
  explicit_conj <- grep(
    "^conjfdr_",
    nms,
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(explicit_conj)) {
    conj_col <- explicit_conj[1]
  }

  # -------------------------
  # QC and filtering
  # -------------------------
  dt[, CONJFDR := as.numeric(get(conj_col))]
  dt[, CHR := as.integer(get(chr_col))]
  dt[, BP  := as.integer(get(bp_col))]
  dt[, SNP := as.character(get(snp_col))]

  before_n <- nrow(dt)

  # The source CSV should already contain conjFDR < 0.05,
  # but we enforce it again here.
  dt <- dt[
    is.finite(CONJFDR) &
    CONJFDR < 0.05 &
    !is.na(SNP) &
    !is.na(CHR) &
    !is.na(BP)
  ]

  # Exclude MHC again for consistency with pleioFDR analysis
  # hg19 chr6:26–34 Mb
  dt <- dt[
    !(CHR == 6 & BP >= 26000000 & BP <= 34000000)
  ]

  dt <- unique(
    dt,
    by = c("SNP", "CHR", "BP")
  )

  setorder(dt, CHR, BP)

  after_n <- nrow(dt)

  # -------------------------
  # Construct FUMA input
  # -------------------------
  out <- data.table(
    SNP = dt$SNP,
    CHR = dt$CHR,
    BP  = dt$BP
  )

  if (!is.na(a1_col)) {
    out[, A1 := as.character(dt[[a1_col]])]
  }

  if (!is.na(a2_col)) {
    out[, A2 := as.character(dt[[a2_col]])]
  }

  # FUMA requires a P-value column.
  # Here P is the conjunctional FDR statistic,
  # following the study design for shared-locus definition.
  out[, P := dt$CONJFDR]

  # Ensure P is valid numeric probability
  out <- out[
    is.finite(P) & P > 0 & P < 0.05
  ]

  pair_out_dir <- file.path(OUT_ROOT, pair)
  dir.create(pair_out_dir, recursive = TRUE, showWarnings = FALSE)

  outfile <- file.path(
    pair_out_dir,
    paste0("FUMA_", pair, "_conjFDRlt0.05.tsv.gz")
  )

  fwrite(
    out,
    outfile,
    sep = "\t",
    quote = FALSE
  )

  # Save a human-readable README per pair
  readme <- c(
    paste0("PAIR: ", pair),
    "",
    "UPLOAD TO: FUMA SNP2GENE",
    "",
    "IMPORTANT:",
    "- This is a standardized conjFDR subset.",
    "- Column P contains conjFDR, not the original single-trait GWAS P-value.",
    "- Do NOT interpret MAGMA from this truncated input.",
    "- Do NOT upload pre-defined lead SNPs for the standardized six-pair run.",
    "",
    "FUMA SETTINGS:",
    "- Input build: GRCh37/hg19 (do NOT tick GRCh38)",
    "- Reference population: 1000G Phase 3 EUR",
    "- Maximum lead SNP P-value: 0.05",
    "- Maximum GWAS P-value: 0.05",
    "- Independent significant SNP r2: 0.6",
    "- Lead SNP r2: 0.1",
    "- Maximum distance to merge LD blocks: 250 kb",
    "- Minimum MAF: 0.01",
    "- Include variants from reference panel: YES",
    "- Identify additional independent lead SNPs: ON",
    "- MAGMA: OFF",
    "",
    "GENE MAPPING (recommended for comparability):",
    "- Positional mapping: ON, 10 kb",
    "- eQTL: GTEx v8 Heart Atrial Appendage",
    "- eQTL: GTEx v8 Heart Left Ventricle",
    "- eQTL: GTEx v8 Whole Blood",
    "- eQTL significant SNP-gene pairs only (FDR <= 0.05)",
    "- Chromatin interaction: Left Ventricle, Right Ventricle, Aorta",
    "",
    paste0("Input file: ", basename(outfile)),
    paste0("Number of SNPs: ", nrow(out))
  )

  writeLines(
    readme,
    file.path(pair_out_dir, paste0("README_FUMA_", pair, ".txt"))
  )

  manifest[[pair]] <- data.table(
    pair = pair,
    source_file = infile,
    source_rows = before_n,
    fuma_rows = nrow(out),
    conjfdr_column = conj_col,
    output_file = outfile
  )

  cat("Source rows :", before_n, "\n")
  cat("FUMA SNPs   :", nrow(out), "\n")
  cat("Output      :", outfile, "\n")
}

# =========================
# 5. MANIFEST
# =========================

manifest_dt <- rbindlist(manifest, fill = TRUE)

manifest_file <- file.path(
  OUT_ROOT,
  "STEP6A_FUMA_6PAIRS_manifest.csv"
)

fwrite(
  manifest_dt,
  manifest_file
)

cat("\n====================================================\n")
cat("STEP6A COMPLETE\n")
cat("====================================================\n")
print(manifest_dt)

cat("\nAll standardized FUMA files are under:\n")
cat(OUT_ROOT, "\n")

cat("\nNEXT:\n")
cat("Run all 6 pairs in FUMA using the SAME settings.\n")
cat("Download each complete FUMA result ZIP.\n")
cat("Do NOT proceed to final recurrent-locus matrix until all 6 FUMA jobs are complete.\n")
cat("====================================================\n")
