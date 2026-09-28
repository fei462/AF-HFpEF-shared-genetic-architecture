# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP6B_PREPARE_FUMA_PREDEFINED_LEADS_6PAIRS.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP6B — PREPARE PRE-DEFINED LEAD SNP FILES FOR 6 FUMA JOBS
#
# Why:
# Current FUMA caps "Maximum P-value of lead SNPs" at 1e-5.
# Our shared-locus statistic is conjFDR (<0.05), not a GWAS P-value.
# Therefore we DO NOT rescale/forge P-values.
# Instead, use the official pleioFDR lead SNPs as pre-defined lead SNPs
# and let FUMA perform standardized LD expansion + functional mapping.
#
# FUMA settings after this:
# - Upload standardized conjFDR SNP file
# - Upload matching pre-defined lead SNP file
# - Identify additional independent lead SNPs: OFF
# - Maximum lead SNP P-value: leave default 5e-8
#   (irrelevant for the supplied pre-defined lead SNPs)
# - Maximum GWAS P-value cutoff: 0.05
# - r2 independent SNP: 0.6
# - r2 lead SNP: 0.1
# - 1000G Phase 3 EUR
# - MAF 0.01
# - merge distance 250 kb
# - MAGMA OFF
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- "D:/A/data/STEP5_PLEIOFDR"

RUN_ROOT <- file.path(ROOT, "07_official_runs")
FUMA_ROOT <- file.path(ROOT, "09_FUMA_STANDARDIZED_6PAIRS")

pairs <- c(
  "AF_HFpEF",
  "AF_BMI",
  "AF_OSA",
  "HFpEF_BMI",
  "BMI_OSA",
  "HFpEF_OSA"
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
library(data.table)

pick_col <- function(nms, patterns, required = TRUE) {
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  if (required) {
    stop("Could not identify column using: ", paste(patterns, collapse = ", "))
  }
  NA_character_
}

manifest <- list()

for (pair in pairs) {

  cat("\n====================================================\n")
  cat("PAIR:", pair, "\n")
  cat("====================================================\n")

  pair_dir <- file.path(RUN_ROOT, pair)

  loci_hits <- list.files(
    pair_dir,
    pattern = "_zscore_conjfdr_0\\.05_loci\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(loci_hits) != 1) {
    stop(
      "Expected exactly one zscore loci CSV for ", pair,
      "\nFound: ", length(loci_hits),
      "\n", paste(loci_hits, collapse = "\n")
    )
  }

  infile <- loci_hits[1]
  dt <- fread(infile)
  nms <- names(dt)

  snp_col <- pick_col(nms, c("^snpid$", "^rsid$", "^snp$", "lead.*snp"))
  chr_col <- pick_col(nms, c("^chrnum$", "^chr$", "chromosome"))
  pos_col <- pick_col(nms, c("^chrpos$", "^bp$", "^pos$", "position"))

  conj_col <- pick_col(
    nms,
    c("^conjfdr_", "^min_conjfdr$"),
    required = FALSE
  )

  out <- data.table(
    rsID = as.character(dt[[snp_col]]),
    chr  = as.integer(dt[[chr_col]]),
    pos  = as.integer(dt[[pos_col]])
  )

  out <- out[
    !is.na(rsID) &
    nzchar(rsID) &
    !is.na(chr) &
    !is.na(pos)
  ]

  # Keep only autosomes and remove MHC again for consistency
  out <- out[
    chr >= 1 & chr <= 22 &
    !(chr == 6 & pos >= 26000000 & pos <= 34000000)
  ]

  out <- unique(out, by = c("rsID", "chr", "pos"))
  setorder(out, chr, pos)

  pair_out_dir <- file.path(FUMA_ROOT, pair)
  dir.create(pair_out_dir, recursive = TRUE, showWarnings = FALSE)

  outfile <- file.path(
    pair_out_dir,
    paste0("FUMA_predefined_leadSNPs_", pair, ".tsv")
  )

  fwrite(
    out,
    outfile,
    sep = "\t",
    quote = FALSE
  )

  manifest[[pair]] <- data.table(
    pair = pair,
    source_loci_file = infile,
    predefined_lead_snps = nrow(out),
    output_file = outfile
  )

  cat("Pre-defined lead SNPs:", nrow(out), "\n")
  cat("Output:", outfile, "\n")
}

manifest_dt <- rbindlist(manifest)

manifest_file <- file.path(
  FUMA_ROOT,
  "STEP6B_predefined_leadSNP_manifest.csv"
)

fwrite(manifest_dt, manifest_file)

cat("\n====================================================\n")
cat("STEP6B COMPLETE\n")
cat("====================================================\n")
print(manifest_dt)
cat("\nManifest:\n", manifest_file, "\n")
cat("====================================================\n")
