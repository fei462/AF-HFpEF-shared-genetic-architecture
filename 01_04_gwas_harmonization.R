# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP1C_GWAS_harmonization_V2_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP 1C — GWAS harmonization + HapMap3/no-MHC QC (V2 FIXED)
# AF–HFpEF–Obesity/OSA project
#
# Confirmed primary datasets:
#   AF    = GCST90624412 (European AF)
#   HFpEF = HERMES FORMAT-METAL_Pheno4_EUR.tsv.gz
#   BMI   = GIANT Locke 2015 SNP_gwas_mc_merge_nogc.tbl.uniq.gz
#   OSA   = FinnGen R9 G6_SLEEPAPNO
#
# Output:
#   Four harmonized GWAS files for LDSC preparation
#   One compact QC report to send back to ChatGPT
#
# IMPORTANT:
#   - This step does NOT run LDSC yet.
#   - It filters to HapMap3/no-MHC SNPs and MAF >= 0.05,
#     matching the cross-trait LDSC strategy of the reference paper.
#   - Raw files are never overwritten.
# ============================================================

# =========================
# [USER EDIT ONLY IF PATH CHANGED]
# =========================
DATA_DIR <- "D:/A/data"

OUT_DIR <- file.path(DATA_DIR, "STEP1C_GWAS_HARMONIZED")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# Packages
# =========================
pkgs <- c("data.table")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(data.table)

# =========================
# Fixed project settings
# =========================
OSA_N_CASE <- 38998
OSA_N_CTRL <- 336659
OSA_N_TOTAL <- OSA_N_CASE + OSA_N_CTRL   # 375657

VALID_BASES <- c("A", "C", "G", "T")

# =========================
# Helpers
# =========================
find_one <- function(pattern, recursive = TRUE) {
  fs <- list.files(DATA_DIR, full.names = TRUE, recursive = recursive)
  hit <- fs[grepl(pattern, basename(fs), ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_exact_basename <- function(name) {
  fs <- list.files(DATA_DIR, full.names = TRUE, recursive = TRUE)
  hit <- fs[basename(fs) == name]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_hm3 <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    stop("Cannot find hm3_no_MHC.list in DATA_DIR.")
  }

  x <- fread(path, header = FALSE, fill = TRUE, showProgress = FALSE)
  vals <- unique(unlist(x, use.names = FALSE))
  vals <- as.character(vals)
  vals <- vals[grepl("^rs[0-9]+$", vals)]
  unique(vals)
}

is_ambiguous <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

clean_common <- function(dt, trait) {
  n0 <- nrow(dt)

  # Uppercase alleles
  dt[, A1 := toupper(A1)]
  dt[, A2 := toupper(A2)]

  # Valid rsIDs and SNPs only
  dt <- dt[!is.na(SNP) & grepl("^rs[0-9]+$", SNP)]
  n_rsid <- nrow(dt)

  # Restrict to bi-allelic SNVs
  dt <- dt[A1 %in% VALID_BASES & A2 %in% VALID_BASES & A1 != A2]
  n_snv <- nrow(dt)

  # Remove strand-ambiguous SNPs
  dt <- dt[!is_ambiguous(A1, A2)]
  n_nonamb <- nrow(dt)

  # Numeric sanity
  dt <- dt[
    is.finite(BETA) &
    is.finite(SE) & SE > 0 &
    is.finite(P) & P > 0 & P <= 1 &
    is.finite(N) & N > 0 &
    is.finite(EAF) & EAF > 0 & EAF < 1
  ]
  n_numeric <- nrow(dt)

  # MAF >= 0.05 for LDSC, matching the reference-paper cross-trait LDSC QC
  dt[, MAF := pmin(EAF, 1 - EAF)]
  dt <- dt[MAF >= 0.05]
  n_maf <- nrow(dt)

  # HapMap3/no-MHC list
  dt <- dt[SNP %in% HM3]
  n_hm3 <- nrow(dt)

  # Remove duplicate rsIDs; keep the row with smallest P
  setorder(dt, SNP, P)
  dt <- dt[!duplicated(SNP)]
  n_unique <- nrow(dt)

  # Z is useful for QC; keep BETA/SE for later MR and direction checks
  dt[, Z := BETA / SE]

  # Sort by SNP for reproducibility
  setorder(dt, SNP)

  qc <- data.table(
    trait = trait,
    raw_rows = n0,
    valid_rsid = n_rsid,
    biallelic_SNV = n_snv,
    non_ambiguous = n_nonamb,
    numeric_valid = n_numeric,
    MAF_ge_0.05 = n_maf,
    HapMap3_noMHC = n_hm3,
    unique_final = n_unique,
    median_N = if (n_unique > 0) median(dt$N, na.rm = TRUE) else NA_real_,
    min_N = if (n_unique > 0) min(dt$N, na.rm = TRUE) else NA_real_,
    max_N = if (n_unique > 0) max(dt$N, na.rm = TRUE) else NA_real_,
    median_MAF = if (n_unique > 0) median(dt$MAF, na.rm = TRUE) else NA_real_,
    min_P = if (n_unique > 0) min(dt$P, na.rm = TRUE) else NA_real_,
    max_abs_Z = if (n_unique > 0) max(abs(dt$Z), na.rm = TRUE) else NA_real_
  )

  keep <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N", "EAF", "MAF", "Z")
  dt <- dt[, ..keep]

  list(data = dt, qc = qc)
}

save_trait <- function(obj, trait) {
  outfile <- file.path(OUT_DIR, paste0(trait, "_HapMap3_noMHC_LDSC.tsv.gz"))
  fwrite(
    obj$data,
    outfile,
    sep = "\t",
    quote = FALSE,
    compress = "gzip",
    na = "NA"
  )

  preview <- file.path(OUT_DIR, paste0(trait, "_preview_first20.tsv"))
  fwrite(head(obj$data, 20), preview, sep = "\t", quote = FALSE, na = "NA")

  cat(trait, "saved:", outfile, "\n")
}

# =========================
# 1) HapMap3/no-MHC list — ROBUST RESOLUTION
# =========================
# Windows may hide the final .txt extension, so a file displayed as
# "hm3_no_MHC.list" can actually be "hm3_no_MHC.list.txt".
# We therefore:
#   (A) search several common SNP-list names;
#   (B) if no list is found, reconstruct the HM3/no-MHC SNP set directly
#       from the extracted LDSC weights.hm3_noMHC.*.l2.ldscore.gz files.
#
all_project_files <- list.files(
  DATA_DIR,
  full.names = TRUE,
  recursive = TRUE,
  all.files = FALSE
)

preferred_hm3_names <- c(
  "hm3_no_MHC.list",
  "hm3_no_MHC.list.txt",
  "w_hm3.snplist",
  "w_hm3.snplist.txt"
)

HM3_PATH <- NA_character_

for (nm in preferred_hm3_names) {
  hit <- all_project_files[
    tolower(basename(all_project_files)) == tolower(nm)
  ]
  if (length(hit) > 0) {
    HM3_PATH <- hit[1]
    break
  }
}

# Broader fallback for oddly named text files
if (is.na(HM3_PATH)) {
  hit <- all_project_files[
    grepl("hm3", basename(all_project_files), ignore.case = TRUE) &
    grepl("MHC|snplist", basename(all_project_files), ignore.case = TRUE) &
    !grepl("\\.l2\\.|\\.annot\\.|\\.M$|\\.M_5_50$", basename(all_project_files), ignore.case = TRUE)
  ]
  if (length(hit) > 0) HM3_PATH <- hit[1]
}

if (!is.na(HM3_PATH) && file.exists(HM3_PATH)) {
  cat("Using HapMap3 SNP list:\n", HM3_PATH, "\n")
  HM3 <- read_hm3(HM3_PATH)
} else {
  cat(
    "Standalone hm3 SNP list not found.\n",
    "Reconstructing HapMap3/no-MHC SNPs from extracted LDSC weights files...\n"
  )

  weight_files <- all_project_files[
    grepl("weights.hm3", basename(all_project_files), fixed = TRUE, ignore.case = TRUE) &
    grepl("\\.l2\\.ldscore\\.gz$", basename(all_project_files), ignore.case = TRUE)
  ]

  # Extra fallback: use full path because the parent directory itself may contain the identifying name
  if (length(weight_files) == 0) {
    weight_files <- all_project_files[
      grepl("weights_hm3_no_MHC|weights.hm3_noMHC", all_project_files, ignore.case = TRUE) &
      grepl("\\.l2\\.ldscore\\.gz$", basename(all_project_files), ignore.case = TRUE)
    ]
  }

  if (length(weight_files) == 0) {
    stop(
      "Could not find either an hm3 SNP-list file OR the extracted ",
      "weights.hm3_noMHC.*.l2.ldscore.gz files.\n",
      "Please send a screenshot of D:/A/data/STEP1B_HFpEF_PREP/LDSC_refs/"
    )
  }

  cat("Found ", length(weight_files), " HM3/no-MHC LDSC weight files.\n", sep = "")

  hm3_parts <- lapply(weight_files, function(f) {
    x <- fread(f, select = "SNP", showProgress = FALSE)
    as.character(x$SNP)
  })

  HM3 <- unique(unlist(hm3_parts, use.names = FALSE))
  HM3 <- HM3[grepl("^rs[0-9]+$", HM3)]

  # Save the resolved list for reproducibility and all later steps
  HM3_PATH <- file.path(OUT_DIR, "hm3_no_MHC_resolved_from_LDSC_weights.list")
  fwrite(
    data.table(SNP = HM3),
    HM3_PATH,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  cat("Resolved HM3 list saved to:\n", HM3_PATH, "\n")
}

HM3 <- unique(as.character(HM3))
cat("HapMap3/no-MHC SNPs loaded:", length(HM3), "\n\n")

if (length(HM3) < 500000) {
  stop(
    "Resolved HapMap3/no-MHC set contains fewer than 500,000 rsIDs (n=",
    length(HM3),
    "). Stop here and inspect the reference files before GWAS filtering."
  )
}

# =========================
# 2) AF — GCST90624412
# =========================
AF_PATH <- find_one("^GCST90624412\\.tsv\\.gz$")
if (is.na(AF_PATH)) stop("GCST90624412.tsv.gz not found.")

cat("Reading AF...\n")
af <- fread(
  AF_PATH,
  select = c(
    "effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "p_value", "rs_id", "n"
  ),
  showProgress = TRUE
)

setnames(
  af,
  c("effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "p_value", "rs_id", "n"),
  c("A1", "A2", "BETA", "SE", "EAF", "P", "SNP", "N")
)

AF <- clean_common(af, "AF_GCST90624412_EUR")
rm(af); gc()
save_trait(AF, "AF_GCST90624412_EUR")

# =========================
# 3) HFpEF — HERMES Pheno4 EUR
# =========================
HF_PATH <- find_exact_basename("FORMAT-METAL_Pheno4_EUR.tsv.gz")
if (is.na(HF_PATH)) stop("FORMAT-METAL_Pheno4_EUR.tsv.gz not found.")

cat("Reading HFpEF...\n")
hf <- fread(
  HF_PATH,
  select = c("rsID", "A1", "A2", "A1_beta", "A1_freq", "se", "pval", "N_total"),
  showProgress = TRUE
)

setnames(
  hf,
  c("rsID", "A1", "A2", "A1_beta", "A1_freq", "se", "pval", "N_total"),
  c("SNP", "A1", "A2", "BETA", "EAF", "SE", "P", "N")
)

HF <- clean_common(hf, "HFpEF_HERMES_Pheno4_EUR")
rm(hf); gc()
save_trait(HF, "HFpEF_HERMES_Pheno4_EUR")

# =========================
# 4) BMI — GIANT Locke 2015
# =========================
BMI_PATH <- find_one("^SNP_gwas_mc_merge_nogc\\.tbl\\.uniq\\.gz$")
if (is.na(BMI_PATH)) stop("GIANT BMI file not found.")

cat("Reading BMI...\n")
bmi <- fread(
  BMI_PATH,
  select = c("SNP", "A1", "A2", "Freq1.Hapmap", "b", "se", "p", "N"),
  showProgress = TRUE
)

setnames(
  bmi,
  c("SNP", "A1", "A2", "Freq1.Hapmap", "b", "se", "p", "N"),
  c("SNP", "A1", "A2", "EAF", "BETA", "SE", "P", "N")
)

BMI <- clean_common(bmi, "BMI_GIANT_Locke2015_EUR")
rm(bmi); gc()
save_trait(BMI, "BMI_GIANT_Locke2015_EUR")

# =========================
# 5) OSA — FinnGen R9
# =========================
OSA_PATH <- find_one("^finngen_R9_G6_SLEEPAPNO\\.gz$")
if (is.na(OSA_PATH)) stop("FinnGen R9 OSA file not found.")

cat("Reading OSA...\n")
osa <- fread(
  OSA_PATH,
  select = c("ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"),
  showProgress = TRUE
)

# FinnGen beta is for ALT allele; therefore:
# A1 = ALT (effect allele), A2 = REF
setnames(
  osa,
  c("ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"),
  c("A2", "A1", "SNP", "P", "BETA", "SE", "EAF")
)

# If multiple rsIDs occur, retain the first valid rsID token.
osa[, SNP := sub("[,; ].*$", "", SNP)]
osa[, N := OSA_N_TOTAL]

OSA <- clean_common(osa, "OSA_FinnGen_R9_G6_SLEEPAPNO")
rm(osa); gc()
save_trait(OSA, "OSA_FinnGen_R9_G6_SLEEPAPNO")

# =========================
# 6) Combined QC
# =========================
QC <- rbindlist(list(AF$qc, HF$qc, BMI$qc, OSA$qc), fill = TRUE)

write.csv(
  QC,
  file.path(OUT_DIR, "STEP1C_01_GWAS_harmonization_QC.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Cross-trait SNP overlap
sets <- list(
  AF = AF$data$SNP,
  HFpEF = HF$data$SNP,
  BMI = BMI$data$SNP,
  OSA = OSA$data$SNP
)

pairs <- combn(names(sets), 2, simplify = FALSE)
overlap <- rbindlist(lapply(pairs, function(x) {
  data.table(
    trait1 = x[1],
    trait2 = x[2],
    n_trait1 = length(sets[[x[1]]]),
    n_trait2 = length(sets[[x[2]]]),
    overlap_SNPs = length(intersect(sets[[x[1]]], sets[[x[2]]]))
  )
}))

write.csv(
  overlap,
  file.path(OUT_DIR, "STEP1C_02_cross_trait_HapMap3_overlap.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Basic allele distribution / extreme Z audit
audit <- rbindlist(list(
  AF$data[, .(trait = "AF", n = .N,
              n_absZ_gt_40 = sum(abs(Z) > 40),
              n_P_lt_5e8 = sum(P < 5e-8))],
  HF$data[, .(trait = "HFpEF", n = .N,
              n_absZ_gt_40 = sum(abs(Z) > 40),
              n_P_lt_5e8 = sum(P < 5e-8))],
  BMI$data[, .(trait = "BMI", n = .N,
              n_absZ_gt_40 = sum(abs(Z) > 40),
              n_P_lt_5e8 = sum(P < 5e-8))],
  OSA$data[, .(trait = "OSA", n = .N,
              n_absZ_gt_40 = sum(abs(Z) > 40),
              n_P_lt_5e8 = sum(P < 5e-8))]
))

write.csv(
  audit,
  file.path(OUT_DIR, "STEP1C_03_extremeZ_and_GWS_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Save machine-readable manifest
manifest <- data.table(
  trait = c("AF", "HFpEF", "BMI", "OSA"),
  source = c(
    "GCST90624412",
    "HERMES Pheno4 EUR",
    "GIANT Locke 2015",
    "FinnGen R9 G6_SLEEPAPNO"
  ),
  standardized_file = c(
    "AF_GCST90624412_EUR_HapMap3_noMHC_LDSC.tsv.gz",
    "HFpEF_HERMES_Pheno4_EUR_HapMap3_noMHC_LDSC.tsv.gz",
    "BMI_GIANT_Locke2015_EUR_HapMap3_noMHC_LDSC.tsv.gz",
    "OSA_FinnGen_R9_G6_SLEEPAPNO_HapMap3_noMHC_LDSC.tsv.gz"
  ),
  A1_definition = c(
    "effect_allele",
    "A1 (A1_beta effect)",
    "A1 (b effect)",
    "ALT (beta effect)"
  ),
  N_definition = c(
    "per-SNP n",
    "per-SNP N_total",
    "per-SNP N",
    paste0("constant ", OSA_N_TOTAL, " = 38,998 cases + 336,659 controls")
  )
)

manifest[, HM3_source := HM3_PATH]

write.csv(
  manifest,
  file.path(OUT_DIR, "STEP1C_04_analysis_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# =========================
# 7) Final checks
# =========================
critical_fail <- any(QC$unique_final < 500000)
low_overlap <- any(overlap$overlap_SNPs < 400000)

readme <- c(
  paste0("STEP 1C completed: ", Sys.time()),
  "",
  "Primary datasets:",
  "AF    = GCST90624412 European",
  "HFpEF = HERMES Pheno4 European non-ischaemic HFpEF",
  "BMI   = GIANT Locke 2015 European",
  "OSA   = FinnGen R9 G6_SLEEPAPNO",
  "",
  paste0("Any trait with <500,000 final HapMap3 SNPs: ", critical_fail),
  paste0("Any pair with <400,000 overlapping HapMap3 SNPs: ", low_overlap),
  "",
  "Please upload ONLY these small files:",
  "STEP1C_01_GWAS_harmonization_QC.csv",
  "STEP1C_02_cross_trait_HapMap3_overlap.csv",
  "STEP1C_03_extremeZ_and_GWS_audit.csv",
  "STEP1C_04_analysis_manifest.csv",
  "",
  "Do NOT upload the four *.tsv.gz harmonized GWAS files.",
  "",
  "After these QC tables are checked, the next step is STEP 1D:",
  "munge_sumstats.py -> univariate LDSC h2 -> cross-trait LDSC rg."
)

writeLines(
  readme,
  file.path(OUT_DIR, "STEP1C_05_README.txt"),
  useBytes = TRUE
)

cat("\n====================================================\n")
cat("STEP 1C COMPLETE\n")
cat("Output folder:", OUT_DIR, "\n")
cat("Upload the four small CSV reports only.\n")
cat("====================================================\n")
