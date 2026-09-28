# CODE RELEASE v1.0
# Derived from the validated full-GWAS standardization stage; MiXeR modeling is not used in the final manuscript.

# ============================================================
# STEP 5A0 — Full-GWAS standardized input preparation for pleioFDR/conjFDR
# AF–HFpEF–BMI–OSA project
#
# IMPORTANT:
# 1) This step does NOT run pleioFDR yet.
# 2) It creates standardized FULL-GWAS inputs from the raw GWAS,
#    rather than reusing the STEP1C HapMap3-only files.
# 3) For case-control traits, N is converted to EFFECTIVE sample size,
#    for balanced case-control scale handling used in the project.
# 4) Processing is performed with DuckDB to avoid loading 20–36 million
#    rows into R memory at once.
#
# Output columns:
#   SNP  A1  A2  N  Z
#
# Output:
#   D:/A/data/STEP5_PLEIOFDR/
#       01_inputs/
#       02_qc/
# ============================================================

# =========================
# [USER EDIT ONLY IF PATH CHANGED]
# =========================
DATA_DIR <- "D:/A/data"

OUT_ROOT <- file.path(DATA_DIR, "STEP5_PLEIOFDR", "00_fullGWAS_standardized")
INPUT_DIR <- file.path(OUT_ROOT, "01_inputs")
QC_DIR <- file.path(OUT_ROOT, "02_qc")

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# Packages
# =========================
pkgs <- c("DBI", "duckdb", "data.table")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(DBI)
library(duckdb)
library(data.table)

# =========================
# Fixed study metadata
# =========================

# AF GCST90624412 = European AF GWAS:
# 228,926 cases + 1,611,415 controls = 1,840,341 total.
AF_CASE <- 228926
AF_CTRL <- 1611415
AF_TOTAL <- AF_CASE + AF_CTRL
AF_CASE_FRAC <- AF_CASE / AF_TOTAL
AF_NEFF_FACTOR <- 4 * AF_CASE_FRAC * (1 - AF_CASE_FRAC)
AF_NEFF_TOTAL <- 4 / (1 / AF_CASE + 1 / AF_CTRL)

# FinnGen R9 OSA
OSA_CASE <- 38998
OSA_CTRL <- 336659
OSA_NEFF <- 4 / (1 / OSA_CASE + 1 / OSA_CTRL)

cat("AF effective-N factor =", AF_NEFF_FACTOR, "\n")
cat("AF overall effective N =", AF_NEFF_TOTAL, "\n")
cat("OSA effective N =", OSA_NEFF, "\n\n")

# =========================
# File discovery
# =========================
all_files <- list.files(
  DATA_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

choose_shortest <- function(x) {
  if (length(x) == 0) return(NA_character_)
  x[which.min(nchar(x))]
}

find_exact <- function(fname, prefer = NULL) {
  hit <- all_files[tolower(basename(all_files)) == tolower(fname)]

  if (!is.null(prefer) && length(hit) > 0) {
    h2 <- hit[grepl(prefer, hit, ignore.case = TRUE)]
    if (length(h2) > 0) return(choose_shortest(h2))
  }

  choose_shortest(hit)
}

AF_PATH <- find_exact("GCST90624412.tsv.gz")
HF_PATH <- find_exact(
  "FORMAT-METAL_Pheno4_EUR.tsv.gz",
  prefer = "STEP1B_HFpEF_PREP.*HERMES_HFpEF"
)
BMI_PATH <- find_exact("SNP_gwas_mc_merge_nogc.tbl.uniq.gz")
OSA_PATH <- find_exact("finngen_R9_G6_SLEEPAPNO.gz")

paths <- c(AF = AF_PATH, HFpEF = HF_PATH, BMI = BMI_PATH, OSA = OSA_PATH)

cat("Resolved raw files:\n")
print(paths)

if (any(is.na(paths)) || any(!file.exists(paths))) {
  stop(
    "One or more raw GWAS files could not be found.\n",
    paste(names(paths), paths, sep = " = ", collapse = "\n")
  )
}

# =========================
# DuckDB connection
# =========================
con <- dbConnect(
  duckdb::duckdb(),
  dbdir = file.path(OUT_ROOT, "STEP2A_duckdb_temp.db"),
  read_only = FALSE
)

on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

# Keep temporary work on D drive
dbExecute(
  con,
  sprintf(
    "SET temp_directory='%s';",
    gsub("'", "''", normalizePath(OUT_ROOT, winslash = "/", mustWork = FALSE))
  )
)

# Limit thread count modestly; edit if desired
dbExecute(con, "SET threads=6;")

# =========================
# SQL helpers
# =========================
sql_path <- function(x) {
  paste0("'", gsub("'", "''", normalizePath(x, winslash = "/", mustWork = TRUE)), "'")
}

qid <- function(x) {
  paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
}

get_schema <- function(path) {
  q <- sprintf(
    "DESCRIBE SELECT * FROM read_csv_auto(%s, header=true, sample_size=100000);",
    sql_path(path)
  )
  dbGetQuery(con, q)
}

pick_col <- function(cols, candidates, required = TRUE) {
  low <- tolower(cols)

  for (cand in candidates) {
    i <- which(low == tolower(cand))
    if (length(i) > 0) return(cols[i[1]])
  }

  if (required) {
    stop(
      "Required column not found. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(cols, collapse = " | ")
    )
  }

  NA_character_
}

is_available <- function(x) {
  length(x) == 1 && !is.na(x) && nzchar(x)
}

# INFO-like column: only use high-confidence names.
pick_info <- function(cols) {
  pick_col(
    cols,
    c(
      "info",
      "imputation_info",
      "imputationinfo",
      "mach_r2",
      "imputation_r2",
      "rsq"
    ),
    required = FALSE
  )
}

# =========================
# Build SQL for each trait
# =========================
make_trait_sql <- function(
  trait,
  path,
  snp_col,
  a1_col,
  a2_col,
  beta_col,
  se_col,
  p_col,
  eaf_col,
  n_expr,
  chr_col = NA_character_,
  bp_col = NA_character_,
  info_col = NA_character_,
  snp_expr_custom = NULL
) {

  src <- sprintf(
    "read_csv_auto(%s, header=true, sample_size=100000)",
    sql_path(path)
  )

  snp_expr <- if (!is.null(snp_expr_custom)) {
    snp_expr_custom
  } else {
    sprintf("CAST(%s AS VARCHAR)", qid(snp_col))
  }

  chr_expr <- if (is_available(chr_col)) {
    sprintf("TRY_CAST(%s AS INTEGER)", qid(chr_col))
  } else {
    "NULL::INTEGER"
  }

  bp_expr <- if (is_available(bp_col)) {
    sprintf("TRY_CAST(%s AS BIGINT)", qid(bp_col))
  } else {
    "NULL::BIGINT"
  }

  info_expr <- if (is_available(info_col)) {
    sprintf("TRY_CAST(%s AS DOUBLE)", qid(info_col))
  } else {
    "NULL::DOUBLE"
  }

  base_sql <- sprintf(
"
WITH x AS (
  SELECT
    %s AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    TRY_CAST(%s AS DOUBLE) AS EAF,
    %s AS N,
    %s AS CHR,
    %s AS BP,
    %s AS INFO
  FROM %s
),
q AS (
  SELECT *,
         BETA / SE AS Z
  FROM x
  WHERE
    SNP IS NOT NULL
    AND regexp_matches(SNP, '^rs[0-9]+$')
    AND A1 IN ('A','C','G','T')
    AND A2 IN ('A','C','G','T')
    AND A1 <> A2
    AND NOT (
      (A1='A' AND A2='T') OR
      (A1='T' AND A2='A') OR
      (A1='C' AND A2='G') OR
      (A1='G' AND A2='C')
    )
    AND BETA IS NOT NULL
    AND SE IS NOT NULL
    AND SE > 0
    AND P IS NOT NULL
    AND P > 0
    AND P <= 1
    AND EAF IS NOT NULL
    AND EAF > 0.01
    AND EAF < 0.99
    AND N IS NOT NULL
    AND N > 0
    %s
    %s
),
d AS (
  SELECT *,
         row_number() OVER (
           PARTITION BY SNP
           ORDER BY P ASC
         ) AS rn
  FROM q
)
SELECT SNP, A1, A2, N, Z
FROM d
WHERE rn = 1
",
    snp_expr,
    qid(a1_col),
    qid(a2_col),
    qid(beta_col),
    qid(se_col),
    qid(p_col),
    qid(eaf_col),
    n_expr,
    chr_expr,
    bp_expr,
    info_expr,
    src,
    if (is_available(info_col)) "AND INFO > 0.9" else "",
    if (is_available(chr_col) && is_available(bp_col))
      "AND NOT (CHR = 6 AND BP BETWEEN 25000000 AND 34000000)"
    else ""
  )

  base_sql
}

# =========================
# AF
# =========================
af_schema <- get_schema(AF_PATH)
af_cols <- af_schema$column_name

AF_SNP <- pick_col(af_cols, c("rs_id", "rsid", "SNP"))
AF_A1 <- pick_col(af_cols, c("effect_allele", "A1"))
AF_A2 <- pick_col(af_cols, c("other_allele", "A2"))
AF_BETA <- pick_col(af_cols, c("beta", "BETA"))
AF_SE <- pick_col(af_cols, c("standard_error", "se", "SE"))
AF_P <- pick_col(af_cols, c("p_value", "pval", "p"))
AF_EAF <- pick_col(af_cols, c("effect_allele_frequency", "eaf", "af"))
AF_N <- pick_col(af_cols, c("n", "N"))
AF_CHR <- pick_col(af_cols, c("chromosome", "chr", "#chrom"), required = FALSE)
AF_BP <- pick_col(af_cols, c("base_pair_location", "position", "pos", "bp"), required = FALSE)
AF_INFO <- pick_info(af_cols)

af_n_expr <- sprintf(
  "(TRY_CAST(%s AS DOUBLE) * %.15f)",
  qid(AF_N),
  AF_NEFF_FACTOR
)

AF_SQL <- make_trait_sql(
  trait = "AF",
  path = AF_PATH,
  snp_col = AF_SNP,
  a1_col = AF_A1,
  a2_col = AF_A2,
  beta_col = AF_BETA,
  se_col = AF_SE,
  p_col = AF_P,
  eaf_col = AF_EAF,
  n_expr = af_n_expr,
  chr_col = AF_CHR,
  bp_col = AF_BP,
  info_col = AF_INFO
)

# =========================
# HFpEF
# =========================
hf_schema <- get_schema(HF_PATH)
hf_cols <- hf_schema$column_name

HF_SNP <- pick_col(hf_cols, c("rsID", "rsid", "SNP"))
HF_A1 <- pick_col(hf_cols, c("A1"))
HF_A2 <- pick_col(hf_cols, c("A2"))
HF_BETA <- pick_col(hf_cols, c("A1_beta", "beta"))
HF_SE <- pick_col(hf_cols, c("se", "SE"))
HF_P <- pick_col(hf_cols, c("pval", "p"))
HF_EAF <- pick_col(hf_cols, c("A1_freq", "eaf"))
HF_NCASE <- pick_col(hf_cols, c("N_case", "n_case"))
HF_NTOTAL <- pick_col(hf_cols, c("N_total", "n_total"))
HF_CHR <- pick_col(hf_cols, c("chr", "chromosome"))
HF_BP <- pick_col(hf_cols, c("pos_b37", "bp", "position"))
HF_INFO <- pick_info(hf_cols)

hf_n_expr <- sprintf(
  "(4.0 / (1.0 / TRY_CAST(%s AS DOUBLE) + 1.0 / (TRY_CAST(%s AS DOUBLE) - TRY_CAST(%s AS DOUBLE))))",
  qid(HF_NCASE),
  qid(HF_NTOTAL),
  qid(HF_NCASE)
)

HF_SQL <- make_trait_sql(
  trait = "HFpEF",
  path = HF_PATH,
  snp_col = HF_SNP,
  a1_col = HF_A1,
  a2_col = HF_A2,
  beta_col = HF_BETA,
  se_col = HF_SE,
  p_col = HF_P,
  eaf_col = HF_EAF,
  n_expr = hf_n_expr,
  chr_col = HF_CHR,
  bp_col = HF_BP,
  info_col = HF_INFO
)

# =========================
# BMI (quantitative trait; use per-SNP N directly)
# =========================
bmi_schema <- get_schema(BMI_PATH)
bmi_cols <- bmi_schema$column_name

BMI_SNP <- pick_col(bmi_cols, c("SNP", "rsid"))
BMI_A1 <- pick_col(bmi_cols, c("A1"))
BMI_A2 <- pick_col(bmi_cols, c("A2"))
BMI_BETA <- pick_col(bmi_cols, c("b", "beta"))
BMI_SE <- pick_col(bmi_cols, c("se", "SE"))
BMI_P <- pick_col(bmi_cols, c("p", "pval"))
BMI_EAF <- pick_col(bmi_cols, c("Freq1.Hapmap", "Freq1", "eaf"))
BMI_N <- pick_col(bmi_cols, c("N", "n"))
BMI_CHR <- pick_col(bmi_cols, c("chr", "chromosome"), required = FALSE)
BMI_BP <- pick_col(bmi_cols, c("bp", "position", "pos"), required = FALSE)
BMI_INFO <- pick_info(bmi_cols)

bmi_n_expr <- sprintf("TRY_CAST(%s AS DOUBLE)", qid(BMI_N))

BMI_SQL <- make_trait_sql(
  trait = "BMI",
  path = BMI_PATH,
  snp_col = BMI_SNP,
  a1_col = BMI_A1,
  a2_col = BMI_A2,
  beta_col = BMI_BETA,
  se_col = BMI_SE,
  p_col = BMI_P,
  eaf_col = BMI_EAF,
  n_expr = bmi_n_expr,
  chr_col = BMI_CHR,
  bp_col = BMI_BP,
  info_col = BMI_INFO
)

# =========================
# OSA
# =========================
osa_schema <- get_schema(OSA_PATH)
osa_cols <- osa_schema$column_name

OSA_RSID <- pick_col(osa_cols, c("rsids", "rsid", "SNP"))
OSA_A2 <- pick_col(osa_cols, c("ref", "A2"))
OSA_A1 <- pick_col(osa_cols, c("alt", "A1"))
OSA_BETA <- pick_col(osa_cols, c("beta"))
OSA_SE <- pick_col(osa_cols, c("sebeta", "se"))
OSA_P <- pick_col(osa_cols, c("pval", "p"))
OSA_EAF <- pick_col(osa_cols, c("af_alt", "eaf"))
OSA_CHR <- pick_col(osa_cols, c("#chrom", "chrom", "chr"), required = FALSE)
OSA_BP <- pick_col(osa_cols, c("pos", "position", "bp"), required = FALSE)
OSA_INFO <- pick_info(osa_cols)

osa_snp_expr <- sprintf(
  "regexp_extract(CAST(%s AS VARCHAR), 'rs[0-9]+')",
  qid(OSA_RSID)
)

OSA_SQL <- make_trait_sql(
  trait = "OSA",
  path = OSA_PATH,
  snp_col = OSA_RSID,
  a1_col = OSA_A1,
  a2_col = OSA_A2,
  beta_col = OSA_BETA,
  se_col = OSA_SE,
  p_col = OSA_P,
  eaf_col = OSA_EAF,
  n_expr = sprintf("%.15f", OSA_NEFF),
  chr_col = OSA_CHR,
  bp_col = OSA_BP,
  info_col = OSA_INFO,
  snp_expr_custom = osa_snp_expr
)

TRAIT_SQL <- list(
  AF = AF_SQL,
  HFpEF = HF_SQL,
  BMI = BMI_SQL,
  OSA = OSA_SQL
)

TRAIT_PATH <- list(
  AF = AF_PATH,
  HFpEF = HF_PATH,
  BMI = BMI_PATH,
  OSA = OSA_PATH
)

INFO_USED <- c(
  AF = is_available(AF_INFO),
  HFpEF = is_available(HF_INFO),
  BMI = is_available(BMI_INFO),
  OSA = is_available(OSA_INFO)
)

MHC_INPUT_FILTER <- c(
  AF = is_available(AF_CHR) && is_available(AF_BP),
  HFpEF = TRUE,
  BMI = is_available(BMI_CHR) && is_available(BMI_BP),
  OSA = is_available(OSA_CHR) && is_available(OSA_BP)
)

# =========================
# Export full-GWAS inputs
# =========================
qc_list <- list()

for (trait in names(TRAIT_SQL)) {

  cat("\n====================================================\n")
  cat("Preparing full-GWAS input:", trait, "\n")
  cat("====================================================\n")

  raw_path <- TRAIT_PATH[[trait]]
  out_path <- file.path(INPUT_DIR, paste0(trait, "_fullGWAS.sumstats.gz"))

  if (file.exists(out_path)) {
    unlink(out_path)
  }

  # Raw row count
  raw_count <- dbGetQuery(
    con,
    sprintf(
      "SELECT count(*) AS n FROM read_csv_auto(%s, header=true, sample_size=100000);",
      sql_path(raw_path)
    )
  )$n[1]

  # Final count/stats (one scan)
  stats_sql <- sprintf(
    "
    SELECT
      count(*) AS final_rows,
      min(N) AS min_N,
      median(N) AS median_N,
      max(N) AS max_N,
      avg(Z*Z) AS mean_Z2,
      max(abs(Z)) AS max_abs_Z
    FROM (%s) s;
    ",
    TRAIT_SQL[[trait]]
  )

  st <- dbGetQuery(con, stats_sql)

  # Export gzip TSV
  copy_sql <- sprintf(
    "
    COPY (%s)
    TO '%s'
    (
      FORMAT CSV,
      DELIMITER '\t',
      HEADER TRUE,
      COMPRESSION GZIP
    );
    ",
    TRAIT_SQL[[trait]],
    gsub("'", "''", normalizePath(out_path, winslash = "/", mustWork = FALSE))
  )

  dbExecute(con, copy_sql)

  if (!file.exists(out_path) || file.info(out_path)$size == 0) {
    stop("Output was not created for ", trait, ": ", out_path)
  }

  qc_list[[trait]] <- data.frame(
    trait = trait,
    raw_file = raw_path,
    raw_rows = raw_count,
    final_rows = st$final_rows,
    retained_pct = 100 * st$final_rows / raw_count,
    min_N = st$min_N,
    median_N = st$median_N,
    max_N = st$max_N,
    mean_Z2 = st$mean_Z2,
    max_abs_Z = st$max_abs_Z,
    info_gt_0.9_applied = INFO_USED[trait],
    mhc_removed_in_input = MHC_INPUT_FILTER[trait],
    output_file = out_path,
    output_size_MB = file.info(out_path)$size / 1024^2,
    stringsAsFactors = FALSE
  )

  cat("Raw rows:   ", format(raw_count, big.mark = ","), "\n")
  cat("Final rows: ", format(st$final_rows, big.mark = ","), "\n")
  cat("Median N:   ", round(st$median_N, 1), "\n")
  cat("Output:     ", out_path, "\n")
}

qc <- rbindlist(qc_list, fill = TRUE)

fwrite(
  qc,
  file.path(QC_DIR, "STEP2A_pleioFDR_input_QC.csv")
)

# =========================
# Readiness / power gate
# =========================

# STEP1D h2 Z scores based on the completed LDSC results.
readiness <- data.table(
  trait = c("AF", "HFpEF", "BMI", "OSA"),
  LDSC_h2 = c(0.0479, 0.0013, 0.1123, 0.0322),
  LDSC_h2_SE = c(0.0050, 0.0005, 0.0061, 0.0025)
)

readiness[, h2_Z := LDSC_h2 / LDSC_h2_SE]

readiness[, preliminary_power := fifelse(
  h2_Z >= 4,
  "robust",
  fifelse(
    h2_Z >= 2,
    "borderline",
    "low"
  )
)]

# The definitive pleioFDR gate is NOT h2 Z:
# it will be univariate pleioFDR AIC/BIC + likelihood/residual fit.
readiness[, mixer_decision := fifelse(
  trait == "HFpEF",
  "RUN fit1/test1, but require acceptable AIC/BIC + fit before any bivariate HFpEF model",
  "RUN fit1/test1; proceed to bivariate only after acceptable univariate fit"
)]

fwrite(
  readiness,
  file.path(QC_DIR, "STEP2A_pleioFDR_readiness.csv")
)

# =========================
# Save column/schema audit
# =========================
schema_audit <- rbindlist(
  list(
    data.table(trait = "AF", column = af_cols),
    data.table(trait = "HFpEF", column = hf_cols),
    data.table(trait = "BMI", column = bmi_cols),
    data.table(trait = "OSA", column = osa_cols)
  )
)

fwrite(
  schema_audit,
  file.path(QC_DIR, "STEP2A_raw_header_audit.csv")
)

# =========================
# Final README
# =========================
readme <- c(
  paste0("STEP 2A completed: ", Sys.time()),
  "",
  "PURPOSE",
  "Prepare full-GWAS full-GWAS inputs from raw AF/HFpEF/BMI/OSA summary statistics.",
  "",
  "IMPORTANT DIFFERENCE FROM STEP1C/LDSC",
  "full-GWAS inputs are NOT restricted to HapMap3.",
  "Current pleioFDR documentation specifically recommends regenerating inputs without HapMap3 restriction for cross-trait pleioFDR.",
  "",
  "QC APPLIED",
  "- rsID required",
  "- biallelic A/C/G/T SNPs",
  "- strand ambiguous A/T and C/G variants removed",
  "- MAF > 0.01",
  "- INFO > 0.9 only when an INFO-like field is actually available",
  "- MHC chr6:25-34 Mb removed when CHR/BP are available in the source",
  "- duplicate rsIDs: retain smallest P",
  "",
  "EFFECTIVE N",
  paste0("AF: case-control effective-N factor = ", signif(AF_NEFF_FACTOR, 8),
         "; overall Neff ~ ", round(AF_NEFF_TOTAL)),
  "HFpEF: exact per-SNP Neff calculated from N_case and N_total.",
  paste0("OSA: constant Neff ~ ", round(OSA_NEFF),
         " from 38,998 cases and 336,659 controls."),
  "BMI: quantitative trait; use original per-SNP N.",
  "",
  "NEXT",
  "Do NOT run bivariate pleioFDR yet.",
  "First run univariate pleioFDR fit1 + test1 for AF, HFpEF, BMI, OSA.",
  "Inspect polygenicity, discoverability, AIC/BIC, likelihood profile and residual/QQ fit.",
  "HFpEF is a priori borderline from LDSC h2 Z ~2.6; if univariate pleioFDR is not supported by AIC/fit, stop HFpEF bivariate pleioFDR rather than tuning until it looks favorable.",
  "",
  "UPLOAD BACK",
  "1) STEP2A_pleioFDR_input_QC.csv",
  "2) STEP2A_pleioFDR_readiness.csv",
  "3) STEP2A_raw_header_audit.csv"
)

writeLines(
  readme,
  file.path(QC_DIR, "STEP2A_README.txt")
)

# Close DB explicitly
dbDisconnect(con, shutdown = TRUE)

cat("\n====================================================\n")
cat("STEP 2A COMPLETE\n")
cat("Output root:\n", OUT_ROOT, "\n")
cat("\nPlease upload ONLY the three small QC CSV files from:\n")
cat(QC_DIR, "\n")
cat("====================================================\n")
