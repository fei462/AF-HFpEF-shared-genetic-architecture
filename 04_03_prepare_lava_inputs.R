# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7B_V2_PREPARE_LAVA_INPUTS_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7B_V2 — PREPARE LAVA INPUTS (ROBUST / LOW-MEMORY)
# AF / HFpEF / BMI / OSA
#
# Fixes:
#   - Uses the CONFIRMED exact raw-GWAS column names.
#   - HFpEF effect column = A1_beta (not beta).
#   - Uses DuckDB so very large GWAS files are not retained in RAM.
#   - Uses TOTAL per-SNP N for LAVA (NOT MiXeR/pleioFDR Neff).
#   - BMI is continuous: cases=NA, controls=NA in input.info.
#
# Outputs:
#   01_sumstats/AF_LAVA.tsv.gz
#   01_sumstats/HFpEF_LAVA.tsv.gz
#   01_sumstats/BMI_LAVA.tsv.gz
#   01_sumstats/OSA_LAVA.tsv.gz
#   02_input/input.info.txt
#   02_input/STEP7B_sample_overlap_TEMPLATE.txt
#   07_tables/STEP7B_V2_sumstats_QC.csv
#   07_tables/STEP7B_V2_input_info_QC.csv
#   07_tables/STEP7B_V2_source_manifest.csv
#
# IMPORTANT:
#   - Run this from a CLEAN R session.
#   - Do NOT run local h2 yet.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")

SUM_DIR <- file.path(LAVA_ROOT, "01_sumstats")
IN_DIR  <- file.path(LAVA_ROOT, "02_input")
TAB_DIR <- file.path(LAVA_ROOT, "07_tables")
LOG_DIR <- file.path(LAVA_ROOT, "08_logs")
TMP_DIR <- file.path(LAVA_ROOT, "09_duckdb_tmp")

for (d in c(SUM_DIR, IN_DIR, TAB_DIR, LOG_DIR, TMP_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c("DBI", "duckdb", "data.table")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(DBI)
library(duckdb)
library(data.table)

# ============================================================
# 2. STUDY METADATA
# ============================================================

AF_CASES    <- 228926
AF_CONTROLS <- 1611415

OSA_CASES    <- 38998
OSA_CONTROLS <- 336659
OSA_TOTAL_N  <- OSA_CASES + OSA_CONTROLS   # 375657

# ============================================================
# 3. RESOLVE RAW GWAS FILES
# ============================================================

all_files <- list.files(
  DATA_ROOT,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

# Avoid accidentally selecting STEP7 outputs
all_files <- all_files[
  !grepl("STEP7_LAVA", all_files, fixed = TRUE)
]

choose_shortest <- function(x) {
  if (!length(x)) return(NA_character_)
  x[which.min(nchar(x))]
}

find_exact <- function(fname, prefer = NULL) {
  hit <- all_files[tolower(basename(all_files)) == tolower(fname)]

  if (!is.null(prefer) && length(hit)) {
    h2 <- hit[grepl(prefer, hit, ignore.case = TRUE)]
    if (length(h2)) return(choose_shortest(h2))
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

SOURCE <- c(
  AF = AF_PATH,
  HFpEF = HF_PATH,
  BMI = BMI_PATH,
  OSA = OSA_PATH
)

cat("\n================ SOURCE FILES ================\n")
print(SOURCE)

if (any(is.na(SOURCE)) || any(!file.exists(SOURCE))) {
  stop(
    "One or more raw GWAS files were not found:\n",
    paste(names(SOURCE), SOURCE, sep = " = ", collapse = "\n")
  )
}

# ============================================================
# 4. DUCKDB CONNECTION
# ============================================================

DB_FILE <- file.path(TMP_DIR, "STEP7B_V2_temp.duckdb")

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_FILE,
  read_only = FALSE
)

on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(
  con,
  sprintf(
    "SET temp_directory='%s';",
    gsub(
      "'",
      "''",
      normalizePath(TMP_DIR, winslash = "/", mustWork = FALSE)
    )
  )
)

dbExecute(con, "SET threads=6;")

sql_path <- function(x) {
  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(x, winslash = "/", mustWork = TRUE)
    ),
    "'"
  )
}

sql_out_path <- function(x) {
  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(x, winslash = "/", mustWork = FALSE)
    ),
    "'"
  )
}

qid <- function(x) {
  paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
}

get_schema <- function(path) {
  dbGetQuery(
    con,
    sprintf(
      "DESCRIBE SELECT * FROM read_csv_auto(%s, header=true, sample_size=100000);",
      sql_path(path)
    )
  )
}

assert_cols <- function(path, required_cols, trait) {
  schema <- get_schema(path)
  cols <- schema$column_name

  missing <- setdiff(required_cols, cols)

  if (length(missing)) {
    stop(
      "\n", trait, ": required columns missing:\n",
      paste(missing, collapse = ", "),
      "\n\nAvailable columns:\n",
      paste(cols, collapse = " | ")
    )
  }

  data.table(
    phenotype = trait,
    required_columns = paste(required_cols, collapse = " | "),
    status = "PASS"
  )
}

# ============================================================
# 5. HARD-CODED, CONFIRMED RAW COLUMN MAPS
# ============================================================

# These exact columns are taken from the project's earlier
# successful STEP1C / STEP2A processing.

AF_COLS <- c(
  "effect_allele",
  "other_allele",
  "beta",
  "standard_error",
  "p_value",
  "rs_id",
  "n"
)

HF_COLS <- c(
  "rsID",
  "A1",
  "A2",
  "A1_beta",
  "se",
  "pval",
  "N_total",
  "N_case"
)

BMI_COLS <- c(
  "SNP",
  "A1",
  "A2",
  "b",
  "se",
  "p",
  "N"
)

OSA_COLS <- c(
  "ref",
  "alt",
  "rsids",
  "pval",
  "beta",
  "sebeta"
)

source_qc <- rbindlist(list(
  assert_cols(AF_PATH, AF_COLS, "AF"),
  assert_cols(HF_PATH, HF_COLS, "HFpEF"),
  assert_cols(BMI_PATH, BMI_COLS, "BMI"),
  assert_cols(OSA_PATH, OSA_COLS, "OSA")
))

# ============================================================
# 6. BUILD LAVA QUERIES
#
# LAVA needs:
#   SNP  A1  A2  N  Z
#
# Z is always calculated from the confirmed signed effect / SE.
# ============================================================

AF_SQL <- sprintf(
"
WITH x AS (
  SELECT
    CAST(%s AS VARCHAR) AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    TRY_CAST(%s AS DOUBLE) AS N
  FROM read_csv_auto(%s, header=true, sample_size=100000)
),
q AS (
  SELECT
    SNP, A1, A2, N,
    BETA / SE AS Z,
    P
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
    AND SE IS NOT NULL AND SE > 0
    AND P IS NOT NULL AND P > 0 AND P <= 1
    AND N IS NOT NULL AND N > 0
),
d AS (
  SELECT *,
         row_number() OVER (
           PARTITION BY SNP
           ORDER BY P ASC, N DESC
         ) AS rn
  FROM q
)
SELECT SNP, A1, A2, N, Z
FROM d
WHERE rn = 1
  AND isfinite(Z)
  AND abs(Z) < 100
",
  qid("rs_id"),
  qid("effect_allele"),
  qid("other_allele"),
  qid("beta"),
  qid("standard_error"),
  qid("p_value"),
  qid("n"),
  sql_path(AF_PATH)
)

HF_SQL <- sprintf(
"
WITH x AS (
  SELECT
    CAST(%s AS VARCHAR) AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    TRY_CAST(%s AS DOUBLE) AS N
  FROM read_csv_auto(%s, header=true, sample_size=100000)
),
q AS (
  SELECT
    SNP, A1, A2, N,
    BETA / SE AS Z,
    P
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
    AND SE IS NOT NULL AND SE > 0
    AND P IS NOT NULL AND P > 0 AND P <= 1
    AND N IS NOT NULL AND N > 0
),
d AS (
  SELECT *,
         row_number() OVER (
           PARTITION BY SNP
           ORDER BY P ASC, N DESC
         ) AS rn
  FROM q
)
SELECT SNP, A1, A2, N, Z
FROM d
WHERE rn = 1
  AND isfinite(Z)
  AND abs(Z) < 100
",
  qid("rsID"),
  qid("A1"),
  qid("A2"),
  qid("A1_beta"),
  qid("se"),
  qid("pval"),
  qid("N_total"),
  sql_path(HF_PATH)
)

BMI_SQL <- sprintf(
"
WITH x AS (
  SELECT
    CAST(%s AS VARCHAR) AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    TRY_CAST(%s AS DOUBLE) AS N
  FROM read_csv_auto(%s, header=true, sample_size=100000)
),
q AS (
  SELECT
    SNP, A1, A2, N,
    BETA / SE AS Z,
    P
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
    AND SE IS NOT NULL AND SE > 0
    AND P IS NOT NULL AND P > 0 AND P <= 1
    AND N IS NOT NULL AND N > 0
),
d AS (
  SELECT *,
         row_number() OVER (
           PARTITION BY SNP
           ORDER BY P ASC, N DESC
         ) AS rn
  FROM q
)
SELECT SNP, A1, A2, N, Z
FROM d
WHERE rn = 1
  AND isfinite(Z)
  AND abs(Z) < 100
",
  qid("SNP"),
  qid("A1"),
  qid("A2"),
  qid("b"),
  qid("se"),
  qid("p"),
  qid("N"),
  sql_path(BMI_PATH)
)

OSA_SQL <- sprintf(
"
WITH x AS (
  SELECT
    regexp_extract(
      CAST(%s AS VARCHAR),
      'rs[0-9]+'
    ) AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    %.0f::DOUBLE AS N
  FROM read_csv_auto(%s, header=true, sample_size=100000)
),
q AS (
  SELECT
    SNP, A1, A2, N,
    BETA / SE AS Z,
    P
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
    AND SE IS NOT NULL AND SE > 0
    AND P IS NOT NULL AND P > 0 AND P <= 1
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
  AND isfinite(Z)
  AND abs(Z) < 100
",
  qid("rsids"),
  qid("alt"),   # FinnGen beta is for ALT
  qid("ref"),
  qid("beta"),
  qid("sebeta"),
  qid("pval"),
  OSA_TOTAL_N,
  sql_path(OSA_PATH)
)

TRAIT_SQL <- list(
  AF = AF_SQL,
  HFpEF = HF_SQL,
  BMI = BMI_SQL,
  OSA = OSA_SQL
)

# ============================================================
# 7. EXPORT EACH LAVA SUMSTATS FILE
# ============================================================

export_trait <- function(trait, sql) {

  outfile <- file.path(
    SUM_DIR,
    paste0(trait, "_LAVA.tsv.gz")
  )

  # Remove old output before overwriting
  if (file.exists(outfile)) file.remove(outfile)

  copy_sql <- sprintf(
    "COPY (%s) TO %s (HEADER, DELIMITER '\t', COMPRESSION GZIP);",
    sql,
    sql_out_path(outfile)
  )

  cat("\nExporting ", trait, " ...\n", sep = "")
  dbExecute(con, copy_sql)

  if (!file.exists(outfile)) {
    stop(trait, " output was not created.")
  }

  outfile
}

OUTPUTS <- mapply(
  export_trait,
  names(TRAIT_SQL),
  TRAIT_SQL,
  SIMPLIFY = TRUE,
  USE.NAMES = TRUE
)

cat("\n================ LAVA FILES CREATED ================\n")
print(OUTPUTS)

# ============================================================
# 8. QC THE EXPORTED FILES WITH DUCKDB
# ============================================================

qc_one <- function(trait, path) {

  src <- sprintf(
    "read_csv_auto(%s, header=true, sample_size=100000)",
    sql_path(path)
  )

  q <- sprintf(
"
SELECT
  count(*) AS lava_rows,
  median(N) AS median_N,
  min(N) AS min_N,
  max(N) AS max_N,
  avg(Z * Z) AS mean_Z2,
  max(abs(Z)) AS max_abs_Z,
  sum(CASE WHEN SNP IS NULL THEN 1 ELSE 0 END) AS missing_SNP,
  sum(CASE WHEN N IS NULL THEN 1 ELSE 0 END) AS missing_N,
  sum(CASE WHEN Z IS NULL THEN 1 ELSE 0 END) AS missing_Z
FROM %s
",
    src
  )

  z <- as.data.table(dbGetQuery(con, q))
  z[, phenotype := trait]
  z[, file := path]

  setcolorder(
    z,
    c(
      "phenotype", "lava_rows",
      "median_N", "min_N", "max_N",
      "mean_Z2", "max_abs_Z",
      "missing_SNP", "missing_N", "missing_Z",
      "file"
    )
  )

  z
}

QC <- rbindlist(
  lapply(names(OUTPUTS), function(tr) {
    qc_one(tr, OUTPUTS[[tr]])
  }),
  fill = TRUE
)

QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B_V2_sumstats_QC.csv"
)

fwrite(QC, QC_FILE)

cat("\n================ SUMSTATS QC ================\n")
print(QC)

# ============================================================
# 9. DERIVE HFpEF CASE/CONTROL METADATA
#
# input.info uses cases/controls to identify binary phenotypes
# and calculate sample case fraction.
# For HFpEF, use medians of per-SNP N_case and N_total-N_case.
# ============================================================

hf_meta_sql <- sprintf(
"
SELECT
  median(TRY_CAST(%s AS DOUBLE)) AS cases_median,
  median(
    TRY_CAST(%s AS DOUBLE) -
    TRY_CAST(%s AS DOUBLE)
  ) AS controls_median,
  median(TRY_CAST(%s AS DOUBLE)) AS total_median
FROM read_csv_auto(%s, header=true, sample_size=100000)
WHERE
  TRY_CAST(%s AS DOUBLE) > 0
  AND TRY_CAST(%s AS DOUBLE) > TRY_CAST(%s AS DOUBLE)
",
  qid("N_case"),
  qid("N_total"),
  qid("N_case"),
  qid("N_total"),
  sql_path(HF_PATH),
  qid("N_case"),
  qid("N_total"),
  qid("N_case")
)

HF_META <- dbGetQuery(con, hf_meta_sql)

HF_CASES <- round(HF_META$cases_median[1])
HF_CONTROLS <- round(HF_META$controls_median[1])

if (
  !is.finite(HF_CASES) ||
  !is.finite(HF_CONTROLS) ||
  HF_CASES <= 0 ||
  HF_CONTROLS <= 0
) {
  stop("Failed to derive valid HFpEF case/control metadata.")
}

# ============================================================
# 10. BUILD OFFICIAL LAVA input.info
#
# Official LAVA docs:
#   continuous phenotype -> cases=NA, controls=NA
# ============================================================

INFO <- data.table(
  phenotype = c("AF", "HFpEF", "BMI", "OSA"),
  cases = c(
    AF_CASES,
    HF_CASES,
    NA,
    OSA_CASES
  ),
  controls = c(
    AF_CONTROLS,
    HF_CONTROLS,
    NA,
    OSA_CONTROLS
  ),
  filename = unname(OUTPUTS[
    c("AF", "HFpEF", "BMI", "OSA")
  ])
)

INFO_FILE <- file.path(
  IN_DIR,
  "input.info.txt"
)

fwrite(
  INFO,
  INFO_FILE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

INFO_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B_V2_input_info_QC.csv"
)

fwrite(INFO, INFO_QC_FILE)

cat("\n================ input.info ================\n")
print(INFO)

# ============================================================
# 11. SAMPLE-OVERLAP TEMPLATE
#
# DO NOT invent zero sample overlap.
# Final matrix will be built from cross-trait LDSC gcov_int.
# ============================================================

phen <- c("AF", "HFpEF", "BMI", "OSA")

M <- matrix(
  NA_real_,
  nrow = length(phen),
  ncol = length(phen),
  dimnames = list(phen, phen)
)

diag(M) <- 1

OVERLAP_TEMPLATE <- file.path(
  IN_DIR,
  "STEP7B_sample_overlap_TEMPLATE.txt"
)

write.table(
  M,
  OVERLAP_TEMPLATE,
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

# ============================================================
# 12. REFERENCE PREFIX HARD CHECK
# ============================================================

REF_PREFIX <- file.path(
  LAVA_ROOT,
  "00_reference",
  "UKB_v1.1",
  "lava-ukb-v1.1"
)

ref_rows <- rbindlist(
  lapply(1:22, function(chr) {
    data.table(
      chr = chr,
      info_exists = file.exists(
        paste0(REF_PREFIX, "_chr", chr, ".info")
      ),
      bcor_exists = file.exists(
        paste0(REF_PREFIX, "_chr", chr, ".bcor")
      )
    )
  })
)

REF_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B_V2_reference_QC.csv"
)

fwrite(ref_rows, REF_QC_FILE)

if (
  any(!ref_rows$info_exists) ||
  any(!ref_rows$bcor_exists)
) {
  print(ref_rows[!info_exists | !bcor_exists])

  stop(
    "\nUKB v1.1 reference is still incomplete/misnamed.\n",
    "Fix the listed chromosomes before running LAVA."
  )
}

# ============================================================
# 13. LAVA PACKAGE CHECK
# ============================================================

lava_installed <- requireNamespace("LAVA", quietly = TRUE)

LAVA_QC <- data.table(
  installed = lava_installed,
  version = if (lava_installed) {
    as.character(packageVersion("LAVA"))
  } else {
    NA_character_
  }
)

fwrite(
  LAVA_QC,
  file.path(TAB_DIR, "STEP7B_V2_LAVA_package_QC.csv")
)

# ============================================================
# 14. SOURCE MANIFEST
# ============================================================

MANIFEST <- data.table(
  phenotype = c("AF", "HFpEF", "BMI", "OSA"),
  source_file = unname(SOURCE[c("AF", "HFpEF", "BMI", "OSA")]),
  effect_definition = c(
    "beta / standard_error; A1=effect_allele",
    "A1_beta / se; A1=A1",
    "b / se; A1=A1",
    "beta / sebeta; A1=ALT"
  ),
  N_definition = c(
    "per-SNP raw n (TOTAL sample size)",
    "per-SNP raw N_total",
    "per-SNP raw N",
    paste0(
      "constant TOTAL N=",
      OSA_TOTAL_N,
      " (",
      OSA_CASES,
      " cases + ",
      OSA_CONTROLS,
      " controls)"
    )
  )
)

MANIFEST_FILE <- file.path(
  TAB_DIR,
  "STEP7B_V2_source_manifest.csv"
)

fwrite(MANIFEST, MANIFEST_FILE)

# ============================================================
# 15. SANITY WINDOWS
# ============================================================

check_window <- function(tr, lo, hi) {
  x <- QC[phenotype == tr]$median_N

  if (
    length(x) != 1 ||
    !is.finite(x) ||
    x < lo ||
    x > hi
  ) {
    warning(
      tr,
      ": median N=",
      x,
      " is outside expected range [",
      lo,
      ", ",
      hi,
      "]"
    )
  }
}

check_window("AF",    1000000, 2000000)
check_window("HFpEF", 5000, 30000)
check_window("BMI",   100000, 500000)
check_window("OSA",   300000, 450000)

# ============================================================
# 16. SAVE SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(LOG_DIR, "STEP7B_V2_sessionInfo.txt")
)

# ============================================================
# 17. FINAL MESSAGE
# ============================================================

cat("\n====================================================\n")
cat("STEP7B_V2 COMPLETE\n")
cat("====================================================\n")

cat("\nQC files to upload to ChatGPT:\n")
cat("1) ", QC_FILE, "\n", sep = "")
cat("2) ", INFO_QC_FILE, "\n", sep = "")
cat("3) ", MANIFEST_FILE, "\n", sep = "")

cat("\nDO NOT run STEP7C yet.\n")
cat(
  "Next step after QC = build the LDSC-derived ",
  "sample.overlap.txt + test LAVA process.input().\n",
  sep = ""
)

cat("====================================================\n")
