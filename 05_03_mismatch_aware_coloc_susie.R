# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP8C_V3_1000G_EUR_MISMATCH_AWARE_COLOC_SUSIE(1).R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP8C_V3 — MISMATCH-AWARE MULTI-SIGNAL COLOCALIZATION WITH SuSiE + coloc.susie
#              USING 1000 GENOMES EUR PLINK GENOTYPES FOR LD
#
# Project: AF / HFpEF / BMI / OSA
#
# WHY V2
# ------
# The previous STEP8C used LAVA UKB .bcor LD. One regional matrix
# contained |r| > 1, which is unacceptable for SuSiE's correlation
# matrix input. V2 therefore DOES NOT use .bcor LD at all.
#
# Instead V2:
#   1) auto-detects the downloaded 1000G EUR PLINK .bed/.bim/.fam
#      reference under D:/A/data;
#   2) auto-detects PLINK 1.9;
#   3) extracts every target region from 1000G EUR;
#   4) aligns each GWAS signed Z to PLINK BIM A1/A2;
#   5) calculates signed LD directly from genotypes with
#         plink --r square bin4
#      while preserving A1/A2 order;
#   6) runs SuSiE with R_finite = actual 1000G EUR reference N
#      when the installed susieR supports it;
#   7) runs coloc.susie at p12=1e-5 and p12=1e-6.
#
# EIGHT PRE-SPECIFIED TESTS (unchanged)
# -------------------------------------
# chr16_locus2135_FTO:
#   AF-BMI, BMI-OSA, AF-OSA
# chr3_BMI_centered:
#   HFpEF-BMI, AF-BMI
# chr11_recurrent:
#   AF-OSA, AF-HFpEF
# chr2_AF_HFpEF:
#   AF-HFpEF
#
# HARD RULES
# ----------
# - No p-value SNP selection.
# - Reference MAF >= 0.01 only (1000G EUR is small, n~503).
# - Ambiguous A/T and C/G SNPs are excluded from signed alignment.
# - Non-ambiguous strand complements are allowed.
# - No nearPD / forced positive-definite repair.
# - Only tiny floating-point cleanup (<=1e-5) is permitted.
# - External reference LD => estimate_residual_variance = FALSE.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. USER SETTINGS
# ============================================================

DATA_ROOT <- "D:/A/data"

# Usually leave these as NA. The script auto-detects them.
# If auto-detection ever fails, set the exact values manually.
PLINK_EXE_MANUAL <- NA_character_
REF_ROOT_MANUAL  <- NA_character_

RESUME_IF_POSSIBLE <- TRUE
SAVE_SUSIE_FITS <- TRUE

SUSIE_L <- 10L

P1 <- 1e-4
P2 <- 1e-4
P12_DEFAULT <- 1e-5
P12_CONSERVATIVE <- 1e-6

# Official mismatch-aware SuSiE-RSS settings for finite external LD.
R_MISMATCH_MODE <- "eb"
SUSIE_FIRST_MAXIT <- 200L

REF_MAF_MIN <- 0.01
MIN_SNPS <- 100L

RUN_MISMATCH_DIAGNOSTIC <- TRUE
MISMATCH_MAX_SNPS <- 2500L

# ============================================================
# 1. PATHS
# ============================================================

LAVA_ROOT <- file.path(
  DATA_ROOT,
  "STEP7_LAVA"
)

COLOC_ROOT <- file.path(
  DATA_ROOT,
  "STEP8_COLOC"
)

REGION_FILE <- file.path(
  COLOC_ROOT,
  "03_tables",
  "STEP8A_target_regions.csv"
)

SUM_FILES <- c(
  AF = file.path(
    LAVA_ROOT,
    "01_sumstats",
    "AF_LAVA.tsv.gz"
  ),
  HFpEF = file.path(
    LAVA_ROOT,
    "01_sumstats",
    "HFpEF_LAVA.tsv.gz"
  ),
  BMI = file.path(
    LAVA_ROOT,
    "01_sumstats",
    "BMI_LAVA.tsv.gz"
  ),
  OSA = file.path(
    LAVA_ROOT,
    "01_sumstats",
    "OSA_LAVA.tsv.gz"
  )
)

INFO_FILE <- file.path(
  LAVA_ROOT,
  "02_input",
  "input.info.txt"
)

OUT_ROOT <- file.path(
  COLOC_ROOT,
  "07_coloc_susie_1000G_EUR_mismatch_aware"
)

WORK_DIR <- file.path(
  OUT_ROOT,
  "00_plink_work"
)

FIT_DIR <- file.path(
  OUT_ROOT,
  "01_susie_fits"
)

RESULT_DIR <- file.path(
  OUT_ROOT,
  "02_coloc_results"
)

TABLE_DIR <- file.path(
  OUT_ROOT,
  "03_tables"
)

LOG_DIR <- file.path(
  OUT_ROOT,
  "04_logs"
)

TMP_DIR <- file.path(
  OUT_ROOT,
  "99_duckdb_tmp"
)

for (d in c(
  OUT_ROOT,
  WORK_DIR,
  FIT_DIR,
  RESULT_DIR,
  TABLE_DIR,
  LOG_DIR,
  TMP_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

STATE_FILE <- file.path(
  OUT_ROOT,
  "STEP8C_V3_checkpoint_state.rds"
)

# ============================================================
# 2. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "DBI",
  "duckdb",
  "coloc",
  "susieR"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(data.table)
library(DBI)
library(duckdb)
library(coloc)
library(susieR)

COLOC_VERSION <- as.character(
  packageVersion("coloc")
)

SUSIER_VERSION <- as.character(
  packageVersion("susieR")
)

cat(
  "\nPackage versions:\n",
  "  coloc  = ", COLOC_VERSION, "\n",
  "  susieR = ", SUSIER_VERSION, "\n",
  sep = ""
)

# ============================================================
# 3. HARD FILE CHECK
# ============================================================

required_files <- c(
  REGION_FILE,
  INFO_FILE,
  SUM_FILES
)

if (any(!file.exists(required_files))) {
  stop(
    "Missing required file(s):\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    )
  )
}

# ============================================================
# 4. FROZEN TARGET REGIONS + 8 TESTS
# ============================================================

REGIONS <- fread(
  REGION_FILE
)

required_region_cols <- c(
  "region",
  "chr",
  "start",
  "stop"
)

if (!all(required_region_cols %in% names(REGIONS))) {
  stop(
    "STEP8A target-region file lacks: ",
    paste(
      setdiff(required_region_cols, names(REGIONS)),
      collapse = ", "
    )
  )
}

REGIONS[, region := as.character(region)]

TESTS <- data.table(
  test_id = c(
    "chr16_locus2135_FTO__AF__BMI",
    "chr16_locus2135_FTO__BMI__OSA",
    "chr16_locus2135_FTO__AF__OSA",
    "chr3_BMI_centered__HFpEF__BMI",
    "chr3_BMI_centered__AF__BMI",
    "chr11_recurrent__AF__OSA",
    "chr11_recurrent__AF__HFpEF",
    "chr2_AF_HFpEF__AF__HFpEF"
  ),
  region = c(
    "chr16_locus2135_FTO",
    "chr16_locus2135_FTO",
    "chr16_locus2135_FTO",
    "chr3_BMI_centered",
    "chr3_BMI_centered",
    "chr11_recurrent",
    "chr11_recurrent",
    "chr2_AF_HFpEF"
  ),
  trait1 = c(
    "AF",
    "BMI",
    "AF",
    "HFpEF",
    "AF",
    "AF",
    "AF",
    "AF"
  ),
  trait2 = c(
    "BMI",
    "OSA",
    "OSA",
    "BMI",
    "BMI",
    "OSA",
    "HFpEF",
    "HFpEF"
  ),
  STEP8B_class = c(
    "robust_H4",
    "robust_H4",
    "prior_sensitive_H4",
    "prior_sensitive_H4",
    "strong_H3",
    "strong_H3",
    "non_H4_mixed",
    "H3_H4_unresolved"
  )
)

TESTS <- merge(
  TESTS,
  REGIONS[
    ,
    .(
      region,
      chr,
      start,
      stop
    )
  ],
  by = "region",
  all.x = TRUE,
  sort = FALSE
)

test_order <- c(
  "chr16_locus2135_FTO__AF__BMI",
  "chr16_locus2135_FTO__BMI__OSA",
  "chr16_locus2135_FTO__AF__OSA",
  "chr3_BMI_centered__HFpEF__BMI",
  "chr3_BMI_centered__AF__BMI",
  "chr11_recurrent__AF__OSA",
  "chr11_recurrent__AF__HFpEF",
  "chr2_AF_HFpEF__AF__HFpEF"
)

TESTS[
  ,
  ord := match(
    test_id,
    test_order
  )
]

setorder(
  TESTS,
  ord
)

TESTS[, ord := NULL]

if (
  nrow(TESTS) != 8L ||
  anyNA(TESTS$chr) ||
  anyNA(TESTS$start) ||
  anyNA(TESTS$stop)
) {
  stop(
    "The 8 STEP8C tests do not map cleanly to frozen STEP8A regions."
  )
}

MANIFEST_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_test_manifest.csv"
)

fwrite(
  TESTS,
  MANIFEST_FILE
)

# ============================================================
# 5. TRAIT METADATA
# ============================================================

INFO <- fread(
  INFO_FILE,
  na.strings = c("NA", "")
)

TRAIT_TYPE <- c(
  AF = "cc",
  HFpEF = "cc",
  BMI = "quant",
  OSA = "cc"
)

CASE_PROP <- setNames(
  rep(NA_real_, length(TRAIT_TYPE)),
  names(TRAIT_TYPE)
)

for (tr in names(TRAIT_TYPE)) {

  if (TRAIT_TYPE[[tr]] == "cc") {

    ca <- as.numeric(
      INFO[
        phenotype == tr,
        cases
      ]
    )

    co <- as.numeric(
      INFO[
        phenotype == tr,
        controls
      ]
    )

    if (
      length(ca) != 1L ||
      length(co) != 1L ||
      !is.finite(ca) ||
      !is.finite(co) ||
      ca <= 0 ||
      co <= 0
    ) {
      stop(
        "Invalid case/control metadata for ",
        tr
      )
    }

    CASE_PROP[[tr]] <- ca / (
      ca + co
    )
  }
}

# ============================================================
# 6. AUTO-DETECT PLINK 1.9
# ============================================================

find_plink <- function() {

  if (
    !is.na(PLINK_EXE_MANUAL) &&
    nzchar(PLINK_EXE_MANUAL) &&
    file.exists(PLINK_EXE_MANUAL)
  ) {
    return(
      normalizePath(
        PLINK_EXE_MANUAL,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  path_hits <- c(
    Sys.which("plink"),
    Sys.which("plink.exe")
  )

  path_hits <- path_hits[
    nzchar(path_hits) &
    file.exists(path_hits)
  ]

  if (length(path_hits)) {
    return(
      normalizePath(
        path_hits[1],
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  roots <- c(
    DATA_ROOT,
    "D:/A",
    "E:/app"
  )

  hits <- character(0)

  for (rr in roots[dir.exists(roots)]) {

    h <- tryCatch(
      list.files(
        rr,
        pattern = "^plink(1\\.9|_1\\.9)?(\\.exe)?$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      ),
      error = function(e) character(0)
    )

    hits <- unique(
      c(
        hits,
        h
      )
    )

    if (length(hits)) break
  }

  if (!length(hits)) {
    stop(
      "\nPLINK 1.9 executable was not found.\n",
      "Place plink.exe somewhere under D:/A/data, or set:\n",
      "PLINK_EXE_MANUAL <- 'D:/.../plink.exe'\n",
      "at the top of this script."
    )
  }

  normalizePath(
    hits[1],
    winslash = "/",
    mustWork = TRUE
  )
}

PLINK_EXE <- find_plink()

cat(
  "\nPLINK executable:\n",
  PLINK_EXE,
  "\n",
  sep = ""
)

plink_version_log <- tryCatch(
  system2(
    PLINK_EXE,
    args = "--version",
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) conditionMessage(e)
)

cat(
  paste(
    plink_version_log,
    collapse = "\n"
  ),
  "\n"
)

if (
  any(
    grepl(
      "PLINK v2",
      plink_version_log,
      ignore.case = TRUE
    )
  )
) {
  stop(
    "This script requires PLINK 1.9, not PLINK 2.x, ",
    "because the binary signed-LD command/output is fixed to PLINK 1.9."
  )
}

# ============================================================
# 7. AUTO-DETECT 1000G EUR PLINK BFILES
# ============================================================

count_fam <- function(famfile) {

  con <- file(
    famfile,
    open = "r"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  n <- 0L

  repeat {

    z <- readLines(
      con,
      n = 10000L,
      warn = FALSE
    )

    if (!length(z)) break

    n <- n + length(z)
  }

  n
}

find_bfile_prefixes <- function() {

  roots <- if (
    !is.na(REF_ROOT_MANUAL) &&
    nzchar(REF_ROOT_MANUAL) &&
    dir.exists(REF_ROOT_MANUAL)
  ) {
    REF_ROOT_MANUAL
  } else {
    DATA_ROOT
  }

  beds <- list.files(
    roots,
    pattern = "\\.bed$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(beds)) {
    stop(
      "\nNo .bed files found under: ",
      roots,
      "\nThe downloaded 1000G EUR archive must be EXTRACTED first."
    )
  }

  prefixes <- sub(
    "\\.bed$",
    "",
    beds,
    ignore.case = TRUE
  )

  keep <- (
    file.exists(
      paste0(prefixes, ".bim")
    ) &
    file.exists(
      paste0(prefixes, ".fam")
    )
  )

  prefixes <- unique(
    prefixes[keep]
  )

  if (!length(prefixes)) {
    stop(
      "No complete PLINK .bed/.bim/.fam trio was found."
    )
  }

  meta <- rbindlist(
    lapply(
      prefixes,
      function(pr) {

        fam <- paste0(
          pr,
          ".fam"
        )

        n_fam <- tryCatch(
          count_fam(fam),
          error = function(e) NA_integer_
        )

        data.table(
          prefix = normalizePath(
            pr,
            winslash = "/",
            mustWork = FALSE
          ),
          n_fam = n_fam,
          name_has_EUR = grepl(
            "EUR|europe",
            pr,
            ignore.case = TRUE
          )
        )
      }
    )
  )

  # Official 1000G EUR has ~503 samples.
  eur <- meta[
    is.finite(n_fam) &
    n_fam >= 450 &
    n_fam <= 550
  ]

  if (!nrow(eur)) {

    fwrite(
      meta,
      file.path(
        TABLE_DIR,
        "STEP8C_V3_detected_PLINK_candidates_ALL.csv"
      )
    )

    stop(
      "\nPLINK bfiles were found, but none had ~503 EUR samples.\n",
      "Inspect STEP8C_V3_detected_PLINK_candidates_ALL.csv.\n",
      "If your EUR reference is in another folder, set REF_ROOT_MANUAL."
    )
  }

  eur[
    ,
    score := (
      100 * as.integer(name_has_EUR) -
      abs(n_fam - 503)
    )
  ]

  setorder(
    eur,
    -score
  )

  eur
}

BFILE_CANDIDATES <- find_bfile_prefixes()

BFILE_CANDIDATE_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_1000G_EUR_bfile_candidates.csv"
)

fwrite(
  BFILE_CANDIDATES,
  BFILE_CANDIDATE_FILE
)

cat(
  "\n1000G EUR candidate PLINK prefixes:\n"
)

print(
  BFILE_CANDIDATES[
    ,
    .(
      prefix,
      n_fam,
      score
    )
  ]
)

# ============================================================
# 8. BIM READER + CHROMOSOME PREFIX SELECTION
# ============================================================

BIM_CACHE <- new.env(
  parent = emptyenv()
)

read_bim <- function(prefix) {

  key <- gsub(
    "[^A-Za-z0-9]",
    "_",
    prefix
  )

  if (
    exists(
      key,
      envir = BIM_CACHE,
      inherits = FALSE
    )
  ) {
    return(
      get(
        key,
        envir = BIM_CACHE,
        inherits = FALSE
      )
    )
  }

  x <- fread(
    paste0(
      prefix,
      ".bim"
    ),
    header = FALSE,
    col.names = c(
      "CHR",
      "SNP",
      "CM",
      "BP",
      "A1",
      "A2"
    ),
    colClasses = list(
      character = c(
        "CHR",
        "SNP",
        "A1",
        "A2"
      ),
      numeric = "CM",
      integer = "BP"
    ),
    showProgress = FALSE
  )

  x[
    ,
    `:=`(
      SNP = tolower(
        as.character(SNP)
      ),
      A1 = toupper(
        as.character(A1)
      ),
      A2 = toupper(
        as.character(A2)
      )
    )
  ]

  assign(
    key,
    x,
    envir = BIM_CACHE
  )

  x
}

choose_prefix_for_region <- function(
  chr_i,
  start_i,
  stop_i
) {

  candidates <- list()

  for (
    k in seq_len(
      nrow(BFILE_CANDIDATES)
    )
  ) {

    pr <- BFILE_CANDIDATES$prefix[k]

    bim <- tryCatch(
      read_bim(pr),
      error = function(e) NULL
    )

    if (is.null(bim)) next

    region_n <- bim[
      as.integer(CHR) == chr_i &
      BP >= start_i &
      BP <= stop_i,
      .N
    ]

    if (region_n <= 0) next

    n_chr <- uniqueN(
      bim$CHR
    )

    name_chr_bonus <- as.integer(
      grepl(
        paste0(
          "(chr)?",
          chr_i,
          "([^0-9]|$)"
        ),
        basename(pr),
        ignore.case = TRUE,
        perl = TRUE
      )
    )

    candidates[[
      length(candidates) + 1L
    ]] <- data.table(
      prefix = pr,
      n_fam = BFILE_CANDIDATES$n_fam[k],
      region_snps = region_n,
      n_chromosomes = n_chr,
      score = (
        BFILE_CANDIDATES$score[k] +
        20 * name_chr_bonus -
        0.1 * n_chr
      )
    )
  }

  if (!length(candidates)) {
    stop(
      "No detected 1000G EUR PLINK bfile contains chr",
      chr_i,
      ":",
      start_i,
      "-",
      stop_i,
      ". Check reference genome build / extracted files."
    )
  }

  z <- rbindlist(
    candidates
  )

  setorder(
    z,
    -score,
    -region_snps
  )

  z[1]
}

REGION_BFILE <- rbindlist(
  lapply(
    seq_len(
      nrow(REGIONS)
    ),
    function(i) {

      r <- REGIONS[i]

      z <- choose_prefix_for_region(
        chr_i = as.integer(r$chr),
        start_i = as.integer(r$start),
        stop_i = as.integer(r$stop)
      )

      data.table(
        region = r$region,
        chr = r$chr,
        start = r$start,
        stop = r$stop,
        prefix = z$prefix,
        reference_N = z$n_fam,
        reference_region_snps = z$region_snps
      )
    }
  )
)

REGION_BFILE_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_region_reference_bfiles.csv"
)

fwrite(
  REGION_BFILE,
  REGION_BFILE_FILE
)

cat(
  "\nReference bfile selected per region:\n"
)

print(
  REGION_BFILE
)

# ============================================================
# 9. BUILD MASTER REFERENCE TARGET TABLE FROM BIM
# ============================================================

valid_base <- function(x) {
  x %in% c(
    "A",
    "C",
    "G",
    "T"
  )
}

is_ambiguous <- function(a1, a2) {
  paste0(
    a1,
    a2
  ) %in% c(
    "AT",
    "TA",
    "CG",
    "GC"
  )
}

REF_TARGET_LIST <- list()

for (
  i in seq_len(
    nrow(REGION_BFILE)
  )
) {

  r <- REGION_BFILE[i]

  bim <- read_bim(
    r$prefix
  )

  x <- bim[
    as.integer(CHR) == as.integer(r$chr) &
    BP >= as.integer(r$start) &
    BP <= as.integer(r$stop) &
    valid_base(A1) &
    valid_base(A2) &
    A1 != A2 &
    !is_ambiguous(A1, A2),
    .(
      region = as.character(r$region),
      CHR = as.integer(CHR),
      SNP = as.character(SNP),
      BP = as.integer(BP),
      REF_A1 = A1,
      REF_A2 = A2,
      bfile_prefix = as.character(r$prefix),
      reference_N = as.integer(r$reference_N)
    )
  ]

  # PLINK --extract is ID-based; duplicated rsIDs are unsafe.
  x <- x[
    !duplicated(SNP) &
    !duplicated(
      SNP,
      fromLast = TRUE
    )
  ]

  if (nrow(x) < MIN_SNPS) {
    stop(
      "Too few usable non-ambiguous 1000G EUR SNPs in ",
      r$region,
      ": ",
      nrow(x)
    )
  }

  REF_TARGET_LIST[[
    length(REF_TARGET_LIST) + 1L
  ]] <- x
}

REF_TARGETS <- rbindlist(
  REF_TARGET_LIST,
  fill = TRUE
)

REF_TARGET_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_reference_target_SNPs.csv.gz"
)

fwrite(
  REF_TARGETS,
  REF_TARGET_FILE
)

# ============================================================
# 10. EXTRACT GWAS A1/A2/N/Z FOR TARGET REFERENCE SNPs
#     USING DUCKDB — ONLY 4 LARGE FILE SCANS
# ============================================================

DB_FILE <- file.path(
  TMP_DIR,
  "STEP8C_V3_temp.duckdb"
)

if (file.exists(DB_FILE)) {
  try(
    unlink(
      DB_FILE,
      force = TRUE
    ),
    silent = TRUE
  )
}

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_FILE,
  read_only = FALSE
)

dbExecute(
  con,
  sprintf(
    "SET temp_directory='%s';",
    gsub(
      "'",
      "''",
      normalizePath(
        TMP_DIR,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
)

dbExecute(
  con,
  "SET threads=6;"
)

dbWriteTable(
  con,
  "ref_targets",
  as.data.frame(
    REF_TARGETS
  ),
  overwrite = TRUE
)

sql_path <- function(x) {

  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(
        x,
        winslash = "/",
        mustWork = TRUE
      )
    ),
    "'"
  )
}

GWAS_REGION <- list()

for (tr in names(SUM_FILES)) {

  cat(
    "\nExtracting A1/A2/Z for ",
    tr,
    " target regions ...\n",
    sep = ""
  )

  q <- sprintf(
"
SELECT
  r.region,
  r.CHR,
  r.BP,
  r.REF_A1,
  r.REF_A2,
  r.bfile_prefix,
  r.reference_N,
  lower(CAST(s.SNP AS VARCHAR)) AS SNP,
  upper(CAST(s.A1 AS VARCHAR)) AS GWAS_A1,
  upper(CAST(s.A2 AS VARCHAR)) AS GWAS_A2,
  TRY_CAST(s.N AS DOUBLE) AS N,
  TRY_CAST(s.Z AS DOUBLE) AS Z
FROM read_csv_auto(
  %s,
  header=true,
  sample_size=100000
) AS s
INNER JOIN ref_targets AS r
  ON lower(CAST(s.SNP AS VARCHAR)) = r.SNP
WHERE
  TRY_CAST(s.N AS DOUBLE) > 0
  AND TRY_CAST(s.Z AS DOUBLE) IS NOT NULL
",
    sql_path(
      SUM_FILES[[tr]]
    )
  )

  x <- as.data.table(
    dbGetQuery(
      con,
      q
    )
  )

  x[
    ,
    phenotype := tr
  ]

  x <- unique(
    x,
    by = c(
      "region",
      "SNP"
    )
  )

  GWAS_REGION[[tr]] <- x

  cat(
    "  extracted rows = ",
    nrow(x),
    "\n",
    sep = ""
  )
}

try(
  dbDisconnect(
    con,
    shutdown = TRUE
  ),
  silent = TRUE
)

rm(con)
gc()

# ============================================================
# 11. ALLELE ALIGNMENT TO PLINK BIM A1
# ============================================================

comp_allele <- function(x) {

  y <- rep(
    NA_character_,
    length(x)
  )

  y[x == "A"] <- "T"
  y[x == "T"] <- "A"
  y[x == "C"] <- "G"
  y[x == "G"] <- "C"

  y
}

align_to_ref <- function(x) {

  g1 <- toupper(
    x$GWAS_A1
  )

  g2 <- toupper(
    x$GWAS_A2
  )

  r1 <- toupper(
    x$REF_A1
  )

  r2 <- toupper(
    x$REF_A2
  )

  c1 <- comp_allele(
    g1
  )

  c2 <- comp_allele(
    g2
  )

  orientation <- rep(
    NA_integer_,
    nrow(x)
  )

  # Same strand
  orientation[
    g1 == r1 &
    g2 == r2
  ] <- 1L

  orientation[
    g1 == r2 &
    g2 == r1
  ] <- -1L

  # Reverse complement; safe because ambiguous SNPs were removed.
  orientation[
    is.na(orientation) &
    c1 == r1 &
    c2 == r2
  ] <- 1L

  orientation[
    is.na(orientation) &
    c1 == r2 &
    c2 == r1
  ] <- -1L

  x[
    ,
    orientation := orientation
  ]

  x[
    ,
    Z_ref := Z * orientation
  ]

  x
}

for (tr in names(GWAS_REGION)) {
  GWAS_REGION[[tr]] <- align_to_ref(
    GWAS_REGION[[tr]]
  )
}

ALIGN_QC <- rbindlist(
  lapply(
    names(GWAS_REGION),
    function(tr) {

      x <- GWAS_REGION[[tr]]

      x[
        ,
        .(
          total_reference_matched = .N,
          aligned = sum(
            !is.na(orientation)
          ),
          incompatible = sum(
            is.na(orientation)
          ),
          flipped = sum(
            orientation == -1L,
            na.rm = TRUE
          ),
          complement_matches = sum(
            !is.na(orientation) &
            (
              comp_allele(GWAS_A1) == REF_A1 |
              comp_allele(GWAS_A1) == REF_A2
            )
          )
        ),
        by = region
      ][
        ,
        phenotype := tr
      ]
    }
  ),
  fill = TRUE
)

ALIGN_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_allele_alignment_QC.csv"
)

fwrite(
  ALIGN_QC,
  ALIGN_QC_FILE
)

# ============================================================
# 12. PLINK HELPERS
# ============================================================

run_plink <- function(
  args,
  logfile
) {

  dir.create(
    dirname(logfile),
    recursive = TRUE,
    showWarnings = FALSE
  )

  status <- system2(
    PLINK_EXE,
    args = args,
    stdout = logfile,
    stderr = logfile
  )

  if (
    !identical(
      as.integer(status),
      0L
    )
  ) {

    tail_log <- tryCatch(
      tail(
        readLines(
          logfile,
          warn = FALSE
        ),
        40
      ),
      error = function(e) character(0)
    )

    stop(
      "PLINK failed with status ",
      status,
      ". Last log lines:\n",
      paste(
        tail_log,
        collapse = "\n"
      )
    )
  }

  invisible(
    TRUE
  )
}

quote_arg <- function(x) {
  shQuote(
    normalizePath(
      x,
      winslash = "/",
      mustWork = FALSE
    )
  )
}

read_plink_bin4_ld <- function(
  ld_file,
  n_snps
) {

  con <- file(
    ld_file,
    open = "rb"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  vals <- readBin(
    con,
    what = "numeric",
    n = n_snps * n_snps,
    size = 4L,
    signed = TRUE,
    endian = "little"
  )

  if (
    length(vals) !=
    n_snps * n_snps
  ) {
    stop(
      "PLINK LD binary size mismatch: expected ",
      n_snps * n_snps,
      " values, got ",
      length(vals)
    )
  }

  matrix(
    vals,
    nrow = n_snps,
    ncol = n_snps,
    byrow = TRUE
  )
}

# ============================================================
# 13. SuSiE CAPABILITY DETECTION
# ============================================================

SUSIE_RSS_FORMALS <- names(
  formals(
    susieR::susie_rss
  )
)

SUPPORTS_R_FINITE <- (
  "R_finite" %in%
    SUSIE_RSS_FORMALS
)

SUPPORTS_R_MISMATCH <- (
  "R_mismatch" %in%
    SUSIE_RSS_FORMALS
)

SUPPORTS_Z_METHOD <- (
  "z_method" %in%
    SUSIE_RSS_FORMALS
)

cat(
  "\nSuSiE capabilities:\n",
  "  R_finite   = ",
  SUPPORTS_R_FINITE,
  "\n",
  "  R_mismatch = ",
  SUPPORTS_R_MISMATCH,
  "\n",
  "  z_method   = ",
  SUPPORTS_Z_METHOD,
  "\n",
  sep = ""
)

if (!SUPPORTS_R_FINITE || !SUPPORTS_R_MISMATCH) {
  stop(
    "\nYour installed susieR is too old for STEP8C_V3.\n",
    "Please update susieR, restart R, and verify:\n",
    "'R_finite' %in% names(formals(susieR::susie_rss))\n",
    "'R_mismatch' %in% names(formals(susieR::susie_rss))\n",
    "Both must be TRUE before running this analysis."
  )
}

# ============================================================
# 14. GENERIC SAFE CALL
# ============================================================

safe_call <- function(expr) {

  warnings_vec <- character(0)
  messages_vec <- character(0)
  err <- NULL

  value <- tryCatch(
    withCallingHandlers(
      withCallingHandlers(
        eval.parent(
          substitute(expr)
        ),
        message = function(m) {
          messages_vec <<- c(
            messages_vec,
            conditionMessage(m)
          )
          invokeRestart(
            "muffleMessage"
          )
        }
      ),
      warning = function(w) {
        warnings_vec <<- c(
          warnings_vec,
          conditionMessage(w)
        )
        invokeRestart(
          "muffleWarning"
        )
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(
    ok = is.null(err),
    value = value,
    error = err,
    warnings = unique(warnings_vec),
    messages = unique(messages_vec)
  )
}

compact_messages <- function(
  x,
  max_n = 10L
) {

  x <- trimws(
    as.character(x)
  )

  x <- unique(
    x[nzchar(x)]
  )

  if (!length(x)) {
    return("")
  }

  if (length(x) > max_n) {
    x <- c(
      x[seq_len(max_n)],
      paste0(
        "... [",
        length(x) - max_n,
        " more]"
      )
    )
  }

  paste(
    x,
    collapse = " || "
  )
}

safe_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_\\-]+",
    "_",
    x
  )
}

# ============================================================
# 15. BUILD coloc/SuSiE DATASET
#
# runsusie() requires beta/varbeta for validation, but the SuSiE
# RSS model is driven by signed z + LD when z is supplied.
# Standardized proxy values keep beta/SE exactly equal to Z.
# ============================================================

build_susie_dataset <- function(
  tr,
  ss,
  R
) {

  n_med <- round(
    median(
      ss$N,
      na.rm = TRUE
    )
  )

  if (
    !is.finite(n_med) ||
    n_med <= 1
  ) {
    stop(
      "Invalid median N for ",
      tr
    )
  }

  z <- as.numeric(
    ss$Z_ref
  )

  se_std <- rep(
    1 / sqrt(n_med),
    length(z)
  )

  d <- list(
    beta = z * se_std,
    varbeta = se_std^2,
    z = z,
    snp = as.character(
      ss$SNP
    ),
    N = n_med,
    type = TRAIT_TYPE[[tr]],
    LD = R
  )

  if (
    TRAIT_TYPE[[tr]] == "cc"
  ) {
    d$s <- CASE_PROP[[tr]]
  } else {
    # The beta/varbeta pair above is deliberately on a standardized
    # phenotype scale only to satisfy coloc::runsusie() validation.
    # SuSiE itself uses the supplied signed z vector. Therefore sdY=1
    # is the internally consistent quantitative-trait scale here.
    d$sdY <- 1
  }

  rownames(
    d$LD
  ) <- d$snp

  colnames(
    d$LD
  ) <- d$snp

  coloc::check_dataset(
    d,
    req = c(
      "beta",
      "varbeta",
      "LD",
      "snp",
      "N"
    ),
    warn.minp = 1
  )

  d
}

run_susie_safe <- function(
  d,
  label,
  ref_n
) {

  if (
    !is.finite(ref_n) ||
    ref_n <= 1
  ) {
    stop(
      "Invalid external LD reference sample size for ",
      label
    )
  }

  args <- list(
    d = d,
    suffix = label,
    maxit = SUSIE_FIRST_MAXIT,
    repeat_until_convergence = TRUE,
    L = SUSIE_L,
    estimate_residual_variance = FALSE,

    # Core V3 correction:
    # account for finite reference sampling uncertainty plus
    # broader summary-statistics / LD population mismatch.
    R_finite = ref_n,
    R_mismatch = R_MISMATCH_MODE,

    # Keep diagnostics in the returned fit.
    track_fit = TRUE
  )

  if (
    SUPPORTS_Z_METHOD
  ) {
    args$z_method <- "wald"
  }

  safe_call(
    do.call(
      coloc::runsusie,
      args
    )
  )
}

# Extract official mismatch-aware diagnostics from susie_rss output.
diag_scalar <- function(
  fit,
  name
) {

  d <- fit$R_finite_diagnostics

  if (
    is.null(d) ||
    is.null(d[[name]])
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    as.numeric(
      d[[name]][1]
    )
  )
}

diag_flag <- function(
  fit,
  name
) {

  d <- fit$R_finite_diagnostics

  if (
    is.null(d) ||
    is.null(d[[name]])
  ) {
    return(
      NA
    )
  }

  as.logical(
    d[[name]][1]
  )
}

# ============================================================
# 16. LD / SUMMARY MISMATCH DIAGNOSTIC
# ============================================================

estimate_s_safe <- function(
  z,
  R,
  n
) {

  if (
    !RUN_MISMATCH_DIAGNOSTIC
  ) {
    return(
      list(
        value = NA_real_,
        status = "disabled"
      )
    )
  }

  if (
    length(z) >
    MISMATCH_MAX_SNPS
  ) {
    return(
      list(
        value = NA_real_,
        status = paste0(
          "skipped_nsnps>",
          MISMATCH_MAX_SNPS
        )
      )
    )
  }

  if (
    !"estimate_s_rss" %in%
    getNamespaceExports(
      "susieR"
    )
  ) {
    return(
      list(
        value = NA_real_,
        status = "function_unavailable"
      )
    )
  }

  ans <- tryCatch(
    susieR::estimate_s_rss(
      z = z,
      R = R,
      n = n,
      method = "null-mle"
    ),
    error = function(e) e
  )

  if (
    inherits(
      ans,
      "error"
    )
  ) {
    return(
      list(
        value = NA_real_,
        status = paste0(
          "error:",
          conditionMessage(ans)
        )
      )
    )
  }

  list(
    value = as.numeric(ans),
    status = "PASS"
  )
}

# ============================================================
# 17. FINEMAPPING HELPERS
# ============================================================

extract_pip <- function(
  fit,
  test_id,
  trait
) {

  pip <- as.numeric(
    fit$pip
  )

  snps <- names(
    fit$pip
  )

  if (
    is.null(snps) ||
    length(snps) != length(pip)
  ) {
    snps <- colnames(
      fit$alpha
    )
  }

  data.table(
    test_id = test_id,
    trait = trait,
    SNP = snps,
    PIP = pip
  )[
    order(
      -PIP
    )
  ]
}

extract_cs <- function(
  fit,
  test_id,
  trait
) {

  cs_list <- fit$sets$cs

  if (
    is.null(cs_list) ||
    !length(cs_list)
  ) {
    return(
      data.table()
    )
  }

  snp_order <- colnames(
    fit$alpha
  )

  pip_dt <- extract_pip(
    fit,
    test_id,
    trait
  )

  out <- list()

  for (
    k in seq_along(
      cs_list
    )
  ) {

    idx <- as.integer(
      cs_list[[k]]
    )

    idx <- idx[
      idx >= 1 &
      idx <= length(snp_order)
    ]

    cs_snps <- snp_order[
      idx
    ]

    pp <- pip_dt[
      SNP %in% cs_snps
    ]

    setorder(
      pp,
      -PIP
    )

    out[[
      length(out) + 1L
    ]] <- data.table(
      test_id = test_id,
      trait = trait,
      cs = k,
      cs_size = length(cs_snps),
      lead_snp = if (
        nrow(pp)
      ) pp$SNP[1] else NA_character_,
      lead_PIP = if (
        nrow(pp)
      ) pp$PIP[1] else NA_real_,
      snps = paste(
        cs_snps,
        collapse = ";"
      )
    )
  }

  rbindlist(
    out,
    fill = TRUE
  )
}

finemap_summary_row <- function(
  fit,
  test_id,
  trait
) {

  pip_dt <- extract_pip(
    fit,
    test_id,
    trait
  )

  cs_dt <- extract_cs(
    fit,
    test_id,
    trait
  )

  data.table(
    test_id = test_id,
    trait = trait,
    converged = isTRUE(
      fit$converged
    ),
    n_variables = length(
      fit$pip
    ),
    n_credible_sets = if (
      is.null(
        fit$sets$cs
      )
    ) 0L else length(
      fit$sets$cs
    ),
    top_snp = if (
      nrow(pip_dt)
    ) pip_dt$SNP[1] else NA_character_,
    top_PIP = if (
      nrow(pip_dt)
    ) pip_dt$PIP[1] else NA_real_,
    largest_CS_size = if (
      nrow(cs_dt)
    ) max(
      cs_dt$cs_size
    ) else NA_integer_,

    # Official finite-LD / mismatch diagnostics.
    Q_art = diag_scalar(
      fit,
      "Q_art"
    ),
    r_over_B = diag_scalar(
      fit,
      "r_over_B"
    ),
    B_corrected = diag_scalar(
      fit,
      "B_corrected"
    ),
    lambda_bias = diag_scalar(
      fit,
      "lambda_bias"
    ),
    R_sensitivity_flag = diag_flag(
      fit,
      "R_sensitivity_flag"
    ),
    R_reliability_flag = diag_flag(
      fit,
      "R_reliability_flag"
    )
  )
}

# ============================================================
# 18. coloc.susie HELPERS
# ============================================================

run_coloc_susie_safe <- function(
  fit1,
  fit2,
  p12_i
) {

  safe_call(
    coloc::coloc.susie(
      dataset1 = fit1,
      dataset2 = fit2,
      back_calculate_lbf = FALSE,
      p1 = P1,
      p2 = P2,
      p12 = p12_i
    )
  )
}

summarise_coloc_susie <- function(
  res,
  test_id,
  prior_label,
  p12_i
) {

  if (
    is.null(res) ||
    is.null(res$summary) ||
    !nrow(res$summary)
  ) {
    return(
      list(
        rows = data.table(),
        test = data.table(
          test_id = test_id,
          prior = prior_label,
          p12 = p12_i,
          n_signal_pairs = 0L,
          max_PP_H4 = NA_real_,
          max_PP_H3 = NA_real_,
          best_hit1 = NA_character_,
          best_hit2 = NA_character_,
          n_H4_ge_0_8 = 0L,
          n_H4_ge_0_9 = 0L
        )
      )
    )
  }

  s <- as.data.table(
    res$summary
  )

  if (
    !all(
      c(
        "PP.H3.abf",
        "PP.H4.abf"
      ) %in%
      names(s)
    )
  ) {
    stop(
      "Unexpected coloc.susie summary format."
    )
  }

  s[
    ,
    `:=`(
      test_id = test_id,
      prior = prior_label,
      p12 = p12_i
    )
  ]

  best_i <- which.max(
    s$PP.H4.abf
  )

  list(
    rows = s,
    test = data.table(
      test_id = test_id,
      prior = prior_label,
      p12 = p12_i,
      n_signal_pairs = nrow(s),
      max_PP_H4 = max(
        s$PP.H4.abf,
        na.rm = TRUE
      ),
      max_PP_H3 = max(
        s$PP.H3.abf,
        na.rm = TRUE
      ),
      best_hit1 = if (
        "hit1" %in%
        names(s)
      ) as.character(
        s$hit1[best_i]
      ) else NA_character_,
      best_hit2 = if (
        "hit2" %in%
        names(s)
      ) as.character(
        s$hit2[best_i]
      ) else NA_character_,
      n_H4_ge_0_8 = sum(
        s$PP.H4.abf >= 0.80,
        na.rm = TRUE
      ),
      n_H4_ge_0_9 = sum(
        s$PP.H4.abf >= 0.90,
        na.rm = TRUE
      )
    )
  )
}

# ============================================================
# 19. CHECKPOINT
# ============================================================

completed_ids <- character(0)

TEST_SUMMARY_LIST <- list()
SIGNAL_ROWS_LIST <- list()
FINEMAP_LIST <- list()
CS_LIST <- list()
PIP_LIST <- list()
LDQC_LIST <- list()
QC_LIST <- list()

if (
  RESUME_IF_POSSIBLE &&
  file.exists(STATE_FILE)
) {

  st <- readRDS(
    STATE_FILE
  )

  if (
    all(
      c(
        "completed_ids",
        "test_summary_list",
        "signal_rows_list",
        "finemap_list",
        "cs_list",
        "pip_list",
        "ldqc_list",
        "qc_list"
      ) %in%
      names(st)
    )
  ) {

    completed_ids <- st$completed_ids
    TEST_SUMMARY_LIST <- st$test_summary_list
    SIGNAL_ROWS_LIST <- st$signal_rows_list
    FINEMAP_LIST <- st$finemap_list
    CS_LIST <- st$cs_list
    PIP_LIST <- st$pip_list
    LDQC_LIST <- st$ldqc_list
    QC_LIST <- st$qc_list

    cat(
      "\nCheckpoint loaded: ",
      length(unique(completed_ids)),
      "/8 completed.\n",
      sep = ""
    )
  }
}

save_state <- function() {

  saveRDS(
    list(
      completed_ids = completed_ids,
      test_summary_list = TEST_SUMMARY_LIST,
      signal_rows_list = SIGNAL_ROWS_LIST,
      finemap_list = FINEMAP_LIST,
      cs_list = CS_LIST,
      pip_list = PIP_LIST,
      ldqc_list = LDQC_LIST,
      qc_list = QC_LIST
    ),
    STATE_FILE,
    compress = FALSE
  )
}

# ============================================================
# 20. MAIN LOOP
# ============================================================

for (
  i in seq_len(
    nrow(TESTS)
  )
) {

  tt <- TESTS[i]

  test_id_i <- as.character(
    tt$test_id
  )

  if (
    test_id_i %in%
    completed_ids
  ) next

  region_i <- as.character(
    tt$region
  )

  trait1_i <- as.character(
    tt$trait1
  )

  trait2_i <- as.character(
    tt$trait2
  )

  chr_i <- as.integer(
    tt$chr
  )

  start_i <- as.integer(
    tt$start
  )

  stop_i <- as.integer(
    tt$stop
  )

  cat(
    "\n====================================================\n",
    "[", i, "/8] ",
    test_id_i,
    "\n====================================================\n",
    sep = ""
  )

  # ----------------------------------------------------------
  # A. TRAIT DATA: aligned to original 1000G BIM A1
  # ----------------------------------------------------------

  g1 <- GWAS_REGION[[
    trait1_i
  ]][
    region == region_i &
    !is.na(orientation)
  ]

  g2 <- GWAS_REGION[[
    trait2_i
  ]][
    region == region_i &
    !is.na(orientation)
  ]

  common <- intersect(
    g1$SNP,
    g2$SNP
  )

  if (
    length(common) < MIN_SNPS
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      status = "TOO_FEW_ALIGNED_COMMON_SNPS",
      message = paste0(
        "n=",
        length(common)
      )
    )

    completed_ids <- c(
      completed_ids,
      test_id_i
    )

    save_state()
    next
  }

  # ----------------------------------------------------------
  # B. REFERENCE PREFIX
  # ----------------------------------------------------------

  ref_row <- REGION_BFILE[
    region == region_i
  ]

  if (
    nrow(ref_row) != 1L
  ) {
    stop(
      "Reference-prefix lookup failed for ",
      region_i
    )
  }

  ref_prefix <- as.character(
    ref_row$prefix
  )

  ref_n <- as.integer(
    ref_row$reference_N
  )

  # ----------------------------------------------------------
  # C. PLINK REGION SUBSET
  # ----------------------------------------------------------

  test_work <- file.path(
    WORK_DIR,
    safe_filename(
      test_id_i
    )
  )

  dir.create(
    test_work,
    recursive = TRUE,
    showWarnings = FALSE
  )

  extract_file <- file.path(
    test_work,
    "extract_snps.txt"
  )

  fwrite(
    data.table(
      SNP = common
    ),
    extract_file,
    col.names = FALSE,
    quote = FALSE
  )

  subset_prefix <- file.path(
    test_work,
    "ref_subset"
  )

  subset_log <- file.path(
    test_work,
    "plink_subset.log.txt"
  )

  run_plink(
    args = c(
      "--bfile",
      shQuote(ref_prefix),
      "--extract",
      shQuote(extract_file),
      "--chr",
      as.character(chr_i),
      "--from-bp",
      as.character(start_i),
      "--to-bp",
      as.character(stop_i),
      "--maf",
      as.character(REF_MAF_MIN),
      "--keep-allele-order",
      "--make-bed",
      "--out",
      shQuote(subset_prefix)
    ),
    logfile = subset_log
  )

  subset_bim_file <- paste0(
    subset_prefix,
    ".bim"
  )

  subset_fam_file <- paste0(
    subset_prefix,
    ".fam"
  )

  if (
    !file.exists(subset_bim_file) ||
    !file.exists(subset_fam_file)
  ) {
    stop(
      "PLINK subset files were not created for ",
      test_id_i
    )
  }

  sbim <- fread(
    subset_bim_file,
    header = FALSE,
    col.names = c(
      "CHR",
      "SNP",
      "CM",
      "BP",
      "A1",
      "A2"
    ),
    showProgress = FALSE
  )

  sbim[
    ,
    `:=`(
      SNP = tolower(
        as.character(SNP)
      ),
      A1 = toupper(
        as.character(A1)
      ),
      A2 = toupper(
        as.character(A2)
      )
    )
  ]

  if (
    nrow(sbim) < MIN_SNPS
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      status = "TOO_FEW_SNPS_AFTER_REFERENCE_MAF",
      message = paste0(
        "n=",
        nrow(sbim)
      )
    )

    completed_ids <- c(
      completed_ids,
      test_id_i
    )

    save_state()
    next
  }

  if (
    anyDuplicated(
      sbim$SNP
    )
  ) {
    stop(
      "Duplicated SNP IDs remain in PLINK subset for ",
      test_id_i
    )
  }

  # ----------------------------------------------------------
  # D. REALIGN BOTH TRAITS TO FINAL SUBSET BIM A1/A2
  # ----------------------------------------------------------

  final_ref <- sbim[
    ,
    .(
      SNP,
      FINAL_A1 = A1,
      FINAL_A2 = A2,
      BP
    )
  ]

  z1 <- merge(
    final_ref,
    g1[
      ,
      .(
        SNP,
        GWAS_A1,
        GWAS_A2,
        N,
        Z
      )
    ],
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
  )

  z2 <- merge(
    final_ref,
    g2[
      ,
      .(
        SNP,
        GWAS_A1,
        GWAS_A2,
        N,
        Z
      )
    ],
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
  )

  # Preserve exact subset BIM order.
  z1 <- z1[
    match(
      final_ref$SNP,
      SNP
    )
  ]

  z2 <- z2[
    match(
      final_ref$SNP,
      SNP
    )
  ]

  # Reuse aligner by renaming final BIM alleles to REF_A1/REF_A2.
  z1[
    ,
    `:=`(
      REF_A1 = FINAL_A1,
      REF_A2 = FINAL_A2
    )
  ]

  z2[
    ,
    `:=`(
      REF_A1 = FINAL_A1,
      REF_A2 = FINAL_A2
    )
  ]

  z1 <- align_to_ref(
    z1
  )

  z2 <- align_to_ref(
    z2
  )

  keep <- (
    !is.na(
      z1$orientation
    ) &
    !is.na(
      z2$orientation
    ) &
    is.finite(
      z1$Z_ref
    ) &
    is.finite(
      z2$Z_ref
    ) &
    is.finite(
      z1$N
    ) &
    is.finite(
      z2$N
    ) &
    z1$N > 0 &
    z2$N > 0
  )

  # If a few SNPs became incompatible, rebuild the PLINK subset
  # so LD and Z use exactly the same variants.
  if (
    !all(keep)
  ) {

    compatible_snps <- final_ref$SNP[
      keep
    ]

    if (
      length(compatible_snps) < MIN_SNPS
    ) {
      stop(
        "Too few SNPs remain after final BIM allele alignment for ",
        test_id_i
      )
    }

    fwrite(
      data.table(
        SNP = compatible_snps
      ),
      extract_file,
      col.names = FALSE,
      quote = FALSE
    )

    run_plink(
      args = c(
        "--bfile",
        shQuote(ref_prefix),
        "--extract",
        shQuote(extract_file),
        "--chr",
        as.character(chr_i),
        "--from-bp",
        as.character(start_i),
        "--to-bp",
        as.character(stop_i),
        "--maf",
        as.character(REF_MAF_MIN),
        "--keep-allele-order",
        "--make-bed",
        "--out",
        shQuote(subset_prefix)
      ),
      logfile = subset_log
    )

    sbim <- fread(
      paste0(
        subset_prefix,
        ".bim"
      ),
      header = FALSE,
      col.names = c(
        "CHR",
        "SNP",
        "CM",
        "BP",
        "A1",
        "A2"
      ),
      showProgress = FALSE
    )

    sbim[
      ,
      `:=`(
        SNP = tolower(
          as.character(SNP)
        ),
        A1 = toupper(
          as.character(A1)
        ),
        A2 = toupper(
          as.character(A2)
        )
      )
    ]

    final_ref <- sbim[
      ,
      .(
        SNP,
        FINAL_A1 = A1,
        FINAL_A2 = A2,
        BP
      )
    ]

    # Rebuild final aligned Z tables.
    z1 <- merge(
      final_ref,
      g1[
        ,
        .(
          SNP,
          GWAS_A1,
          GWAS_A2,
          N,
          Z
        )
      ],
      by = "SNP",
      all.x = TRUE,
      sort = FALSE
    )[
      match(
        final_ref$SNP,
        SNP
      )
    ]

    z2 <- merge(
      final_ref,
      g2[
        ,
        .(
          SNP,
          GWAS_A1,
          GWAS_A2,
          N,
          Z
        )
      ],
      by = "SNP",
      all.x = TRUE,
      sort = FALSE
    )[
      match(
        final_ref$SNP,
        SNP
      )
    ]

    z1[
      ,
      `:=`(
        REF_A1 = FINAL_A1,
        REF_A2 = FINAL_A2
      )
    ]

    z2[
      ,
      `:=`(
        REF_A1 = FINAL_A1,
        REF_A2 = FINAL_A2
      )
    ]

    z1 <- align_to_ref(
      z1
    )

    z2 <- align_to_ref(
      z2
    )

    if (
      anyNA(z1$orientation) ||
      anyNA(z2$orientation)
    ) {
      stop(
        "Allele mismatch still remains after PLINK subset rebuild for ",
        test_id_i
      )
    }
  }

  # ----------------------------------------------------------
  # E. SIGNED LD MATRIX FROM 1000G EUR
  # ----------------------------------------------------------

  ld_out <- file.path(
    test_work,
    "signed_ld"
  )

  ld_log <- file.path(
    test_work,
    "plink_ld.log.txt"
  )

  run_plink(
    args = c(
      "--bfile",
      shQuote(subset_prefix),
      "--keep-allele-order",
      "--r",
      "square",
      "bin4",
      "yes-really",
      "--out",
      shQuote(ld_out)
    ),
    logfile = ld_log
  )

  ld_candidates <- c(
    paste0(
      ld_out,
      ".ld.bin"
    ),
    paste0(
      ld_out,
      ".ld.bin4"
    )
  )

  ld_file <- ld_candidates[
    file.exists(
      ld_candidates
    )
  ]

  if (!length(ld_file)) {
    stop(
      "PLINK signed LD binary file was not found for ",
      test_id_i
    )
  }

  n_snp <- nrow(
    sbim
  )

  R <- read_plink_bin4_ld(
    ld_file[1],
    n_snp
  )

  snps <- as.character(
    sbim$SNP
  )

  rownames(R) <- snps
  colnames(R) <- snps

  # ----------------------------------------------------------
  # F. HARD LD QC
  # ----------------------------------------------------------

  if (
    any(
      !is.finite(R)
    )
  ) {
    stop(
      "Non-finite PLINK LD entries for ",
      test_id_i
    )
  }

  max_asym <- max(
    abs(
      R -
      t(R)
    )
  )

  diag_dev <- max(
    abs(
      diag(R) -
      1
    )
  )

  max_abs_r_raw <- max(
    abs(R)
  )

  if (
    max_asym > 1e-5
  ) {
    stop(
      "PLINK LD matrix asymmetric for ",
      test_id_i,
      "; max difference=",
      max_asym
    )
  }

  if (
    diag_dev > 1e-5
  ) {
    stop(
      "PLINK LD diagonal differs from 1 for ",
      test_id_i,
      "; max difference=",
      diag_dev
    )
  }

  if (
    max_abs_r_raw >
    1 + 1e-5
  ) {
    stop(
      "PLINK genotype-derived LD still has |r| materially >1 for ",
      test_id_i,
      "; max=",
      max_abs_r_raw
    )
  }

  # Pure floating-point cleanup only.
  R <- (
    R +
    t(R)
  ) / 2

  R[
    R > 1 &
    R <= 1 + 1e-5
  ] <- 1

  R[
    R < -1 &
    R >= -1 - 1e-5
  ] <- -1

  diag(R) <- 1

  # Final order check.
  if (
    !identical(
      z1$SNP,
      snps
    ) ||
    !identical(
      z2$SNP,
      snps
    )
  ) {
    stop(
      "Final Z / PLINK LD SNP order mismatch for ",
      test_id_i
    )
  }

  # ----------------------------------------------------------
  # G. SuSiE DATASETS
  # ----------------------------------------------------------

  d1 <- build_susie_dataset(
    trait1_i,
    z1,
    R
  )

  d2 <- build_susie_dataset(
    trait2_i,
    z2,
    R
  )

  # ----------------------------------------------------------
  # H. OPTIONAL LD-MISMATCH DIAGNOSTIC
  # ----------------------------------------------------------

  sdiag1 <- estimate_s_safe(
    d1$z,
    R,
    d1$N
  )

  sdiag2 <- estimate_s_safe(
    d2$z,
    R,
    d2$N
  )

  LDQC_LIST[[
    length(LDQC_LIST) + 1L
  ]] <- data.table(
    test_id = test_id_i,
    region = region_i,
    trait1 = trait1_i,
    trait2 = trait2_i,
    reference_prefix = ref_prefix,
    reference_N = ref_n,
    reference_MAF_min = REF_MAF_MIN,
    n_snps = n_snp,
    max_abs_r_raw = max_abs_r_raw,
    max_abs_asymmetry = max_asym,
    max_diag_deviation = diag_dev,
    R_finite_used = TRUE,
    R_mismatch_mode = R_MISMATCH_MODE,
    estimate_s_trait1 = sdiag1$value,
    estimate_s_trait1_status = sdiag1$status,
    estimate_s_trait2 = sdiag2$value,
    estimate_s_trait2_status = sdiag2$status,
    median_N_trait1 = d1$N,
    median_N_trait2 = d2$N,
    flips_trait1 = sum(
      z1$orientation == -1L,
      na.rm = TRUE
    ),
    flips_trait2 = sum(
      z2$orientation == -1L,
      na.rm = TRUE
    )
  )

  # ----------------------------------------------------------
  # I. SuSiE
  # ----------------------------------------------------------

  fit1_call <- run_susie_safe(
    d1,
    paste0(
      test_id_i,
      "_",
      trait1_i
    ),
    ref_n
  )

  fit2_call <- run_susie_safe(
    d2,
    paste0(
      test_id_i,
      "_",
      trait2_i
    ),
    ref_n
  )

  if (
    !fit1_call$ok ||
    is.null(fit1_call$value) ||
    !fit2_call$ok ||
    is.null(fit2_call$value)
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      status = "SUSIE_FAIL",
      message = compact_messages(
        c(
          paste0(
            trait1_i,
            ":",
            fit1_call$error
          ),
          fit1_call$warnings,
          fit1_call$messages,
          paste0(
            trait2_i,
            ":",
            fit2_call$error
          ),
          fit2_call$warnings,
          fit2_call$messages
        )
      )
    )

    completed_ids <- c(
      completed_ids,
      test_id_i
    )

    save_state()

    rm(
      R,
      d1,
      d2
    )

    gc()
    next
  }

  fit1 <- fit1_call$value
  fit2 <- fit2_call$value

  if (
    !isTRUE(
      fit1$converged
    ) ||
    !isTRUE(
      fit2$converged
    )
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      status = "SUSIE_NOT_CONVERGED",
      message = compact_messages(
        c(
          fit1_call$warnings,
          fit1_call$messages,
          fit2_call$warnings,
          fit2_call$messages
        )
      )
    )

    completed_ids <- c(
      completed_ids,
      test_id_i
    )

    save_state()
    next
  }

  if (SAVE_SUSIE_FITS) {

    saveRDS(
      fit1,
      file.path(
        FIT_DIR,
        paste0(
          safe_filename(
            test_id_i
          ),
          "__",
          trait1_i,
          "_susie.rds"
        )
      ),
      compress = "gzip"
    )

    saveRDS(
      fit2,
      file.path(
        FIT_DIR,
        paste0(
          safe_filename(
            test_id_i
          ),
          "__",
          trait2_i,
          "_susie.rds"
        )
      ),
      compress = "gzip"
    )
  }

  FINEMAP_LIST[[
    length(FINEMAP_LIST) + 1L
  ]] <- finemap_summary_row(
    fit1,
    test_id_i,
    trait1_i
  )

  FINEMAP_LIST[[
    length(FINEMAP_LIST) + 1L
  ]] <- finemap_summary_row(
    fit2,
    test_id_i,
    trait2_i
  )

  cs1 <- extract_cs(
    fit1,
    test_id_i,
    trait1_i
  )

  cs2 <- extract_cs(
    fit2,
    test_id_i,
    trait2_i
  )

  if (nrow(cs1)) {
    CS_LIST[[
      length(CS_LIST) + 1L
    ]] <- cs1
  }

  if (nrow(cs2)) {
    CS_LIST[[
      length(CS_LIST) + 1L
    ]] <- cs2
  }

  PIP_LIST[[
    length(PIP_LIST) + 1L
  ]] <- extract_pip(
    fit1,
    test_id_i,
    trait1_i
  )

  PIP_LIST[[
    length(PIP_LIST) + 1L
  ]] <- extract_pip(
    fit2,
    test_id_i,
    trait2_i
  )

  # ----------------------------------------------------------
  # J. coloc.susie — TWO PRIORS
  # ----------------------------------------------------------

  coloc_default <- run_coloc_susie_safe(
    fit1,
    fit2,
    P12_DEFAULT
  )

  coloc_cons <- run_coloc_susie_safe(
    fit1,
    fit2,
    P12_CONSERVATIVE
  )

  if (
    !coloc_default$ok ||
    is.null(coloc_default$value) ||
    !coloc_cons$ok ||
    is.null(coloc_cons$value)
  ) {

    QC_LIST[[
      length(QC_LIST) + 1L
    ]] <- data.table(
      test_id = test_id_i,
      status = "COLOC_SUSIE_FAIL",
      message = compact_messages(
        c(
          coloc_default$error,
          coloc_default$warnings,
          coloc_default$messages,
          coloc_cons$error,
          coloc_cons$warnings,
          coloc_cons$messages
        )
      )
    )

    completed_ids <- c(
      completed_ids,
      test_id_i
    )

    save_state()
    next
  }

  sd <- summarise_coloc_susie(
    coloc_default$value,
    test_id_i,
    "default_p12_1e-5",
    P12_DEFAULT
  )

  sc <- summarise_coloc_susie(
    coloc_cons$value,
    test_id_i,
    "conservative_p12_1e-6",
    P12_CONSERVATIVE
  )

  if (nrow(sd$rows)) {
    SIGNAL_ROWS_LIST[[
      length(SIGNAL_ROWS_LIST) + 1L
    ]] <- sd$rows
  }

  if (nrow(sc$rows)) {
    SIGNAL_ROWS_LIST[[
      length(SIGNAL_ROWS_LIST) + 1L
    ]] <- sc$rows
  }

  row_summary <- merge(
    sd$test[
      ,
      .(
        test_id,
        n_signal_pairs_default = n_signal_pairs,
        max_PP_H4_default = max_PP_H4,
        max_PP_H3_default = max_PP_H3,
        best_hit1_default = best_hit1,
        best_hit2_default = best_hit2,
        n_H4_ge_0_8_default = n_H4_ge_0_8,
        n_H4_ge_0_9_default = n_H4_ge_0_9
      )
    ],
    sc$test[
      ,
      .(
        test_id,
        n_signal_pairs_conservative = n_signal_pairs,
        max_PP_H4_conservative = max_PP_H4,
        max_PP_H3_conservative = max_PP_H3,
        best_hit1_conservative = best_hit1,
        best_hit2_conservative = best_hit2,
        n_H4_ge_0_8_conservative = n_H4_ge_0_8,
        n_H4_ge_0_9_conservative = n_H4_ge_0_9
      )
    ],
    by = "test_id",
    all = TRUE
  )

  row_summary[
    ,
    `:=`(
      region = region_i,
      trait1 = trait1_i,
      trait2 = trait2_i,
      STEP8B_class = as.character(
        tt$STEP8B_class
      ),
      n_snps = n_snp,
      reference_N = ref_n,
      susie_trait1_CS = if (
        is.null(
          fit1$sets$cs
        )
      ) 0L else length(
        fit1$sets$cs
      ),
      susie_trait2_CS = if (
        is.null(
          fit2$sets$cs
        )
      ) 0L else length(
        fit2$sets$cs
      ),
      robust_H4_ge_0_8_both_priors = (
        max_PP_H4_default >= 0.80 &
        max_PP_H4_conservative >= 0.80
      ),
      robust_H4_ge_0_9_both_priors = (
        max_PP_H4_default >= 0.90 &
        max_PP_H4_conservative >= 0.90
      )
    )
  ]

  TEST_SUMMARY_LIST[[
    length(TEST_SUMMARY_LIST) + 1L
  ]] <- row_summary

  QC_LIST[[
    length(QC_LIST) + 1L
  ]] <- data.table(
    test_id = test_id_i,
    status = "PASS",
    message = compact_messages(
      c(
        fit1_call$warnings,
        fit1_call$messages,
        fit2_call$warnings,
        fit2_call$messages,
        coloc_default$warnings,
        coloc_default$messages,
        coloc_cons$warnings,
        coloc_cons$messages
      )
    )
  )

  completed_ids <- unique(
    c(
      completed_ids,
      test_id_i
    )
  )

  save_state()

  cat(
    "PASS\n",
    "  SNPs = ",
    n_snp,
    "\n",
    "  max PP.H4 default = ",
    format(
      sd$test$max_PP_H4,
      digits = 4
    ),
    "\n",
    "  max PP.H4 conservative = ",
    format(
      sc$test$max_PP_H4,
      digits = 4
    ),
    "\n",
    sep = ""
  )

  rm(
    R,
    d1,
    d2,
    fit1,
    fit2
  )

  gc()
}

# ============================================================
# 21. COMBINE OUTPUTS
# ============================================================

TEST_SUMMARY <- if (
  length(TEST_SUMMARY_LIST)
) {
  rbindlist(
    TEST_SUMMARY_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

SIGNAL_ROWS <- if (
  length(SIGNAL_ROWS_LIST)
) {
  rbindlist(
    SIGNAL_ROWS_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

FINEMAP <- if (
  length(FINEMAP_LIST)
) {
  rbindlist(
    FINEMAP_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

CS <- if (
  length(CS_LIST)
) {
  rbindlist(
    CS_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

PIP <- if (
  length(PIP_LIST)
) {
  rbindlist(
    PIP_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

LDQC <- if (
  length(LDQC_LIST)
) {
  rbindlist(
    LDQC_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

QC <- if (
  length(QC_LIST)
) {
  rbindlist(
    QC_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (nrow(TEST_SUMMARY)) {

  TEST_SUMMARY[
    ,
    evidence_class :=
      fifelse(
        robust_H4_ge_0_9_both_priors,
        "Robust multi-signal H4 >=0.9 at both priors",
        fifelse(
          robust_H4_ge_0_8_both_priors,
          "Robust multi-signal H4 >=0.8 at both priors",
          fifelse(
            max_PP_H4_default >= 0.80,
            "H4 >=0.8 only at default prior",
            fifelse(
              max_PP_H3_default >= 0.80,
              "Predominantly H3 among tested SuSiE signal pairs",
              "No dominant H3/H4 signal-pair conclusion"
            )
          )
        )
      )
  ]
}

# ============================================================
# 22. SAVE OUTPUTS
# ============================================================

TEST_SUMMARY_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_coloc_susie_test_summary.csv"
)

SIGNAL_ROWS_FILE <- file.path(
  RESULT_DIR,
  "STEP8C_V3_coloc_susie_signal_pairs.csv"
)

FINEMAP_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_susie_finemap_summary.csv"
)

CS_FILE <- file.path(
  RESULT_DIR,
  "STEP8C_V3_susie_credible_sets.csv"
)

PIP_FILE <- file.path(
  RESULT_DIR,
  "STEP8C_V3_susie_all_PIPs.csv.gz"
)

LDQC_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_LD_QC.csv"
)

QC_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_processing_QC.csv"
)

fwrite(
  TEST_SUMMARY,
  TEST_SUMMARY_FILE
)

fwrite(
  SIGNAL_ROWS,
  SIGNAL_ROWS_FILE
)

fwrite(
  FINEMAP,
  FINEMAP_FILE
)

# Compact mismatch-aware reliability table.
RELIABILITY_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_susie_reliability_QC.csv"
)

RELIABILITY <- FINEMAP[
  ,
  .(
    test_id,
    trait,
    converged,
    n_credible_sets,
    top_snp,
    top_PIP,
    Q_art,
    r_over_B,
    B_corrected,
    lambda_bias,
    R_sensitivity_flag,
    R_reliability_flag
  )
]

fwrite(
  RELIABILITY,
  RELIABILITY_FILE
)

fwrite(
  CS,
  CS_FILE
)

fwrite(
  PIP,
  PIP_FILE
)

fwrite(
  LDQC,
  LDQC_FILE
)

fwrite(
  QC,
  QC_FILE
)

# ============================================================
# 23. METHOD PROVENANCE
# ============================================================

METHOD <- data.table(
  item = c(
    "LD_reference",
    "reference_population",
    "reference_expected_N",
    "LD_software",
    "LD_command",
    "reference_MAF_min",
    "ambiguous_SNPs",
    "allele_orientation",
    "coloc_version",
    "susieR_version",
    "SUSIE_L",
    "estimate_residual_variance",
    "R_finite_supported",
    "R_mismatch_supported",
    "R_mismatch_mode",
    "first_maxit",
    "p12_default",
    "p12_conservative"
  ),
  value = c(
    "1000 Genomes Phase 3 PLINK genotype reference supplied by LAVA",
    "European",
    "approximately 503; exact N read from .fam",
    paste(
      plink_version_log,
      collapse = " "
    ),
    "--r square bin4 --keep-allele-order",
    REF_MAF_MIN,
    "excluded before signed alignment",
    "GWAS Z aligned to final PLINK BIM A1; reversed alleles flip Z",
    COLOC_VERSION,
    SUSIER_VERSION,
    SUSIE_L,
    FALSE,
    SUPPORTS_R_FINITE,
    SUPPORTS_R_MISMATCH,
    R_MISMATCH_MODE,
    SUSIE_FIRST_MAXIT,
    P12_DEFAULT,
    P12_CONSERVATIVE
  )
)

METHOD_FILE <- file.path(
  TABLE_DIR,
  "STEP8C_V3_method_provenance.csv"
)

fwrite(
  METHOD,
  METHOD_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP8C_V3_sessionInfo.txt"
  )
)

if (
  length(
    unique(
      completed_ids
    )
  ) >= 8L &&
  file.exists(
    STATE_FILE
  )
) {
  file.remove(
    STATE_FILE
  )
}

# ============================================================
# 24. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP8C_V3 COMPLETE — mismatch-aware 1000G EUR PLINK LD + SuSiE\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nProcessing status:\n"
)

print(
  QC[
    ,
    .N,
    by = status
  ]
)

if (nrow(TEST_SUMMARY)) {

  cat(
    "\nTest-level result:\n"
  )

  print(
    TEST_SUMMARY[
      ,
      .(
        test_id,
        STEP8B_class,
        n_snps,
        reference_N,
        susie_trait1_CS,
        susie_trait2_CS,
        max_PP_H4_default,
        max_PP_H4_conservative,
        max_PP_H3_default,
        best_hit1_default,
        best_hit2_default,
        evidence_class
      )
    ]
  )
}

cat(
  "\nUPLOAD THESE 8 FILES:\n",
  "1) ", TEST_SUMMARY_FILE, "\n",
  "2) ", SIGNAL_ROWS_FILE, "\n",
  "3) ", FINEMAP_FILE, "\n",
  "4) ", CS_FILE, "\n",
  "5) ", LDQC_FILE, "\n",
  "6) ", QC_FILE, "\n",
  "7) ", ALIGN_QC_FILE, "\n",
  "8) ", RELIABILITY_FILE, "\n",
  sep = ""
)

cat(
  "\nIf processing_QC contains any SUSIE_FAIL or COLOC_SUSIE_FAIL, ",
  "upload that QC file before changing any parameter.\n",
  sep = ""
)

cat(
  "====================================================\n"
)
