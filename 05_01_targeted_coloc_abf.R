# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP8A_V2_TARGETED_COLOC_ABF_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP8A_V2 — TARGETED COLOCALIZATION SCREENING WITH coloc.abf
# Project: AF / HFpEF / BMI / OSA
#
# Scientific design:
#   Regions and trait-pairs are PRE-SPECIFIED from independent
#   upstream analyses (FUMA and LAVA), before coloc is run.
#
#   Region 1: chr11 recurrent FUMA region
#     AF-HFpEF, AF-BMI, AF-OSA, BMI-OSA
#
#   Region 2: chr16 LAVA locus 2135 (FTO-overlapping block)
#     AF-BMI, AF-OSA, BMI-OSA
#
#   Region 3: chr3 BMI-centered recurrent FUMA region
#     AF-BMI, HFpEF-BMI, BMI-OSA
#
#   Region 4: chr2 AF-HFpEF FUMA region
#     AF-HFpEF
#
# Method:
#   - Uses ALL available regional SNPs (NO p-value thresholding).
#   - Uses pairwise intersections, not the four-trait intersection.
#   - Uses Z-derived two-sided p-values from STEP7B LAVA files.
#   - Uses MAF from the UKB v1.1 LAVA reference .info files.
#   - coloc.abf is a screening analysis under the single causal
#     variant assumption.
#   - Regions with convincing/ambiguous H4 vs H3 evidence should
#     subsequently undergo coloc.susie with an LD matrix.
#
# Output:
#   STEP8A_coloc_ABF_summary.csv
#   STEP8A_coloc_ABF_all_snp_posteriors.csv.gz
#   STEP8A_region_extraction_QC.csv
#   STEP8A_test_manifest.csv
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")

SUM_DIR <- file.path(
  LAVA_ROOT,
  "01_sumstats"
)

REF_DIR <- file.path(
  LAVA_ROOT,
  "00_reference",
  "UKB_v1.1"
)

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP8_COLOC"
)

REGION_DIR <- file.path(
  OUT_ROOT,
  "01_regional_sumstats"
)

RESULT_DIR <- file.path(
  OUT_ROOT,
  "02_coloc_abf"
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
  REGION_DIR,
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

TRAIT_FILES <- c(
  AF = file.path(
    SUM_DIR,
    "AF_LAVA.tsv.gz"
  ),
  HFpEF = file.path(
    SUM_DIR,
    "HFpEF_LAVA.tsv.gz"
  ),
  BMI = file.path(
    SUM_DIR,
    "BMI_LAVA.tsv.gz"
  ),
  OSA = file.path(
    SUM_DIR,
    "OSA_LAVA.tsv.gz"
  )
)

if (any(!file.exists(TRAIT_FILES))) {
  stop(
    "Missing STEP7B LAVA summary-stat file(s):\n",
    paste(
      TRAIT_FILES[!file.exists(TRAIT_FILES)],
      collapse = "\n"
    )
  )
}

INFO_FILE <- file.path(
  LAVA_ROOT,
  "02_input",
  "input.info.txt"
)

if (!file.exists(INFO_FILE)) {
  stop("Missing input.info.txt: ", INFO_FILE)
}

# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c(
  "data.table",
  "DBI",
  "duckdb",
  "coloc"
)

for (p in pkgs) {
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

cat(
  "\ncoloc version: ",
  as.character(packageVersion("coloc")),
  "\n",
  sep = ""
)

# ============================================================
# 2. PRE-SPECIFIED TARGET REGIONS
#
# hg19 / GRCh37
#
# chr11 and chr2:
#   ~1 Mb windows centered on the lead signals already identified
#   by AF-HFpEF/recurrent FUMA analysis.
#
# chr16:
#   exact LAVA locus 2135 boundaries.
#
# chr3:
#   region covering the previously frozen recurrent FUMA cluster
#   with flanking sequence; this is fixed before coloc.
# ============================================================

REGIONS <- data.table(
  region = c(
    "chr11_recurrent",
    "chr16_locus2135_FTO",
    "chr3_BMI_centered",
    "chr2_AF_HFpEF"
  ),
  chr = c(
    11L,
    16L,
    3L,
    2L
  ),
  start = c(
    16024798L,
    53393883L,
    185170000L,
    200401427L
  ),
  stop = c(
    17024798L,
    54866095L,
    186170000L,
    201401427L
  ),
  definition = c(
    "rs11529589-centered +/-500 kb window covering recurrent chr11 FUMA overlap",
    "Exact LAVA locus 2135",
    "Fixed 1 Mb window covering recurrent chr3 BMI-centered FUMA overlap",
    "rs203768-centered +/-500 kb window covering AF-HFpEF chr2 FUMA locus"
  )
)

REGION_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_target_regions.csv"
)

fwrite(
  REGIONS,
  REGION_FILE
)

# ============================================================
# 3. PRE-SPECIFIED COLOC TESTS
# ============================================================

TESTS <- rbindlist(list(

  data.table(
    region = "chr11_recurrent",
    trait1 = c(
      "AF",
      "AF",
      "AF",
      "BMI"
    ),
    trait2 = c(
      "HFpEF",
      "BMI",
      "OSA",
      "OSA"
    ),
    upstream_basis = "FUMA recurrent overlap"
  ),

  data.table(
    region = "chr16_locus2135_FTO",
    trait1 = c(
      "AF",
      "AF",
      "BMI"
    ),
    trait2 = c(
      "BMI",
      "OSA",
      "OSA"
    ),
    upstream_basis = "LAVA triad: all three pairwise local rg significant"
  ),

  data.table(
    region = "chr3_BMI_centered",
    trait1 = c(
      "AF",
      "HFpEF",
      "BMI"
    ),
    trait2 = c(
      "BMI",
      "BMI",
      "OSA"
    ),
    upstream_basis = "FUMA recurrent BMI-centered overlap"
  ),

  data.table(
    region = "chr2_AF_HFpEF",
    trait1 = "AF",
    trait2 = "HFpEF",
    upstream_basis = "AF-HFpEF FUMA locus"
  )
))

TESTS[
  ,
  test_id := paste(
    region,
    trait1,
    trait2,
    sep = "__"
  )
]

TEST_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_test_manifest.csv"
)

fwrite(
  TESTS,
  TEST_FILE
)

# ============================================================
# 4. INPUT.INFO METADATA
# ============================================================

INFO <- fread(
  INFO_FILE,
  na.strings = c(
    "NA",
    ""
  )
)

required_info <- c(
  "phenotype",
  "cases",
  "controls"
)

if (!all(required_info %in% names(INFO))) {
  stop(
    "input.info lacks: ",
    paste(
      setdiff(
        required_info,
        names(INFO)
      ),
      collapse = ", "
    )
  )
}

INFO <- INFO[
  phenotype %in% names(TRAIT_FILES)
]

if (!setequal(
  INFO$phenotype,
  names(TRAIT_FILES)
)) {
  stop(
    "input.info does not contain exactly AF/HFpEF/BMI/OSA."
  )
}

# Trait types
TRAIT_TYPE <- c(
  AF = "cc",
  HFpEF = "cc",
  BMI = "quant",
  OSA = "cc"
)

CASE_PROP <- setNames(
  rep(
    NA_real_,
    4
  ),
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
      length(ca) != 1 ||
      length(co) != 1 ||
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

    CASE_PROP[[tr]] <- ca / (ca + co)
  }
}

# ============================================================
# 5. READ REFERENCE SNP INFO FOR TARGET CHROMOSOMES
#
# LAVA binary reference .info is expected to contain:
#   SNP, CHR, POS, A1, A2, NOBS, FREQ
# FREQ is needed by LAVA for binary-trait processing and is used
# here to derive MAF=min(FREQ,1-FREQ).
# ============================================================

REF_TARGETS_LIST <- list()

for (chr_i in sort(unique(REGIONS$chr))) {

  f <- file.path(
    REF_DIR,
    paste0(
      "lava-ukb-v1.1_chr",
      chr_i,
      ".info"
    )
  )

  if (!file.exists(f)) {
    stop(
      "Missing UKB reference info: ",
      f
    )
  }

  x <- fread(
    f,
    showProgress = FALSE
  )

  required_ref <- c(
    "SNP",
    "CHR",
    "POS",
    "FREQ"
  )

  if (!all(required_ref %in% names(x))) {
    stop(
      "\nUKB reference .info lacks required columns on chr",
      chr_i,
      ".\nColumns found:\n",
      paste(
        names(x),
        collapse = ", "
      ),
      "\nWe need SNP, CHR, POS, FREQ before coloc can proceed."
    )
  }

  x[
    ,
    `:=`(
      SNP = tolower(
        as.character(SNP)
      ),
      CHR = as.integer(CHR),
      POS = as.integer(POS),
      FREQ = as.numeric(FREQ)
    )
  ]

  x[
    ,
    MAF := pmin(
      FREQ,
      1 - FREQ
    )
  ]

  # IMPORTANT:
  # data.table does NOT use tidyverse-style !! evaluation.
  # Use an explicitly named scalar to avoid column/variable collisions.
  rr <- REGIONS[
    REGIONS$chr == chr_i
  ]

  if (!nrow(rr)) {
    stop(
      "Internal region-selection error: no target region defined for chr",
      chr_i
    )
  }

  for (j in seq_len(nrow(rr))) {

    r <- rr[j]

    region_name  <- as.character(r$region)
    region_chr   <- as.integer(r$chr)
    region_start <- as.integer(r$start)
    region_stop  <- as.integer(r$stop)

    y <- x[
      CHR == region_chr &
      POS >= region_start &
      POS <= region_stop &
      is.finite(MAF) &
      MAF > 0 &
      MAF <= 0.5,
      .(
        region = region_name,
        SNP,
        CHR,
        POS,
        MAF
      )
    ]

    cat(
      "Reference extraction: ",
      region_name,
      " | chr",
      region_chr,
      ":",
      region_start,
      "-",
      region_stop,
      " | SNPs=",
      nrow(y),
      "\n",
      sep = ""
    )

    if (!nrow(y)) {
      stop(
        "No UKB reference SNPs extracted for region ",
        region_name,
        " (chr",
        region_chr,
        ":",
        region_start,
        "-",
        region_stop,
        "). Check build/coordinates/reference."
      )
    }

    REF_TARGETS_LIST[[length(REF_TARGETS_LIST) + 1L]] <- y
  }
}

REF_TARGETS <- rbindlist(
  REF_TARGETS_LIST,
  fill = TRUE
)

REF_TARGETS <- unique(
  REF_TARGETS,
  by = c(
    "region",
    "SNP"
  )
)

REF_REGION_QC <- REF_TARGETS[
  ,
  .(
    reference_snps = uniqueN(SNP),
    min_pos = min(POS),
    max_pos = max(POS),
    median_MAF = median(MAF)
  ),
  by = .(
    region,
    CHR
  )
]

expected_regions <- REGIONS$region
missing_regions <- setdiff(
  expected_regions,
  REF_REGION_QC$region
)

if (length(missing_regions)) {
  stop(
    "UKB reference extraction missed target region(s): ",
    paste(missing_regions, collapse = ", ")
  )
}

if (!nrow(REF_TARGETS)) {
  stop(
    "No target-region SNPs extracted from UKB reference."
  )
}

REF_TARGET_FILE <- file.path(
  REGION_DIR,
  "STEP8A_UKB_reference_target_SNPs.tsv.gz"
)

REF_REGION_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_reference_region_QC.csv"
)

fwrite(
  REF_TARGETS,
  REF_TARGET_FILE,
  sep = "\t"
)

fwrite(
  REF_REGION_QC,
  REF_REGION_QC_FILE
)

cat("\nReference-region QC:\n")
print(REF_REGION_QC)

# ============================================================
# 6. DUCKDB: EXTRACT EACH TRAIT ONLY ON TARGET SNPs
#
# This preserves pairwise SNP coverage and does NOT require a SNP
# to be present in all four traits.
# ============================================================

DB_FILE <- file.path(
  TMP_DIR,
  "STEP8A_temp.duckdb"
)

# Start from a clean temporary DuckDB file.
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

# Write reference targets into DuckDB
dbWriteTable(
  con,
  "ref_targets",
  as.data.frame(REF_TARGETS),
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

EXTRACTED <- list()
EXTRACT_QC <- list()

for (tr in names(TRAIT_FILES)) {

  cat(
    "\nExtracting regional SNPs for ",
    tr,
    " ...\n",
    sep = ""
  )

  q <- sprintf(
"
SELECT
  r.region,
  lower(CAST(s.SNP AS VARCHAR)) AS SNP,
  r.CHR,
  r.POS,
  r.MAF,
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
      TRAIT_FILES[[tr]]
    )
  )

  x <- as.data.table(
    dbGetQuery(
      con,
      q
    )
  )

  if (!nrow(x)) {
    stop(
      "No regional SNPs extracted for ",
      tr
    )
  }

  x <- unique(
    x,
    by = c(
      "region",
      "SNP"
    )
  )

  x[
    ,
    pvalue := pmax(
      2 * pnorm(
        abs(Z),
        lower.tail = FALSE
      ),
      1e-300
    )
  ]

  x[
    ,
    phenotype := tr
  ]

  EXTRACTED[[tr]] <- x

  out <- file.path(
    REGION_DIR,
    paste0(
      "STEP8A_",
      tr,
      "_target_regions.tsv.gz"
    )
  )

  fwrite(
    x,
    out,
    sep = "\t"
  )

  EXTRACT_QC[[tr]] <- x[
    ,
    .(
      n_snps = .N,
      median_N = median(
        N,
        na.rm = TRUE
      ),
      min_N = min(
        N,
        na.rm = TRUE
      ),
      max_N = max(
        N,
        na.rm = TRUE
      ),
      min_p = min(
        pvalue,
        na.rm = TRUE
      ),
      median_MAF = median(
        MAF,
        na.rm = TRUE
      )
    ),
    by = region
  ][
    ,
    phenotype := tr
  ]
}

EXTRACT_QC_DT <- rbindlist(
  EXTRACT_QC,
  fill = TRUE
)

EXTRACT_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_region_extraction_QC.csv"
)

fwrite(
  EXTRACT_QC_DT,
  EXTRACT_QC_FILE
)

# Regional summary statistics are now in memory; DuckDB is no longer needed.
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
# 7. COLOC DATASET BUILDER
# ============================================================

build_coloc_dataset <- function(
  tr,
  dat
) {

  n_med <- round(
    median(
      dat$N,
      na.rm = TRUE
    )
  )

  if (
    !is.finite(n_med) ||
    n_med <= 0
  ) {
    stop(
      "Invalid regional N for ",
      tr
    )
  }

  d <- list(
    pvalues = dat$pvalue,
    MAF = dat$MAF,
    N = n_med,
    type = TRAIT_TYPE[[tr]],
    snp = dat$SNP,
    position = dat$POS
  )

  if (
    TRAIT_TYPE[[tr]] == "cc"
  ) {
    d$s <- CASE_PROP[[tr]]
  }

  # Official coloc format validation
  coloc::check_dataset(
    d,
    warn.minp = 1
  )

  d
}

# ============================================================
# 8. RUN PRE-SPECIFIED coloc.abf TESTS
# ============================================================

SUMMARY_LIST <- list()
SNP_POST_LIST <- list()

for (i in seq_len(nrow(TESTS))) {

  tt <- TESTS[i]

  test_id_i <- as.character(tt$test_id)
  region_i  <- as.character(tt$region)
  trait1_i  <- as.character(tt$trait1)
  trait2_i  <- as.character(tt$trait2)
  basis_i   <- as.character(tt$upstream_basis)

  cat(
    "\n[",
    i,
    "/",
    nrow(TESTS),
    "] ",
    test_id_i,
    "\n",
    sep = ""
  )

  d1 <- EXTRACTED[[trait1_i]][
    region == region_i
  ]

  d2 <- EXTRACTED[[trait2_i]][
    region == region_i
  ]

  # Pairwise intersection
  m <- merge(
    d1[
      ,
      .(
        SNP,
        CHR,
        POS,
        MAF,
        N1 = N,
        Z1 = Z,
        P1 = pvalue
      )
    ],
    d2[
      ,
      .(
        SNP,
        N2 = N,
        Z2 = Z,
        P2 = pvalue
      )
    ],
    by = "SNP",
    all = FALSE
  )

  m <- m[
    is.finite(MAF) &
    MAF > 0 &
    MAF <= 0.5 &
    is.finite(P1) &
    P1 > 0 &
    P1 <= 1 &
    is.finite(P2) &
    P2 > 0 &
    P2 <= 1 &
    is.finite(N1) &
    N1 > 0 &
    is.finite(N2) &
    N2 > 0
  ]

  m <- unique(
    m,
    by = "SNP"
  )

  setorder(
    m,
    POS
  )

  if (nrow(m) < 100) {

    SUMMARY_LIST[[length(SUMMARY_LIST) + 1L]] <- data.table(
      test_id = test_id_i,
      region = region_i,
      trait1 = trait1_i,
      trait2 = trait2_i,
      upstream_basis = basis_i,
      nsnps = nrow(m),
      status = "TOO_FEW_COMMON_SNPS",
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      H4_over_H3H4 = NA_real_,
      top_snp_H4 = NA_character_,
      top_snp_PP_H4 = NA_real_,
      credible95_n = NA_integer_
    )

    next
  }

  dat1 <- data.table(
    SNP = m$SNP,
    POS = m$POS,
    MAF = m$MAF,
    N = m$N1,
    pvalue = m$P1
  )

  dat2 <- data.table(
    SNP = m$SNP,
    POS = m$POS,
    MAF = m$MAF,
    N = m$N2,
    pvalue = m$P2
  )

  D1 <- build_coloc_dataset(
    trait1_i,
    dat1
  )

  D2 <- build_coloc_dataset(
    trait2_i,
    dat2
  )

  # Default coloc priors:
  # p1=1e-4, p2=1e-4, p12=1e-5
  res <- coloc::coloc.abf(
    dataset1 = D1,
    dataset2 = D2
  )

  sm <- res$summary
  rr <- as.data.table(
    res$results
  )

  rr[
    ,
    `:=`(
      test_id = test_id_i,
      region = region_i,
      trait1 = trait1_i,
      trait2 = trait2_i
    )
  ]

  # 95% credible set from SNP.PP.H4
  rr2 <- copy(rr)

  rr2 <- rr2[
    is.finite(SNP.PP.H4)
  ]

  setorder(
    rr2,
    -SNP.PP.H4
  )

  rr2[
    ,
    cum_PP_H4 := cumsum(
      SNP.PP.H4
    )
  ]

  credible95 <- rr2[
    shift(
      cum_PP_H4,
      fill = 0
    ) < 0.95
  ]

  top_snp <- if (nrow(rr2)) {
    as.character(
      rr2$SNP[1]
    )
  } else {
    NA_character_
  }

  top_pp <- if (nrow(rr2)) {
    rr2$SNP.PP.H4[1]
  } else {
    NA_real_
  }

  h3 <- as.numeric(
    sm["PP.H3.abf"]
  )

  h4 <- as.numeric(
    sm["PP.H4.abf"]
  )

  h4_cond <- if (
    is.finite(h3) &&
    is.finite(h4) &&
    (h3 + h4) > 0
  ) {
    h4 / (h3 + h4)
  } else {
    NA_real_
  }

  SUMMARY_LIST[[length(SUMMARY_LIST) + 1L]] <- data.table(
    test_id = test_id_i,
    region = region_i,
    trait1 = trait1_i,
    trait2 = trait2_i,
    upstream_basis = basis_i,
    nsnps = as.integer(
      sm["nsnps"]
    ),
    status = "PASS",
    PP.H0 = as.numeric(
      sm["PP.H0.abf"]
    ),
    PP.H1 = as.numeric(
      sm["PP.H1.abf"]
    ),
    PP.H2 = as.numeric(
      sm["PP.H2.abf"]
    ),
    PP.H3 = h3,
    PP.H4 = h4,
    H4_over_H3H4 = h4_cond,
    top_snp_H4 = top_snp,
    top_snp_PP_H4 = top_pp,
    credible95_n = nrow(
      credible95
    ),
    median_N_trait1 = median(
      m$N1,
      na.rm = TRUE
    ),
    median_N_trait2 = median(
      m$N2,
      na.rm = TRUE
    ),
    min_p_trait1 = min(
      m$P1,
      na.rm = TRUE
    ),
    min_p_trait2 = min(
      m$P2,
      na.rm = TRUE
    )
  )

  SNP_POST_LIST[[length(SNP_POST_LIST) + 1L]] <- rr
}

# ============================================================
# 9. COMBINE + DESCRIPTIVE FLAGS
# ============================================================

SUMMARY <- rbindlist(
  SUMMARY_LIST,
  fill = TRUE
)

SUMMARY[
  ,
  screening_interpretation :=
    fifelse(
      status != "PASS",
      status,
      fifelse(
        PP.H4 >= 0.80 &
        H4_over_H3H4 >= 0.80,
        "Strong H4 support (screening)",
        fifelse(
          PP.H3 >= 0.80,
          "Strong H3 support: distinct variants (screening)",
          fifelse(
            (PP.H3 + PP.H4) >= 0.80,
            "Both traits associated but H3/H4 unresolved",
            "Limited evidence that both traits are associated in region"
          )
        )
      )
    )
]

SNP_POST <- if (
  length(SNP_POST_LIST)
) {
  rbindlist(
    SNP_POST_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

SUMMARY_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_coloc_ABF_summary.csv"
)

SNP_POST_FILE <- file.path(
  RESULT_DIR,
  "STEP8A_coloc_ABF_all_snp_posteriors.csv.gz"
)

fwrite(
  SUMMARY,
  SUMMARY_FILE
)

fwrite(
  SNP_POST,
  SNP_POST_FILE
)

# ============================================================
# 10. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP8A_sessionInfo.txt"
  )
)

# ============================================================
# 11. FINAL REPORT
# ============================================================

cat("\n====================================================\n")
cat("STEP8A_V2 COMPLETE — coloc.abf SCREENING\n")
cat("====================================================\n")

cat("\nRegional extraction QC:\n")
print(
  EXTRACT_QC_DT[
    order(
      region,
      phenotype
    )
  ]
)

cat("\nColocalization summary:\n")
print(
  SUMMARY[
    ,
    .(
      test_id,
      nsnps,
      PP.H3,
      PP.H4,
      H4_over_H3H4,
      top_snp_H4,
      credible95_n,
      screening_interpretation
    )
  ]
)

cat("\nUPLOAD THESE 5 FILES:\n")
cat("1) ", SUMMARY_FILE, "\n", sep = "")
cat("2) ", EXTRACT_QC_FILE, "\n", sep = "")
cat("3) ", REF_REGION_QC_FILE, "\n", sep = "")
cat("4) ", TEST_FILE, "\n", sep = "")
cat("5) ", SNP_POST_FILE, "\n", sep = "")

cat(
  "\nDo NOT yet call any locus colocalized solely from PP.H4. ",
  "We will first inspect H3 vs H4 and then run prior-sensitivity / ",
  "coloc.susie where indicated.\n",
  sep = ""
)

cat("====================================================\n")
