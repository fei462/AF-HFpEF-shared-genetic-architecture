# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP9A_V3_PREPARE_SMR_GWAS_GTEX_QC_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP9A_V3 — SMR/HEIDI ENVIRONMENT + GWAS .ma + GTEx v8 LITE QC
# Project: AF / HFpEF / BMI / OSA
#
# PURPOSE
# -------
# This script PREPARES (but does not yet perform the formal
# genome-wide SMR/HEIDI association analyses).
#
# It will:
#   1) Validate the Windows SMR executable.
#   2) Detect/extract GTEx v8 cis-eQTL LITE BESD datasets.
#   3) Validate 8 pre-specified primary GTEx tissues.
#   4) Locate the existing 1000G EUR PLINK LD reference.
#   5) Locate the four original GWAS summary-statistic files.
#   6) Auto-map their columns and convert them to official SMR/GCTA
#      GWAS .ma format:
#
#         SNP A1 A2 freq b se p n
#
#   7) Preserve the original reported GWAS P value and, by default,
#      re-scale SE so that b/SE reproduces the original two-sided P.
#      This follows the SMR documentation recommendation for GWAS
#      files whose rounded b and SE do not exactly reproduce P.
#   8) Create a frozen 4-trait x 8-tissue STEP9B job plan.
#   9) Save all QC/readiness tables.
#
# IMPORTANT SCIENTIFIC NOTES
# --------------------------
# - GTEx v8 LITE is used as a SCREENING resource. It contains only
#   cis-eQTL SNPs with P < 1e-5. For final key hits with insufficient
#   HEIDI SNP coverage, the corresponding FULL tissue dataset may be
#   needed later.
# - This script does NOT use MiXeR output.
# - This script does NOT run formal SMR/HEIDI yet.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. USER / PROJECT SETTINGS
# ============================================================

ROOT <- "D:/A/data/STEP9_SMR"
DATA_ROOT <- "D:/A/data"

SMR_DIR <- file.path(
  ROOT,
  "00_software",
  "smr-1.3.1-win-x86_64"
)

GTEX_DIR <- file.path(
  ROOT,
  "03_eqtl",
  "GTEx_v8"
)

GWAS_OUT_DIR <- file.path(
  ROOT,
  "02_gwas_ma"
)

TABLE_DIR <- file.path(
  ROOT,
  "05_tables"
)

FIG_DIR <- file.path(
  ROOT,
  "06_figures"
)

LOG_DIR <- file.path(
  ROOT,
  "07_logs"
)

TMP_DIR <- file.path(
  ROOT,
  "99_tmp_STEP9A"
)

for (d in c(
  GWAS_OUT_DIR,
  TABLE_DIR,
  FIG_DIR,
  LOG_DIR,
  TMP_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ----------------------------
# Primary tissues are frozen BEFORE SMR results are observed.
# ----------------------------
PRIMARY_TISSUES <- c(
  "Heart_Atrial_Appendage",
  "Heart_Left_Ventricle",
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum",
  "Artery_Aorta",
  "Artery_Coronary",
  "Whole_Blood",
  "Lung"
)

# SMR/Portal default-like instrument threshold for the top cis-eQTL.
PEQTL_SMR <- 5e-8

# HEIDI method 1 = newer HEIDI approach in current SMR.
HEIDI_METHOD <- 1L

THREAD_NUM <- 8L

# Run --descriptive-cis as an integrity check for all 8 tissues.
RUN_DESCRIPTIVE_CIS_QC <- TRUE

# Re-scale SE to reproduce original P from the original beta.
# This affects the numerical scale of SE but preserves the original
# beta sign and exact reported P; recommended by current SMR docs when
# rounded b/SE cause P inconsistencies.
ADJUST_SE_TO_REPORTED_P <- TRUE

# Number of rows processed per chunk when writing final .ma files.
MA_CHUNK_LINES <- 200000L

# Remove temporary "slim" GWAS TSV after successful final .ma creation.
DELETE_TEMP_SLIM <- TRUE

# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c(
  "data.table",
  "DBI",
  "duckdb"
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

# ============================================================
# 2. GENERIC HELPERS
# ============================================================

norm_path <- function(x, mustWork = FALSE) {
  normalizePath(
    x,
    winslash = "/",
    mustWork = mustWork
  )
}

safe_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_.\\-]+",
    "_",
    x
  )
}

count_lines_fast <- function(f) {

  con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) {
    gzfile(f, open = "rt")
  } else {
    file(f, open = "rt")
  }

  on.exit(
    close(con),
    add = TRUE
  )

  n <- 0L

  repeat {

    z <- readLines(
      con,
      n = 100000L,
      warn = FALSE
    )

    if (!length(z)) break

    n <- n + length(z)
  }

  n
}

quote_ident <- function(x) {
  paste0(
    '"',
    gsub(
      '"',
      '""',
      x,
      fixed = TRUE
    ),
    '"'
  )
}

sql_file <- function(x) {
  paste0(
    "'",
    gsub(
      "'",
      "''",
      norm_path(
        x,
        mustWork = TRUE
      ),
      fixed = TRUE
    ),
    "'"
  )
}

canonical_name <- function(x) {
  # IMPORTANT:
  # Convert to lowercase BEFORE removing non-alphanumeric characters.
  # The previous order incorrectly removed uppercase letters from
  # names such as rsID, N_total, A1_beta, etc.
  x <- tolower(
    as.character(x)
  )

  gsub(
    "[^a-z0-9]+",
    "",
    x
  )
}

# ============================================================
# 3. SMR EXECUTABLE QC
# ============================================================

smr_candidates <- list.files(
  SMR_DIR,
  pattern = "^smr.*\\.exe$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(smr_candidates)) {
  stop(
    "SMR executable not found under:\n",
    SMR_DIR
  )
}

# Prefer the exact filename supplied by the user.
exact_smr <- file.path(
  SMR_DIR,
  "smr-1.3.1-win.exe"
)

if (file.exists(exact_smr)) {
  SMR_EXE <- norm_path(
    exact_smr,
    TRUE
  )
} else {
  SMR_EXE <- norm_path(
    smr_candidates[1],
    TRUE
  )
}

cat(
  "\nSMR executable:\n",
  SMR_EXE,
  "\n",
  sep = ""
)

# Capture banner/help. Some versions return a non-zero status when
# invoked without a complete analysis; therefore existence + later
# descriptive-cis execution are the stronger integrity checks.
SMR_BANNER <- tryCatch(
  system2(
    SMR_EXE,
    stdout = TRUE,
    stderr = TRUE
  ),
  warning = function(w) {
    conditionMessage(w)
  },
  error = function(e) {
    conditionMessage(e)
  }
)

writeLines(
  as.character(SMR_BANNER),
  file.path(
    LOG_DIR,
    "STEP9A_SMR_banner.txt"
  )
)

# ============================================================
# 4. GTEx v8 LITE — AUTO EXTRACT IF NEEDED
# ============================================================

extract_archive_once <- function(archive, outdir) {

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  marker <- file.path(
    outdir,
    ".STEP9A_EXTRACTED"
  )

  if (file.exists(marker)) {
    return(
      invisible(TRUE)
    )
  }

  cat(
    "Extracting: ",
    basename(archive),
    "\n",
    sep = ""
  )

  ok <- TRUE

  tryCatch({

    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {

      unzip(
        archive,
        exdir = outdir
      )

    } else if (
      grepl(
        "\\.(tar\\.gz|tgz|tar)$",
        archive,
        ignore.case = TRUE
      )
    ) {

      untar(
        archive,
        exdir = outdir
      )

    } else {

      ok <- FALSE
    }

  }, error = function(e) {

    ok <<- FALSE

    warning(
      "Failed to extract ",
      archive,
      ": ",
      conditionMessage(e)
    )
  })

  if (ok) {
    writeLines(
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      marker
    )
  }

  invisible(ok)
}

find_besd_files <- function() {
  list.files(
    GTEX_DIR,
    pattern = "\\.besd$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

besd_now <- find_besd_files()

if (!length(besd_now)) {

  archives <- list.files(
    GTEX_DIR,
    pattern = "\\.(zip|tar\\.gz|tgz|tar)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(archives)) {
    stop(
      "No .besd files and no extractable GTEx archive were found under:\n",
      GTEX_DIR
    )
  }

  extracted_root <- file.path(
    GTEX_DIR,
    "extracted"
  )

  # First extraction pass.
  for (a in archives) {

    target <- file.path(
      extracted_root,
      sub(
        "\\.(zip|tar\\.gz|tgz|tar)$",
        "",
        basename(a),
        ignore.case = TRUE
      )
    )

    extract_archive_once(
      a,
      target
    )
  }

  # Some GTEx bundles contain tissue-level ZIP archives inside.
  nested_archives <- list.files(
    extracted_root,
    pattern = "\\.zip$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(nested_archives)) {

    for (a in nested_archives) {

      target <- file.path(
        dirname(a),
        sub(
          "\\.zip$",
          "",
          basename(a),
          ignore.case = TRUE
        )
      )

      extract_archive_once(
        a,
        target
      )
    }
  }
}

BESD_FILES <- find_besd_files()

if (!length(BESD_FILES)) {
  stop(
    "GTEx extraction finished but no .besd files were found."
  )
}

# ============================================================
# 5. GTEx BESD TRIPLET MANIFEST
# ============================================================

gtex_manifest_list <- list()

for (besd in BESD_FILES) {

  prefix <- sub(
    "\\.besd$",
    "",
    besd,
    ignore.case = TRUE
  )

  epi <- paste0(
    prefix,
    ".epi"
  )

  esi <- paste0(
    prefix,
    ".esi"
  )

  if (
    file.exists(epi) &&
    file.exists(esi)
  ) {

    gtex_manifest_list[[
      length(gtex_manifest_list) + 1L
    ]] <- data.table(
      dataset = basename(prefix),
      prefix = norm_path(prefix),
      besd = norm_path(besd),
      epi = norm_path(epi),
      esi = norm_path(esi),
      besd_MB = round(
        file.info(besd)$size /
          1024^2,
        2
      ),
      epi_lines = count_lines_fast(epi),
      esi_lines = count_lines_fast(esi)
    )
  }
}

GTEX_MANIFEST <- if (length(gtex_manifest_list)) {
  rbindlist(
    gtex_manifest_list,
    fill = TRUE
  )
} else {
  data.table()
}

if (!nrow(GTEX_MANIFEST)) {
  stop(
    "No complete GTEx .besd/.epi/.esi triplets were found."
  )
}

GTEX_MANIFEST[
  ,
  path_search_key :=
    tolower(
      paste(
        dataset,
        prefix
      )
    )
]

GTEX_MANIFEST_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GTEx_all_dataset_manifest.csv"
)

fwrite(
  GTEX_MANIFEST,
  GTEX_MANIFEST_FILE
)

# Match frozen primary tissue names.
SELECTED_LIST <- list()

for (tissue in PRIMARY_TISSUES) {

  key <- tolower(tissue)

  hits <- GTEX_MANIFEST[
    grepl(
      key,
      path_search_key,
      fixed = TRUE
    )
  ]

  if (!nrow(hits)) {

    SELECTED_LIST[[
      length(SELECTED_LIST) + 1L
    ]] <- data.table(
      tissue = tissue,
      found = FALSE,
      prefix = NA_character_,
      besd_MB = NA_real_,
      n_probes = NA_integer_,
      n_eqtl_snps = NA_integer_
    )

    next
  }

  # Prefer exact dataset basename if possible, otherwise shortest path.
  hits[
    ,
    exact_score :=
      as.integer(
        tolower(dataset) ==
          key
      )
  ]

  hits[
    ,
    path_len :=
      nchar(prefix)
  ]

  setorder(
    hits,
    -exact_score,
    path_len
  )

  h <- hits[1]

  SELECTED_LIST[[
    length(SELECTED_LIST) + 1L
  ]] <- data.table(
    tissue = tissue,
    found = TRUE,
    prefix = h$prefix,
    besd_MB = h$besd_MB,
    n_probes = h$epi_lines,
    n_eqtl_snps = h$esi_lines
  )
}

SELECTED_GTEX <- rbindlist(
  SELECTED_LIST,
  fill = TRUE
)

if (any(!SELECTED_GTEX$found)) {

  fwrite(
    SELECTED_GTEX,
    file.path(
      TABLE_DIR,
      "STEP9A_GTEx_selected_tissues_INCOMPLETE.csv"
    )
  )

  stop(
    "Missing pre-specified GTEx tissue(s): ",
    paste(
      SELECTED_GTEX[
        found == FALSE,
        tissue
      ],
      collapse = ", "
    ),
    "\nInspect STEP9A_GTEx_all_dataset_manifest.csv."
  )
}

# Per-tissue Bonferroni reference threshold.
SELECTED_GTEX[
  ,
  per_tissue_probe_bonferroni :=
    0.05 /
    n_probes
]

# ============================================================
# 6. OPTIONAL REAL SMR BESD INTEGRITY TEST: --descriptive-cis
# ============================================================

run_command <- function(
  exe,
  args,
  logfile
) {

  status <- suppressWarnings(
    system2(
      exe,
      args = args,
      stdout = logfile,
      stderr = logfile
    )
  )

  as.integer(
    status
  )
}

SELECTED_GTEX[
  ,
  descriptive_cis_status :=
    if (
      RUN_DESCRIPTIVE_CIS_QC
    ) {
      "NOT_RUN"
    } else {
      "DISABLED"
    }
]

if (RUN_DESCRIPTIVE_CIS_QC) {

  desc_dir <- file.path(
    ROOT,
    "03_eqtl",
    "GTEx_v8",
    "STEP9A_descriptive_cis_QC"
  )

  dir.create(
    desc_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in seq_len(nrow(SELECTED_GTEX))) {

    tissue <- SELECTED_GTEX$tissue[i]
    prefix <- SELECTED_GTEX$prefix[i]

    outprefix <- file.path(
      desc_dir,
      tissue
    )

    logfile <- file.path(
      LOG_DIR,
      paste0(
        "STEP9A_descriptive_cis_",
        tissue,
        ".log"
      )
    )

    cat(
      "\nSMR --descriptive-cis QC: ",
      tissue,
      "\n",
      sep = ""
    )

    status <- run_command(
      SMR_EXE,
      c(
        "--descriptive-cis",
        "--beqtl-summary",
        shQuote(prefix),
        "--out",
        shQuote(outprefix)
      ),
      logfile
    )

    SELECTED_GTEX[
      i,
      descriptive_cis_status :=
        if (
          identical(status, 0L)
        ) {
          "PASS"
        } else {
          paste0(
            "FAIL_STATUS_",
            status
          )
        }
    ]
  }
}

SELECTED_GTEX_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GTEx_selected_tissues.csv"
)

fwrite(
  SELECTED_GTEX,
  SELECTED_GTEX_FILE
)

# ============================================================
# 7. LOCATE ORIGINAL GWAS FILES
# ============================================================

find_named_file <- function(
  candidate_basenames,
  label
) {

  exact_top <- file.path(
    DATA_ROOT,
    candidate_basenames
  )

  top_hit <- exact_top[
    file.exists(exact_top)
  ]

  if (length(top_hit)) {
    return(
      norm_path(
        top_hit[1],
        TRUE
      )
    )
  }

  all_files <- list.files(
    DATA_ROOT,
    recursive = TRUE,
    full.names = TRUE
  )

  # Do not accidentally use newly generated STEP9 output.
  all_files <- all_files[
    !grepl(
      "/STEP9_SMR/",
      norm_path(all_files),
      fixed = TRUE
    )
  ]

  bn <- basename(
    all_files
  )

  for (cand in candidate_basenames) {

    hits <- all_files[
      tolower(bn) ==
        tolower(cand)
    ]

    if (length(hits)) {

      # Prefer the shortest path.
      hits <- hits[
        order(
          nchar(hits)
        )
      ]

      return(
        norm_path(
          hits[1],
          TRUE
        )
      )
    }
  }

  stop(
    "Could not find original GWAS file for ",
    label,
    ". Expected one of:\n",
    paste(
      candidate_basenames,
      collapse = "\n"
    )
  )
}

GWAS_FILES <- c(
  AF = find_named_file(
    c(
      "GCST90624412.tsv.gz",
      "GCST90624412.tsv",
      "GCST90624412.txt.gz"
    ),
    "AF"
  ),

  HFpEF = find_named_file(
    c(
      "FORMAT-METAL_Pheno4_EUR.tsv.gz",
      "FORMAT-METAL_Pheno4_EUR.tsv",
      "FORMAT-METAL_Pheno4_EUR.txt.gz"
    ),
    "HFpEF"
  ),

  BMI = find_named_file(
    c(
      "SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
      "SNP_gwas_mc_merge_nogc.tbl.uniq"
    ),
    "BMI"
  ),

  OSA = find_named_file(
    c(
      "finngen_R9_G6_SLEEPAPNO.gz",
      "finngen_R9_G6_SLEEPAPNO.tsv.gz",
      "finngen_R9_G6_SLEEPAPNO.txt.gz"
    ),
    "OSA"
  )
)

cat(
  "\nOriginal GWAS files:\n"
)

print(
  GWAS_FILES
)

# ============================================================
# 8. GWAS COLUMN DETECTION
# ============================================================

# Priority-based matching.
pick_col <- function(
  nms,
  exact = character(0),
  contains = character(0),
  required = TRUE,
  label = ""
) {

  canon <- canonical_name(
    nms
  )

  for (p in exact) {

    hit <- which(
      canon ==
        canonical_name(p)
    )

    if (length(hit)) {
      return(
        nms[hit[1]]
      )
    }
  }

  for (p in contains) {

    key <- canonical_name(p)

    hit <- which(
      grepl(
        key,
        canon,
        fixed = TRUE
      )
    )

    if (length(hit)) {
      return(
        nms[hit[1]]
      )
    }
  }

  if (required) {
    stop(
      "Could not detect required column: ",
      label,
      "\nAvailable columns:\n",
      paste(
        nms,
        collapse = ", "
      )
    )
  }

  NA_character_
}

# Trait-specific mapping priorities are intentional.
detect_columns <- function(
  trait,
  file
) {

  hdr <- fread(
    file,
    nrows = 0,
    check.names = FALSE,
    showProgress = FALSE
  )

  nms <- names(hdr)

  if (!length(nms)) {
    stop(
      "Could not read header from ",
      file
    )
  }

  if (trait == "AF") {

    snp <- pick_col(
      nms,
      exact = c(
        "SNP",
        "rsid",
        "rsID",
        "variant_id"
      ),
      contains = c(
        "rsid"
      ),
      label = "AF SNP"
    )

    a1 <- pick_col(
      nms,
      exact = c(
        "effect_allele",
        "A1",
        "EA"
      ),
      contains = c(
        "effectallele"
      ),
      label = "AF effect allele"
    )

    a2 <- pick_col(
      nms,
      exact = c(
        "other_allele",
        "A2",
        "NEA",
        "non_effect_allele"
      ),
      contains = c(
        "otherallele",
        "noneffectallele"
      ),
      label = "AF non-effect allele"
    )

    beta <- pick_col(
      nms,
      exact = c(
        "beta",
        "effect",
        "BETA"
      ),
      contains = c(
        "beta"
      ),
      required = FALSE,
      label = "AF beta"
    )

    or_col <- pick_col(
      nms,
      exact = c(
        "odds_ratio",
        "OR"
      ),
      contains = c(
        "oddsratio"
      ),
      required = FALSE,
      label = "AF OR"
    )

    se <- pick_col(
      nms,
      exact = c(
        "standard_error",
        "se",
        "SE",
        "sebeta"
      ),
      contains = c(
        "standarderror"
      ),
      required = FALSE,
      label = "AF SE"
    )

    p <- pick_col(
      nms,
      exact = c(
        "p_value",
        "pval",
        "p",
        "P"
      ),
      contains = c(
        "pvalue"
      ),
      label = "AF P"
    )

    freq <- pick_col(
      nms,
      exact = c(
        "effect_allele_frequency",
        "eaf",
        "EAF",
        "freq"
      ),
      contains = c(
        "effectallelefrequency"
      ),
      required = FALSE,
      label = "AF EAF"
    )

    n <- pick_col(
      nms,
      exact = c(
        "n",
        "N",
        "n_total",
        "N_total",
        "sample_size"
      ),
      contains = c(
        "samplesize"
      ),
      required = FALSE,
      label = "AF N"
    )

    ncase <- pick_col(
      nms,
      exact = c(
        "n_cases",
        "N_cases",
        "cases"
      ),
      contains = c(
        "ncases"
      ),
      required = FALSE
    )

    nctrl <- pick_col(
      nms,
      exact = c(
        "n_controls",
        "N_controls",
        "controls"
      ),
      contains = c(
        "ncontrols"
      ),
      required = FALSE
    )

    N_CONSTANT <- 1840341

  } else if (trait == "HFpEF") {

    snp <- pick_col(
      nms,
      exact = c(
        "rsID",
        "rsid",
        "MarkerName",
        "SNP"
      ),
      contains = c(
        "rsid",
        "markername"
      ),
      label = "HFpEF SNP"
    )

    a1 <- pick_col(
      nms,
      exact = c(
        "Allele1",
        "A1",
        "effect_allele"
      ),
      contains = c(
        "allele1"
      ),
      label = "HFpEF effect allele"
    )

    a2 <- pick_col(
      nms,
      exact = c(
        "Allele2",
        "A2",
        "other_allele"
      ),
      contains = c(
        "allele2"
      ),
      label = "HFpEF other allele"
    )

    beta <- pick_col(
      nms,
      exact = c(
        "A1_beta",
        "beta",
        "BETA",
        "Effect"
      ),
      contains = c(
        "a1beta",
        "beta",
        "effect"
      ),
      required = FALSE,
      label = "HFpEF beta"
    )

    or_col <- pick_col(
      nms,
      exact = c(
        "OR",
        "odds_ratio"
      ),
      required = FALSE
    )

    se <- pick_col(
      nms,
      exact = c(
        "StdErr",
        "SE",
        "se",
        "standard_error"
      ),
      contains = c(
        "stderr"
      ),
      required = FALSE,
      label = "HFpEF SE"
    )

    p <- pick_col(
      nms,
      exact = c(
        "pval",
        "PVAL",
        "P-value",
        "P.value",
        "PValue",
        "p",
        "P"
      ),
      contains = c(
        "pval",
        "pvalue"
      ),
      label = "HFpEF P"
    )

    freq <- pick_col(
      nms,
      exact = c(
        "A1_freq",
        "Freq1",
        "freq",
        "EAF",
        "effect_allele_frequency"
      ),
      contains = c(
        "a1freq",
        "freq1"
      ),
      required = FALSE,
      label = "HFpEF A1 frequency"
    )

    n <- pick_col(
      nms,
      exact = c(
        "N_total",
        "n_total",
        "N",
        "n",
        "TotalSampleSize"
      ),
      contains = c(
        "ntotal",
        "totalsamplesize"
      ),
      required = FALSE,
      label = "HFpEF N"
    )

    ncase <- pick_col(
      nms,
      exact = c(
        "N_case",
        "N_cases",
        "n_case",
        "n_cases"
      ),
      contains = c(
        "ncase"
      ),
      required = FALSE,
      label = "HFpEF case N"
    )

    nctrl <- pick_col(
      nms,
      exact = c(
        "N_control",
        "N_controls",
        "n_control",
        "n_controls"
      ),
      contains = c(
        "ncontrol"
      ),
      required = FALSE
    )

    N_CONSTANT <- NA_real_

  } else if (trait == "BMI") {

    snp <- pick_col(
      nms,
      exact = c(
        "SNP",
        "MarkerName",
        "rsid"
      ),
      label = "BMI SNP"
    )

    a1 <- pick_col(
      nms,
      exact = c(
        "A1",
        "Allele1"
      ),
      label = "BMI A1"
    )

    a2 <- pick_col(
      nms,
      exact = c(
        "A2",
        "Allele2"
      ),
      label = "BMI A2"
    )

    beta <- pick_col(
      nms,
      exact = c(
        "b",
        "beta",
        "Effect"
      ),
      label = "BMI beta"
    )

    or_col <- NA_character_

    se <- pick_col(
      nms,
      exact = c(
        "se",
        "SE",
        "StdErr"
      ),
      required = FALSE,
      label = "BMI SE"
    )

    p <- pick_col(
      nms,
      exact = c(
        "p",
        "P",
        "P.value",
        "P-value"
      ),
      label = "BMI P"
    )

    freq <- pick_col(
      nms,
      exact = c(
        "Freq1",
        "freq",
        "EAF",
        "Freq.Allele1.HapMapCEU",
        "Freq.Allele1.1000G",
        "Freq1.Hapmap",
        "Freq1.1000G"
      ),
      contains = c(
        "freqallele1",
        "freq1"
      ),
      required = FALSE,
      label = "BMI A1 frequency"
    )

    n <- pick_col(
      nms,
      exact = c(
        "N",
        "n",
        "TotalSampleSize"
      ),
      required = FALSE,
      label = "BMI N"
    )

    ncase <- NA_character_
    nctrl <- NA_character_

    N_CONSTANT <- 339224

  } else if (trait == "OSA") {

    snp <- pick_col(
      nms,
      exact = c(
        "rsids",
        "rsid",
        "SNP"
      ),
      contains = c(
        "rsids"
      ),
      label = "OSA SNP"
    )

    # FinnGen beta is coded for ALT.
    a1 <- pick_col(
      nms,
      exact = c(
        "alt",
        "effect_allele",
        "A1"
      ),
      label = "OSA effect allele"
    )

    a2 <- pick_col(
      nms,
      exact = c(
        "ref",
        "other_allele",
        "A2"
      ),
      label = "OSA other allele"
    )

    beta <- pick_col(
      nms,
      exact = c(
        "beta",
        "BETA"
      ),
      label = "OSA beta"
    )

    or_col <- NA_character_

    se <- pick_col(
      nms,
      exact = c(
        "sebeta",
        "se",
        "SE",
        "standard_error"
      ),
      required = FALSE,
      label = "OSA SE"
    )

    p <- pick_col(
      nms,
      exact = c(
        "pval",
        "p_value",
        "p",
        "P"
      ),
      label = "OSA P"
    )

    freq <- pick_col(
      nms,
      exact = c(
        "af_alt",
        "effect_allele_frequency",
        "eaf",
        "EAF"
      ),
      label = "OSA ALT frequency"
    )

    n <- pick_col(
      nms,
      exact = c(
        "n",
        "N",
        "n_total",
        "N_total"
      ),
      required = FALSE
    )

    ncase <- pick_col(
      nms,
      exact = c(
        "n_cases",
        "N_cases"
      ),
      required = FALSE
    )

    nctrl <- pick_col(
      nms,
      exact = c(
        "n_controls",
        "N_controls"
      ),
      required = FALSE
    )

    N_CONSTANT <- 375657

  } else {

    stop(
      "Unsupported trait: ",
      trait
    )
  }

  if (
    is.na(beta) &&
    is.na(or_col)
  ) {
    stop(
      trait,
      ": neither beta nor OR column could be identified."
    )
  }

  if (
    is.na(se) &&
    !ADJUST_SE_TO_REPORTED_P
  ) {
    stop(
      trait,
      ": no SE column and SE-from-P adjustment is disabled."
    )
  }

  list(
    trait = trait,
    source_file = file,
    SNP = snp,
    A1 = a1,
    A2 = a2,
    BETA = beta,
    OR = or_col,
    SE = se,
    P = p,
    FREQ = freq,
    N = n,
    N_CASE = ncase,
    N_CTRL = nctrl,
    N_CONSTANT = N_CONSTANT,
    header = paste(
      nms,
      collapse = "|"
    )
  )
}

COLUMN_MAPS <- lapply(
  names(GWAS_FILES),
  function(tr) {
    detect_columns(
      tr,
      GWAS_FILES[[tr]]
    )
  }
)

names(COLUMN_MAPS) <- names(
  GWAS_FILES
)

COLUMN_MAP_DT <- rbindlist(
  lapply(
    COLUMN_MAPS,
    function(x) {

      data.table(
        trait = x$trait,
        source_file = x$source_file,
        SNP = x$SNP,
        A1 = x$A1,
        A2 = x$A2,
        BETA = x$BETA,
        OR = x$OR,
        SE = x$SE,
        P = x$P,
        FREQ = x$FREQ,
        N = x$N,
        N_CASE = x$N_CASE,
        N_CTRL = x$N_CTRL,
        N_CONSTANT = x$N_CONSTANT
      )
    }
  ),
  fill = TRUE
)

COLUMN_MAP_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GWAS_column_mapping.csv"
)

fwrite(
  COLUMN_MAP_DT,
  COLUMN_MAP_FILE
)

# Save raw source headers as an additional diagnostic. This makes any
# future source-format discrepancy immediately visible without re-reading
# multi-GB GWAS files.
GWAS_HEADER_DIAG <- rbindlist(
  lapply(
    COLUMN_MAPS,
    function(x) {
      data.table(
        trait = x$trait,
        source_file = x$source_file,
        raw_header = x$header
      )
    }
  ),
  fill = TRUE
)

GWAS_HEADER_DIAG_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GWAS_raw_headers.csv"
)

fwrite(
  GWAS_HEADER_DIAG,
  GWAS_HEADER_DIAG_FILE
)

cat(
  "\nDetected GWAS columns:\n"
)

print(
  COLUMN_MAP_DT
)

# ============================================================
# 9. MEMORY-SAFE GWAS -> TEMP SLIM TABLE USING DUCKDB
# ============================================================

DB_FILE <- file.path(
  TMP_DIR,
  "STEP9A_gwas.duckdb"
)

if (file.exists(DB_FILE)) {
  unlink(
    DB_FILE,
    force = TRUE
  )
}

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_FILE,
  read_only = FALSE
)

dbExecute(
  con,
  "SET threads=6;"
)

dbExecute(
  con,
  sprintf(
    "SET temp_directory='%s';",
    gsub(
      "'",
      "''",
      norm_path(
        TMP_DIR
      ),
      fixed = TRUE
    )
  )
)

build_sql_exprs <- function(map) {

  q <- quote_ident

  beta_expr <- if (
    !is.na(map$BETA)
  ) {
    paste0(
      "TRY_CAST(",
      q(map$BETA),
      " AS DOUBLE)"
    )
  } else {
    paste0(
      "LN(TRY_CAST(",
      q(map$OR),
      " AS DOUBLE))"
    )
  }

  se_expr <- if (
    !is.na(map$SE)
  ) {
    paste0(
      "TRY_CAST(",
      q(map$SE),
      " AS DOUBLE)"
    )
  } else {
    "CAST(NULL AS DOUBLE)"
  }

  p_expr <- paste0(
    "TRY_CAST(",
    q(map$P),
    " AS DOUBLE)"
  )

  freq_expr <- if (
    !is.na(map$FREQ)
  ) {
    paste0(
      "TRY_CAST(",
      q(map$FREQ),
      " AS DOUBLE)"
    )
  } else {
    "CAST(NULL AS DOUBLE)"
  }

  if (!is.na(map$N)) {

    n_expr <- paste0(
      "TRY_CAST(",
      q(map$N),
      " AS DOUBLE)"
    )

  } else if (
    !is.na(map$N_CASE) &&
    !is.na(map$N_CTRL)
  ) {

    n_expr <- paste0(
      "TRY_CAST(",
      q(map$N_CASE),
      " AS DOUBLE) + ",
      "TRY_CAST(",
      q(map$N_CTRL),
      " AS DOUBLE)"
    )

  } else if (
    is.finite(
      map$N_CONSTANT
    )
  ) {

    n_expr <- format(
      map$N_CONSTANT,
      scientific = FALSE,
      trim = TRUE
    )

  } else {

    stop(
      map$trait,
      ": sample size cannot be derived."
    )
  }

  list(
    beta = beta_expr,
    se = se_expr,
    p = p_expr,
    freq = freq_expr,
    n = n_expr
  )
}

create_temp_slim <- function(
  map,
  output_file
) {

  e <- build_sql_exprs(
    map
  )

  src <- sql_file(
    map$source_file
  )

  out <- gsub(
    "'",
    "''",
    norm_path(
      output_file
    ),
    fixed = TRUE
  )

  sql <- sprintf(
"
COPY (
  WITH raw AS (
    SELECT *
    FROM read_csv_auto(
      %s,
      header = TRUE,
      sample_size = 100000,
      ignore_errors = TRUE
    )
  ),
  slim AS (
    SELECT
      LOWER(TRIM(CAST(%s AS VARCHAR))) AS SNP,
      UPPER(TRIM(CAST(%s AS VARCHAR))) AS A1,
      UPPER(TRIM(CAST(%s AS VARCHAR))) AS A2,
      %s AS freq,
      %s AS b,
      %s AS se,
      %s AS p,
      %s AS n
    FROM raw
  ),
  ranked AS (
    SELECT *,
      ROW_NUMBER() OVER (
        PARTITION BY SNP
        ORDER BY
          p ASC NULLS LAST,
          n DESC NULLS LAST
      ) AS rn
    FROM slim
    WHERE
      SNP IS NOT NULL
      AND SNP <> ''
      AND A1 IS NOT NULL
      AND A2 IS NOT NULL
      AND b IS NOT NULL
      AND p IS NOT NULL
      AND n IS NOT NULL
      AND n > 0
  )
  SELECT
    SNP, A1, A2, freq, b, se, p, n
  FROM ranked
  WHERE rn = 1
)
TO '%s'
(
  HEADER,
  DELIMITER '\t'
);
",
    src,
    quote_ident(map$SNP),
    quote_ident(map$A1),
    quote_ident(map$A2),
    e$freq,
    e$beta,
    e$se,
    e$p,
    e$n,
    out
  )

  dbExecute(
    con,
    sql
  )

  if (!file.exists(output_file)) {
    stop(
      "DuckDB failed to create temporary slim GWAS file for ",
      map$trait
    )
  }

  invisible(
    output_file
  )
}

# ============================================================
# 10. STREAMING FINAL .ma WRITER
# ============================================================

finalize_ma <- function(
  trait,
  slim_file,
  ma_file
) {

  if (file.exists(ma_file)) {
    unlink(
      ma_file
    )
  }

  con_in <- file(
    slim_file,
    open = "rt"
  )

  on.exit(
    close(con_in),
    add = TRUE
  )

  header <- readLines(
    con_in,
    n = 1L,
    warn = FALSE
  )

  if (!length(header)) {
    stop(
      trait,
      ": slim file has no header."
    )
  }

  n_in <- 0L
  n_written <- 0L
  n_invalid_allele <- 0L
  n_invalid_numeric <- 0L
  n_freq_missing <- 0L
  n_rsid <- 0L
  first_write <- TRUE

  # Store only a bounded sample of pre-adjustment P discrepancies.
  logp_diff_sample <- numeric(0)

  repeat {

    lines <- readLines(
      con_in,
      n = MA_CHUNK_LINES,
      warn = FALSE
    )

    if (!length(lines)) break

    n_in <- n_in + length(lines)

    dt <- fread(
      text = paste(
        c(
          header,
          lines
        ),
        collapse = "\n"
      ),
      sep = "\t",
      na.strings = c(
        "",
        "NA",
        "NaN",
        "nan"
      ),
      showProgress = FALSE
    )

    # Explicit numeric conversion.
    for (v in c(
      "freq",
      "b",
      "se",
      "p",
      "n"
    )) {
      dt[
        ,
        (v) := suppressWarnings(
          as.numeric(
            get(v)
          )
        )
      ]
    }

    dt[
      ,
      `:=`(
        SNP = tolower(
          trimws(
            as.character(SNP)
          )
        ),
        A1 = toupper(
          trimws(
            as.character(A1)
          )
        ),
        A2 = toupper(
          trimws(
            as.character(A2)
          )
        )
      )
    ]

    valid_allele <- (
      nchar(dt$A1) == 1L &
      nchar(dt$A2) == 1L &
      dt$A1 %in% c(
        "A",
        "C",
        "G",
        "T"
      ) &
      dt$A2 %in% c(
        "A",
        "C",
        "G",
        "T"
      ) &
      dt$A1 != dt$A2
    )

    valid_numeric <- (
      nzchar(dt$SNP) &
      is.finite(dt$b) &
      is.finite(dt$p) &
      dt$p >= 0 &
      dt$p <= 1 &
      is.finite(dt$n) &
      dt$n > 0
    )

    n_invalid_allele <-
      n_invalid_allele +
      sum(
        !valid_allele,
        na.rm = TRUE
      )

    n_invalid_numeric <-
      n_invalid_numeric +
      sum(
        !valid_numeric,
        na.rm = TRUE
      )

    keep <- (
      valid_allele &
      valid_numeric
    )

    dt <- dt[
      keep
    ]

    if (!nrow(dt)) {
      next
    }

    # Frequency must correspond to A1. Missing is allowed as NA.
    dt[
      !is.finite(freq) |
      freq < 0 |
      freq > 1,
      freq := NA_real_
    ]

    n_freq_missing <-
      n_freq_missing +
      sum(
        is.na(dt$freq)
      )

    # Keep original P exactly, avoiding literal zero.
    dt[
      ,
      p := pmin(
        pmax(
          p,
          1e-300
        ),
        1
      )
    ]

    # Pre-adjustment QC: how different is P(beta/SE) from reported P?
    qc_idx <- (
      is.finite(dt$se) &
      dt$se > 0 &
      is.finite(dt$b)
    )

    if (
      any(qc_idx) &&
      length(logp_diff_sample) < 200000L
    ) {

      p_recalc <- 2 * pnorm(
        abs(
          dt$b[qc_idx] /
            dt$se[qc_idx]
        ),
        lower.tail = FALSE
      )

      p_recalc <- pmax(
        p_recalc,
        1e-300
      )

      d <- abs(
        log10(
          p_recalc
        ) -
        log10(
          dt$p[qc_idx]
        )
      )

      room <- 200000L -
        length(
          logp_diff_sample
        )

      logp_diff_sample <- c(
        logp_diff_sample,
        head(
          d[
            is.finite(d)
          ],
          room
        )
      )
    }

    # Derive / adjust SE from reported P while preserving beta.
    if (ADJUST_SE_TO_REPORTED_P) {

      zstar <- qnorm(
        dt$p / 2,
        lower.tail = FALSE
      )

      can_adjust <- (
        is.finite(zstar) &
        zstar > 0 &
        is.finite(dt$b) &
        dt$b != 0
      )

      # IMPORTANT:
      # Do NOT combine a subsetted data.table column with an
      # external full-length vector inside :=. In:
      #
      #   dt[can_adjust, se := abs(b / zstar)]
      #
      # 'b' is evaluated only on selected rows, whereas 'zstar'
      # remains full chunk length, causing RHS/LHS length mismatch.
      #
      # Compute the full replacement vector first, then assign by
      # explicit integer indices.
      new_se <- rep(
        NA_real_,
        nrow(dt)
      )

      idx_adjust <- which(
        can_adjust
      )

      if (length(idx_adjust)) {

        new_se[
          idx_adjust
        ] <- abs(
          dt$b[
            idx_adjust
          ] /
          zstar[
            idx_adjust
          ]
        )

        dt[
          idx_adjust,
          se := new_se[
            idx_adjust
          ]
        ]
      }
    }

    # Defensive QC: after optional re-scaling, SE must have
    # exactly one value per retained row.
    if (length(dt$se) != nrow(dt)) {
      stop(
        trait,
        ": internal SE length mismatch after P-based adjustment."
      )
    }

    # If a row could not be re-scaled, a valid original SE is required.
    dt <- dt[
      is.finite(se) &
      se > 0
    ]

    if (!nrow(dt)) {
      next
    }

    dt[
      ,
      n := round(n)
    ]

    n_rsid <-
      n_rsid +
      sum(
        grepl(
          "^rs[0-9]+$",
          dt$SNP,
          ignore.case = TRUE
        )
      )

    out <- dt[
      ,
      .(
        SNP,
        A1,
        A2,
        freq,
        b,
        se,
        p,
        n
      )
    ]

    fwrite(
      out,
      ma_file,
      sep = "\t",
      append = !first_write,
      col.names = first_write,
      quote = FALSE,
      na = "NA"
    )

    first_write <- FALSE

    n_written <-
      n_written +
      nrow(out)

    rm(
      dt,
      out
    )

    gc(
      verbose = FALSE
    )
  }

  if (
    !file.exists(ma_file) ||
    n_written == 0L
  ) {
    stop(
      trait,
      ": no valid rows were written to .ma."
    )
  }

  data.table(
    trait = trait,
    temp_rows = n_in,
    final_rows = n_written,
    invalid_allele_rows = n_invalid_allele,
    invalid_numeric_rows = n_invalid_numeric,
    missing_freq_rows = n_freq_missing,
    missing_freq_fraction =
      n_freq_missing /
      max(
        n_written,
        1
      ),
    rsid_fraction =
      n_rsid /
      max(
        n_written,
        1
      ),
    raw_b_se_vs_reported_p_median_abs_log10_diff =
      if (
        length(logp_diff_sample)
      ) {
        median(
          logp_diff_sample,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    raw_b_se_vs_reported_p_p95_abs_log10_diff =
      if (
        length(logp_diff_sample)
      ) {
        as.numeric(
          quantile(
            logp_diff_sample,
            0.95,
            na.rm = TRUE
          )
        )
      } else {
        NA_real_
      },
    SE_adjusted_to_reported_P =
      ADJUST_SE_TO_REPORTED_P,
    ma_file = norm_path(
      ma_file,
      TRUE
    )
  )
}

GWAS_QC_LIST <- list()

for (tr in names(COLUMN_MAPS)) {

  cat(
    "\n====================================================\n",
    "GWAS -> SMR .ma: ",
    tr,
    "\n====================================================\n",
    sep = ""
  )

  map <- COLUMN_MAPS[[tr]]

  slim <- file.path(
    TMP_DIR,
    paste0(
      tr,
      "_STEP9A_slim.tsv"
    )
  )

  ma <- file.path(
    GWAS_OUT_DIR,
    paste0(
      tr,
      ".ma"
    )
  )

  if (file.exists(slim)) {
    unlink(
      slim
    )
  }

  create_temp_slim(
    map,
    slim
  )

  qc <- finalize_ma(
    tr,
    slim,
    ma
  )

  qc[
    ,
    source_file :=
      map$source_file
  ]

  GWAS_QC_LIST[[
    length(GWAS_QC_LIST) + 1L
  ]] <- qc

  if (
    DELETE_TEMP_SLIM &&
    file.exists(slim)
  ) {
    unlink(
      slim
    )
  }

  gc(
    verbose = FALSE
  )
}

try(
  dbDisconnect(
    con,
    shutdown = TRUE
  ),
  silent = TRUE
)

GWAS_QC <- rbindlist(
  GWAS_QC_LIST,
  fill = TRUE
)

GWAS_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GWAS_conversion_QC.csv"
)

fwrite(
  GWAS_QC,
  GWAS_QC_FILE
)

# ============================================================
# 11. FINAL .ma FORMAT HARD QC
# ============================================================

MA_FORMAT_QC_LIST <- list()

for (tr in GWAS_QC$trait) {

  f <- file.path(
    GWAS_OUT_DIR,
    paste0(
      tr,
      ".ma"
    )
  )

  head_dt <- fread(
    f,
    nrows = 10000L,
    na.strings = "NA",
    showProgress = FALSE
  )

  required <- c(
    "SNP",
    "A1",
    "A2",
    "freq",
    "b",
    "se",
    "p",
    "n"
  )

  exact_header <- identical(
    names(head_dt),
    required
  )

  valid_sample <- (
    exact_header &&
    nrow(head_dt) > 0 &&
    all(
      is.finite(
        head_dt$b
      )
    ) &&
    all(
      is.finite(
        head_dt$se
      ) &
      head_dt$se > 0
    ) &&
    all(
      is.finite(
        head_dt$p
      ) &
      head_dt$p > 0 &
      head_dt$p <= 1
    ) &&
    all(
      is.finite(
        head_dt$n
      ) &
      head_dt$n > 0
    )
  )

  # After adjustment, the recomputed P should be essentially exact.
  p_check <- 2 * pnorm(
    abs(
      head_dt$b /
        head_dt$se
    ),
    lower.tail = FALSE
  )

  p_check <- pmax(
    p_check,
    1e-300
  )

  p_diff <- abs(
    log10(
      p_check
    ) -
    log10(
      pmax(
        head_dt$p,
        1e-300
      )
    )
  )

  MA_FORMAT_QC_LIST[[
    length(MA_FORMAT_QC_LIST) + 1L
  ]] <- data.table(
    trait = tr,
    exact_header = exact_header,
    first_10000_rows_valid = valid_sample,
    median_abs_log10P_difference_after_finalization =
      median(
        p_diff,
        na.rm = TRUE
      ),
    max_abs_log10P_difference_after_finalization =
      max(
        p_diff,
        na.rm = TRUE
      )
  )
}

MA_FORMAT_QC <- rbindlist(
  MA_FORMAT_QC_LIST
)

MA_FORMAT_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_GWAS_ma_format_QC.csv"
)

fwrite(
  MA_FORMAT_QC,
  MA_FORMAT_QC_FILE
)

if (
  any(
    !MA_FORMAT_QC$exact_header |
    !MA_FORMAT_QC$first_10000_rows_valid
  )
) {
  stop(
    "At least one generated .ma file failed hard format QC. ",
    "Inspect STEP9A_GWAS_ma_format_QC.csv."
  )
}

# ============================================================
# 12. LOCATE EXISTING 1000G EUR PLINK REFERENCE
# ============================================================

count_fam_samples <- function(f) {
  count_lines_fast(f)
}

fam_files <- list.files(
  DATA_ROOT,
  pattern = "\\.fam$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude derived test-specific subsets from STEP8C and STEP9 output.
fam_norm <- norm_path(
  fam_files
)

exclude_ref <- (
  grepl(
    "/00_plink_work/",
    fam_norm,
    fixed = TRUE
  ) |
  grepl(
    "/06_coloc_susie",
    fam_norm,
    fixed = TRUE
  ) |
  grepl(
    "/07_coloc_susie",
    fam_norm,
    fixed = TRUE
  ) |
  grepl(
    "/05_SER_reliability_audit/",
    fam_norm,
    fixed = TRUE
  ) |
  grepl(
    "/STEP9_SMR/",
    fam_norm,
    fixed = TRUE
  ) |
  grepl(
    "ref_subset",
    fam_norm,
    ignore.case = TRUE
  )
)

fam_files <- fam_files[
  !exclude_ref
]

REF_LIST <- list()

for (fam in fam_files) {

  prefix <- sub(
    "\\.fam$",
    "",
    fam,
    ignore.case = TRUE
  )

  bed <- paste0(
    prefix,
    ".bed"
  )

  bim <- paste0(
    prefix,
    ".bim"
  )

  if (
    !file.exists(bed) ||
    !file.exists(bim)
  ) {
    next
  }

  nfam <- tryCatch(
    count_fam_samples(fam),
    error = function(e) NA_integer_
  )

  if (
    !is.finite(nfam) ||
    nfam < 450 ||
    nfam > 550
  ) {
    next
  }

  # Only read chromosome column, minimizing memory.
  chr_dt <- tryCatch(
    fread(
      bim,
      header = FALSE,
      select = 1L,
      col.names = "CHR",
      showProgress = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(chr_dt)) {
    next
  }

  chrs <- sort(
    unique(
      as.character(
        chr_dt$CHR
      )
    )
  )

  REF_LIST[[
    length(REF_LIST) + 1L
  ]] <- data.table(
    prefix = norm_path(prefix),
    reference_N = nfam,
    n_snps = nrow(chr_dt),
    n_chromosomes = length(chrs),
    chromosomes = paste(
      chrs,
      collapse = ";"
    ),
    name_has_1000G_EUR =
      grepl(
        "1000|1kg|g1000|eur",
        prefix,
        ignore.case = TRUE
      )
  )

  rm(
    chr_dt
  )

  gc(
    verbose = FALSE
  )
}

REF_MANIFEST <- if (length(REF_LIST)) {
  rbindlist(
    REF_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (!nrow(REF_MANIFEST)) {
  stop(
    "No original ~500-sample 1000G EUR PLINK .bed/.bim/.fam ",
    "reference was found under D:/A/data."
  )
}

REF_MANIFEST[
  ,
  score :=
    100 *
      as.integer(
        name_has_1000G_EUR
      ) +
    5 *
      n_chromosomes +
    log10(
      pmax(
        n_snps,
        1
      )
    )
]

setorder(
  REF_MANIFEST,
  -score
)

# Classify reference structure.
if (
  any(
    REF_MANIFEST$n_chromosomes >= 20
  )
) {

  REFERENCE_MODE <- "whole_genome"

  PRIMARY_REF_PREFIX <-
    REF_MANIFEST[
      n_chromosomes >= 20
    ][
      order(
        -score
      )
    ]$prefix[1]

} else {

  REFERENCE_MODE <- "per_chromosome"

  PRIMARY_REF_PREFIX <- NA_character_
}

REF_MANIFEST[
  ,
  reference_mode :=
    REFERENCE_MODE
]

REF_MANIFEST_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_1000G_EUR_reference_manifest.csv"
)

fwrite(
  REF_MANIFEST,
  REF_MANIFEST_FILE
)

cat(
  "\n1000G EUR reference mode: ",
  REFERENCE_MODE,
  "\n",
  sep = ""
)

if (
  REFERENCE_MODE ==
    "whole_genome"
) {

  cat(
    "Primary whole-genome --bfile prefix:\n",
    PRIMARY_REF_PREFIX,
    "\n",
    sep = ""
  )

} else {

  cat(
    "No single whole-genome prefix detected. ",
    "STEP9B will resolve SMR reference chromosome-wise.\n",
    sep = ""
  )
}

# ============================================================
# 13. FROZEN STEP9B JOB PLAN: 4 TRAITS x 8 TISSUES
# ============================================================

JOB_LIST <- list()

for (tr in c(
  "AF",
  "HFpEF",
  "BMI",
  "OSA"
)) {

  gwas_ma <- file.path(
    GWAS_OUT_DIR,
    paste0(
      tr,
      ".ma"
    )
  )

  for (i in seq_len(nrow(SELECTED_GTEX))) {

    tissue <- SELECTED_GTEX$tissue[i]

    outprefix <- file.path(
      ROOT,
      "04_smr_results",
      tr,
      paste0(
        tr,
        "__",
        tissue
      )
    )

    JOB_LIST[[
      length(JOB_LIST) + 1L
    ]] <- data.table(
      trait = tr,
      tissue = tissue,
      gwas_ma = norm_path(
        gwas_ma,
        TRUE
      ),
      eqtl_prefix =
        SELECTED_GTEX$prefix[i],
      n_probes =
        SELECTED_GTEX$n_probes[i],
      per_tissue_probe_bonferroni =
        SELECTED_GTEX$per_tissue_probe_bonferroni[i],
      peqtl_smr =
        PEQTL_SMR,
      heidi_method =
        HEIDI_METHOD,
      thread_num =
        THREAD_NUM,
      reference_mode =
        REFERENCE_MODE,
      whole_genome_bfile =
        PRIMARY_REF_PREFIX,
      output_prefix =
        norm_path(
          outprefix
        ),
      qtl_data_tier =
        "GTEx_v8_cis-eQTL_LITE_screening"
    )
  }
}

JOB_PLAN <- rbindlist(
  JOB_LIST,
  fill = TRUE
)

JOB_PLAN_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_STEP9B_job_plan_32runs.csv"
)

fwrite(
  JOB_PLAN,
  JOB_PLAN_FILE
)

# ============================================================
# 14. READINESS SUMMARY
# ============================================================

SMR_EXE_OK <- file.exists(
  SMR_EXE
)

GTEX_8_OK <- (
  nrow(SELECTED_GTEX) ==
    length(PRIMARY_TISSUES) &&
  all(
    SELECTED_GTEX$found
  )
)

GTEX_DESCRIPTIVE_OK <- if (
  RUN_DESCRIPTIVE_CIS_QC
) {
  all(
    SELECTED_GTEX$descriptive_cis_status ==
      "PASS"
  )
} else {
  NA
}

GWAS_4_OK <- (
  nrow(GWAS_QC) == 4L &&
  all(
    GWAS_QC$final_rows >
      100000L
  ) &&
  all(
    MA_FORMAT_QC$exact_header
  ) &&
  all(
    MA_FORMAT_QC$first_10000_rows_valid
  )
)

REF_OK <- nrow(
  REF_MANIFEST
) >= 1L

READINESS <- data.table(
  component = c(
    "SMR_executable",
    "GTEx_8_primary_tissues_triplets",
    "GTEx_descriptive_cis_QC",
    "GWAS_4_SMR_ma_files",
    "1000G_EUR_PLINK_reference",
    "STEP9B_32run_job_plan"
  ),
  pass = c(
    SMR_EXE_OK,
    GTEX_8_OK,
    GTEX_DESCRIPTIVE_OK,
    GWAS_4_OK,
    REF_OK,
    nrow(JOB_PLAN) == 32L
  ),
  detail = c(
    SMR_EXE,
    paste(
      SELECTED_GTEX$tissue,
      collapse = ";"
    ),
    paste(
      SELECTED_GTEX$descriptive_cis_status,
      collapse = ";"
    ),
    paste0(
      GWAS_QC$trait,
      "=",
      GWAS_QC$final_rows,
      collapse = ";"
    ),
    paste0(
      "mode=",
      REFERENCE_MODE,
      "; candidates=",
      nrow(
        REF_MANIFEST
      )
    ),
    paste0(
      nrow(JOB_PLAN),
      " planned trait-tissue analyses"
    )
  )
)

READINESS_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_readiness.csv"
)

fwrite(
  READINESS,
  READINESS_FILE
)

# Overall readiness:
# descriptive-cis is an integrity test, so if enabled it must pass.
required_pass <- READINESS[
  component !=
    "GTEx_descriptive_cis_QC" |
    RUN_DESCRIPTIVE_CIS_QC,
  pass
]

OVERALL_READY <- all(
  required_pass,
  na.rm = FALSE
)

# ============================================================
# 15. SESSION + FINAL REPORT
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP9A_sessionInfo.txt"
  )
)

cat(
  "\n====================================================\n",
  "STEP9A_V3 COMPLETE\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

cat(
  "\nGWAS conversion QC:\n"
)

print(
  GWAS_QC[
    ,
    .(
      trait,
      final_rows,
      missing_freq_fraction,
      rsid_fraction,
      raw_b_se_vs_reported_p_median_abs_log10_diff,
      SE_adjusted_to_reported_P
    )
  ]
)

cat(
  "\nGTEx primary tissues:\n"
)

print(
  SELECTED_GTEX[
    ,
    .(
      tissue,
      n_probes,
      n_eqtl_snps,
      besd_MB,
      descriptive_cis_status
    )
  ]
)

cat(
  "\n1000G EUR reference mode: ",
  REFERENCE_MODE,
  "\n",
  sep = ""
)

cat(
  "\nOVERALL_READY_FOR_STEP9B = ",
  OVERALL_READY,
  "\n",
  sep = ""
)

cat(
  "\nUPLOAD THESE 6 FILES:\n",
  "1) ", READINESS_FILE, "\n",
  "2) ", GWAS_QC_FILE, "\n",
  "3) ", MA_FORMAT_QC_FILE, "\n",
  "4) ", SELECTED_GTEX_FILE, "\n",
  "5) ", REF_MANIFEST_FILE, "\n",
  "6) ", COLUMN_MAP_FILE, "\n",
  sep = ""
)

cat(
  "\nIf any future GWAS-column mapping issue occurs, also upload:\n",
  GWAS_HEADER_DIAG_FILE,
  "\n",
  sep = ""
)

cat(
  "\nAlso upload the Console screenshot if ",
  "OVERALL_READY_FOR_STEP9B is FALSE.\n",
  sep = ""
)

cat(
  "\nNOTE: GTEx v8 LITE is the primary SCREENING dataset. ",
  "After STEP9B, any key association with weak/failed HEIDI because ",
  "of insufficient SNP coverage will be re-tested in the FULL ",
  "corresponding tissue dataset rather than rescued by relaxing ",
  "the HEIDI criteria.\n",
  sep = ""
)

cat(
  "====================================================\n"
)
