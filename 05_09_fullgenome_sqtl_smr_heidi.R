# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP9C2_V3_FULLGENOME_GTEX_SQTL_SMR_HEIDI(1).R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP9C2_V3 — FULL-GENOME GTEx v8 cis-sQTL LITE SMR + HEIDI
# Project: AF / HFpEF / BMI / OSA
#
# IMPORTANT FIX:
# GTEx v8 sQTL LITE is split by chromosome:
#   sQTL_<TISSUE>.lite.chr1 ... chr22
#
# The previous STEP9C2_V2 accidentally selected only chr1 for each
# tissue. THIS SCRIPT:
#   A) verifies all 22 chromosome BESD triplets for each of 8 tissues
#   B) merges chr1-22 into ONE full-genome BESD per tissue using:
#        smr --besd-flist <list> --make-besd --out <merged_prefix>
#   C) runs --descriptive-cis QC on each merged tissue
#   D) runs 4 traits x 8 tissues = 32 full-genome SMR+HEIDI jobs
#   E) performs multiple-testing, HEIDI classification, target-region
#      extraction, cross-trait overlap, and eQTL+sQTL integration
#
# Official SMR note:
#   --besd-flist + --make-besd can merge multiple BESD files.
#   We DO NOT use --geno-uni because chromosome-specific .esi files
#   are not identical.
#
# This script uses NEW output directories and will NOT reuse the
# invalid chr1-only STEP9C2_V2 .smr files.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

ROOT <- "D:/A/data/STEP9_SMR"

SMR_EXE <- file.path(
  ROOT,
  "00_software",
  "smr-1.3.1-win-x86_64",
  "smr-1.3.1-win.exe"
)

SQTL_ROOT <- file.path(
  ROOT,
  "03_sqtl",
  "GTEx_v8",
  "GTEx_V8_cis_sqtl_summary_lite",
  "sQTL_besd_lite",
  "sQTL_besd_lite"
)

GWAS_DIR <- file.path(
  ROOT,
  "02_gwas_ma"
)

# Already validated in STEP8/STEP9
BFILE <- "D:/A/data/STEP8_COLOC/00_LD_REFERENCE/g1000_eur/g1000_eur"

# New directories: do not contaminate/reuse chr1-only V2 outputs
OUT_ROOT <- file.path(
  ROOT,
  "06_sqtl_results_fullgenome"
)

MERGED_DIR <- file.path(
  OUT_ROOT,
  "00_merged_sQTL"
)

RAW_RESULT_DIR <- file.path(
  OUT_ROOT,
  "01_raw_SMR"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "02_QC"
)

TABLE_DIR <- file.path(
  ROOT,
  "05_tables",
  "STEP9C2_V3"
)

LOG_DIR <- file.path(
  ROOT,
  "07_logs",
  "STEP9C2_V3"
)

for (d in c(
  OUT_ROOT,
  MERGED_DIR,
  RAW_RESULT_DIR,
  QC_DIR,
  TABLE_DIR,
  LOG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# eQTL master for final regulatory-layer integration
EQTL_MASTER <- file.path(
  ROOT,
  "05_tables",
  "STEP9B",
  "STEP9B_all_32_SMR_HEIDI_results.csv.gz"
)

# ============================================================
# 1. FROZEN ANALYSIS SETTINGS
# ============================================================

TRAITS <- c(
  "AF",
  "HFpEF",
  "BMI",
  "OSA"
)

TISSUES <- c(
  "Heart_Atrial_Appendage",
  "Heart_Left_Ventricle",
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum",
  "Artery_Aorta",
  "Artery_Coronary",
  "Whole_Blood",
  "Lung"
)

GWAS <- c(
  AF = file.path(GWAS_DIR, "AF.ma"),
  HFpEF = file.path(GWAS_DIR, "HFpEF.ma"),
  BMI = file.path(GWAS_DIR, "BMI.ma"),
  OSA = file.path(GWAS_DIR, "OSA.ma")
)

CHR_EXPECTED <- 1:22

THREAD_NUM <- 8L

MAF_THRESHOLD <- 0.01
PEQTL_SMR <- 5e-8
PEQTL_HEIDI <- 1.57e-3

HEIDI_METHOD <- 1L
HEIDI_MIN_M <- 3L
HEIDI_MAX_M <- 20L
HEIDI_P_THRESHOLD <- 0.01

RESUME_MERGES <- TRUE
RESUME_SMR <- TRUE

# ============================================================
# 2. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)

# ============================================================
# 3. HELPERS
# ============================================================

norm_path <- function(
  x,
  mustWork = FALSE
) {
  normalizePath(
    x,
    winslash = "/",
    mustWork = mustWork
  )
}

norm_prefix <- function(prefix) {

  parent <- dirname(prefix)
  stem <- basename(prefix)

  parent_norm <- normalizePath(
    parent,
    winslash = "/",
    mustWork = TRUE
  )

  paste0(
    sub("/$", "", parent_norm),
    "/",
    stem
  )
}

check_besd_prefix <- function(prefix) {

  all(
    file.exists(
      paste0(
        prefix,
        c(
          ".besd",
          ".epi",
          ".esi"
        )
      )
    )
  )
}

check_plink_prefix <- function(prefix) {

  all(
    file.exists(
      paste0(
        prefix,
        c(
          ".bed",
          ".bim",
          ".fam"
        )
      )
    )
  )
}

valid_smr_file <- function(f) {

  if (
    !file.exists(f) ||
    file.info(f)$size <= 0
  ) {
    return(FALSE)
  }

  z <- tryCatch(
    fread(
      f,
      nrows = 3L,
      showProgress = FALSE,
      na.strings = c(
        "NA",
        "NaN",
        "nan",
        "-9"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(FALSE)
  }

  nms <- tolower(
    names(z)
  )

  all(
    c(
      "probeid",
      "p_smr"
    ) %in%
      nms
  )
}

read_tail <- function(
  f,
  n = 60L
) {

  if (!file.exists(f)) {
    return("")
  }

  x <- tryCatch(
    readLines(
      f,
      warn = FALSE
    ),
    error = function(e) character(0)
  )

  paste(
    tail(x, n),
    collapse = " || "
  )
}

as_num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

pick_col <- function(
  nms,
  candidates,
  required = TRUE
) {

  low <- tolower(
    nms
  )

  for (cand in candidates) {

    idx <- which(
      low ==
        tolower(cand)
    )

    if (length(idx)) {
      return(
        nms[idx[1]]
      )
    }
  }

  if (required) {
    stop(
      "Required column missing: ",
      paste(
        candidates,
        collapse = "/"
      ),
      "\nAvailable: ",
      paste(
        nms,
        collapse = ", "
      )
    )
  }

  NA_character_
}

safe_min <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (!length(x)) {
    return(NA_real_)
  }

  min(x)
}

# ============================================================
# 4. HARD INPUT QC
# ============================================================

if (!file.exists(SMR_EXE)) {
  stop(
    "SMR executable missing:\n",
    SMR_EXE
  )
}

if (!dir.exists(SQTL_ROOT)) {
  stop(
    "sQTL root missing:\n",
    SQTL_ROOT
  )
}

if (!check_plink_prefix(BFILE)) {
  stop(
    "Invalid PLINK reference prefix:\n",
    BFILE
  )
}

if (any(!file.exists(GWAS))) {
  stop(
    "Missing GWAS .ma file(s):\n",
    paste(
      GWAS[
        !file.exists(GWAS)
      ],
      collapse = "\n"
    )
  )
}

cat(
  "\nINPUT QC PASS\n",
  "SMR: ", SMR_EXE, "\n",
  "sQTL root: ", SQTL_ROOT, "\n",
  "1000G EUR: ", BFILE, "\n",
  sep = ""
)

# ============================================================
# 5. DISCOVER ALL CHROMOSOME-SPLIT sQTL BESDs
# ============================================================

BESD_FILES <- list.files(
  SQTL_ROOT,
  pattern = "\\.besd$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(BESD_FILES)) {
  stop(
    "No sQTL .besd files found."
  )
}

SOURCE_LIST <- list()

for (f in BESD_FILES) {

  prefix <- sub(
    "\\.besd$",
    "",
    f,
    ignore.case = TRUE
  )

  if (!check_besd_prefix(prefix)) {
    next
  }

  base <- basename(prefix)

  # Expected:
  # sQTL_Heart_Atrial_Appendage.lite.chr1
  m <- regexec(
    "^sQTL_(.+)\\.lite\\.chr([0-9]+)$",
    base,
    ignore.case = FALSE,
    perl = TRUE
  )

  mm <- regmatches(
    base,
    m
  )[[1]]

  if (length(mm) != 3L) {
    next
  }

  tissue <- mm[2]
  chr <- suppressWarnings(
    as.integer(mm[3])
  )

  SOURCE_LIST[[
    length(SOURCE_LIST) + 1L
  ]] <- data.table(
    tissue = tissue,
    chr = chr,
    prefix = norm_prefix(prefix),
    besd_file = norm_path(f, TRUE)
  )
}

SOURCE <- rbindlist(
  SOURCE_LIST,
  fill = TRUE
)

if (!nrow(SOURCE)) {
  stop(
    "Could not parse chromosome-split sQTL filenames."
  )
}

SOURCE_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_all_chr_split_sQTL_manifest.csv"
)

fwrite(
  SOURCE,
  SOURCE_FILE
)

cat(
  "\nParsed chromosome-split BESD datasets: ",
  nrow(SOURCE),
  "\n",
  sep = ""
)

# ============================================================
# 6. VERIFY ALL 22 CHROMOSOMES FOR THE 8 PRE-SPECIFIED TISSUES
# ============================================================

CHR_QC_LIST <- list()

for (tissue_name in TISSUES) {

  # Use explicit vector indexing to avoid data.table name collisions.
  z <- SOURCE[
    SOURCE$tissue == tissue_name
  ]

  found_chr <- sort(
    unique(
      z$chr
    )
  )

  missing_chr <- setdiff(
    CHR_EXPECTED,
    found_chr
  )

  duplicated_chr <- z[
    ,
    .N,
    by = chr
  ][
    N != 1L,
    chr
  ]

  CHR_QC_LIST[[
    length(CHR_QC_LIST) + 1L
  ]] <- data.table(
    tissue = tissue_name,
    n_chr_files = nrow(z),
    n_unique_chr = length(found_chr),
    missing_chr = paste(
      missing_chr,
      collapse = ";"
    ),
    duplicated_chr = paste(
      duplicated_chr,
      collapse = ";"
    ),
    PASS =
      length(missing_chr) == 0L &&
      length(duplicated_chr) == 0L &&
      setequal(
        found_chr,
        CHR_EXPECTED
      )
  )
}

CHR_QC <- rbindlist(
  CHR_QC_LIST
)

CHR_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_22chr_source_QC.csv"
)

fwrite(
  CHR_QC,
  CHR_QC_FILE
)

cat(
  "\n22-chromosome source QC:\n"
)

print(
  CHR_QC
)

if (!all(CHR_QC$PASS)) {
  stop(
    "At least one selected tissue does not have exactly chr1-22. ",
    "Inspect STEP9C2_V3_22chr_source_QC.csv."
  )
}

# ============================================================
# 7. MERGE chr1-22 -> ONE FULL-GENOME BESD PER TISSUE
# ============================================================

MERGE_QC_LIST <- list()

for (ti in seq_along(TISSUES)) {

  tissue_name <- TISSUES[ti]

  z <- SOURCE[
    SOURCE$tissue ==
      tissue_name
  ][
    order(chr)
  ]

  stopifnot(
    nrow(z) == 22L,
    identical(
      z$chr,
      1:22
    )
  )

  tissue_dir <- file.path(
    MERGED_DIR,
    tissue_name
  )

  dir.create(
    tissue_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  flist <- file.path(
    tissue_dir,
    paste0(
      tissue_name,
      "_chr1_22_besd.list"
    )
  )

  # IMPORTANT: list contains BESD PREFIXES, not .besd filenames.
  writeLines(
    z$prefix,
    flist
  )

  merged_prefix <- file.path(
    tissue_dir,
    paste0(
      "sQTL_",
      tissue_name,
      ".lite.fullgenome"
    )
  )

  merged_prefix <- norm_prefix(
    merged_prefix
  )

  merge_log <- file.path(
    LOG_DIR,
    paste0(
      "MERGE__",
      tissue_name,
      ".log"
    )
  )

  existing_valid <- check_besd_prefix(
    merged_prefix
  )

  if (
    RESUME_MERGES &&
    existing_valid
  ) {

    merge_state <- "SKIPPED_EXISTING_VALID"
    merge_status <- 0L

  } else {

    # Remove stale partial triplet if present.
    unlink(
      paste0(
        merged_prefix,
        c(
          ".besd",
          ".epi",
          ".esi"
        )
      )
    )

    args <- c(
      "--besd-flist",
      shQuote(
        norm_path(
          flist,
          TRUE
        )
      ),
      "--make-besd",
      "--out",
      shQuote(
        merged_prefix
      )
    )

    cat(
      "\n[MERGE ",
      ti,
      "/8] ",
      tissue_name,
      "\n",
      sep = ""
    )

    merge_status <- suppressWarnings(
      system2(
        SMR_EXE,
        args = args,
        stdout = merge_log,
        stderr = merge_log
      )
    )

    merge_status <- as.integer(
      merge_status
    )

    if (
      identical(
        merge_status,
        0L
      ) &&
      check_besd_prefix(
        merged_prefix
      )
    ) {
      merge_state <- "PASS"
    } else {
      merge_state <- "FAIL"
    }
  }

  merged_valid <- check_besd_prefix(
    merged_prefix
  )

  n_probes <- NA_integer_
  n_snps <- NA_integer_

  if (merged_valid) {

    n_probes <- length(
      readLines(
        paste0(
          merged_prefix,
          ".epi"
        ),
        warn = FALSE
      )
    )

    n_snps <- length(
      readLines(
        paste0(
          merged_prefix,
          ".esi"
        ),
        warn = FALSE
      )
    )
  }

  MERGE_QC_LIST[[
    length(MERGE_QC_LIST) + 1L
  ]] <- data.table(
    tissue = tissue_name,
    n_input_chromosomes = nrow(z),
    merge_status_code = merge_status,
    merge_state = merge_state,
    merged_valid = merged_valid,
    n_probes_fullgenome = n_probes,
    n_snps_fullgenome = n_snps,
    merged_prefix = merged_prefix,
    flist = norm_path(
      flist,
      TRUE
    ),
    merge_log = merge_log,
    log_tail = if (
      merge_state == "FAIL"
    ) {
      read_tail(
        merge_log
      )
    } else {
      ""
    }
  )
}

MERGE_QC <- rbindlist(
  MERGE_QC_LIST,
  fill = TRUE
)

MERGE_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_merge_QC.csv"
)

fwrite(
  MERGE_QC,
  MERGE_QC_FILE
)

cat(
  "\nFull-genome BESD merge QC:\n"
)

print(
  MERGE_QC[
    ,
    .(
      tissue,
      n_input_chromosomes,
      merge_state,
      merged_valid,
      n_probes_fullgenome,
      n_snps_fullgenome
    )
  ]
)

if (
  any(
    !MERGE_QC$merged_valid
  )
) {
  stop(
    "At least one chr1-22 BESD merge failed. ",
    "Upload STEP9C2_V3_merge_QC.csv before formal SMR."
  )
}

# ============================================================
# 8. --descriptive-cis QC ON MERGED FULL-GENOME BESD
# ============================================================

MERGE_QC[
  ,
  descriptive_cis_status :=
    NA_character_
]

for (i in seq_len(nrow(MERGE_QC))) {

  tissue_name <- MERGE_QC$tissue[i]
  prefix <- MERGE_QC$merged_prefix[i]

  outprefix <- file.path(
    QC_DIR,
    paste0(
      "descriptive_fullgenome__",
      tissue_name
    )
  )

  outprefix <- norm_prefix(
    outprefix
  )

  logf <- file.path(
    LOG_DIR,
    paste0(
      "DESCRIPTIVE__",
      tissue_name,
      ".log"
    )
  )

  args <- c(
    "--descriptive-cis",
    "--beqtl-summary",
    shQuote(
      prefix
    ),
    "--out",
    shQuote(
      outprefix
    )
  )

  cat(
    "\n[DESCRIPTIVE ",
    i,
    "/8] ",
    tissue_name,
    "\n",
    sep = ""
  )

  status <- suppressWarnings(
    system2(
      SMR_EXE,
      args = args,
      stdout = logf,
      stderr = logf
    )
  )

  status <- as.integer(
    status
  )

  MERGE_QC[
    i,
    descriptive_cis_status :=
      if (
        identical(
          status,
          0L
        )
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

fwrite(
  MERGE_QC,
  MERGE_QC_FILE
)

if (
  any(
    MERGE_QC$descriptive_cis_status !=
      "PASS"
  )
) {
  stop(
    "At least one merged full-genome sQTL dataset failed ",
    "--descriptive-cis QC."
  )
}

# ============================================================
# 9. FULL-GENOME TISSUE MANIFEST + BONFERRONI
# ============================================================

SELECTED <- MERGE_QC[
  ,
  .(
    tissue,
    prefix = merged_prefix,
    n_probes = n_probes_fullgenome,
    n_snps = n_snps_fullgenome,
    descriptive_cis_status
  )
]

SELECTED[
  ,
  tissue_Bonferroni :=
    0.05 /
    n_probes
]

SELECTED_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_selected_fullgenome_sQTL_tissues.csv"
)

fwrite(
  SELECTED,
  SELECTED_FILE
)

cat(
  "\nFull-genome tissue manifest:\n"
)

print(
  SELECTED
)

# ============================================================
# 10. BUILD 32 FULL-GENOME SMR JOBS
# ============================================================

JOB_LIST <- list()

for (tr in TRAITS) {

  for (i in seq_len(nrow(SELECTED))) {

    tissue_name <- SELECTED$tissue[i]

    out_dir <- file.path(
      RAW_RESULT_DIR,
      tr
    )

    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    outprefix <- file.path(
      out_dir,
      paste0(
        tr,
        "__",
        tissue_name,
        "__FULLGENOME"
      )
    )

    JOB_LIST[[
      length(JOB_LIST) + 1L
    ]] <- data.table(
      trait = tr,
      tissue = tissue_name,
      gwas_ma = norm_path(
        GWAS[[tr]],
        TRUE
      ),
      sqtl_prefix =
        SELECTED$prefix[i],
      n_probes =
        SELECTED$n_probes[i],
      tissue_Bonferroni =
        SELECTED$tissue_Bonferroni[i],
      output_prefix =
        norm_prefix(
          outprefix
        )
    )
  }
}

JOBS <- rbindlist(
  JOB_LIST
)

stopifnot(
  nrow(JOBS) == 32L
)

JOB_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_job_plan_32_fullgenome_runs.csv"
)

fwrite(
  JOBS,
  JOB_FILE
)

# ============================================================
# 11. RUN 32 FULL-GENOME SMR + HEIDI JOBS
# ============================================================

RUN_LIST <- list()

for (i in seq_len(nrow(JOBS))) {

  tr <- JOBS$trait[i]
  tissue_name <- JOBS$tissue[i]
  outprefix <- JOBS$output_prefix[i]

  smr_file <- paste0(
    outprefix,
    ".smr"
  )

  logf <- file.path(
    LOG_DIR,
    paste0(
      sprintf("%02d", i),
      "__",
      tr,
      "__",
      tissue_name,
      "__FULLGENOME.log"
    )
  )

  cat(
    "\n====================================================\n",
    "[SMR ",
    i,
    "/32] ",
    tr,
    " x ",
    tissue_name,
    "\n",
    sep = ""
  )

  if (
    RESUME_SMR &&
    valid_smr_file(
      smr_file
    )
  ) {

    status <- 0L
    state <- "SKIPPED_EXISTING_VALID"

  } else {

    if (file.exists(smr_file)) {
      unlink(
        smr_file
      )
    }

    args <- c(
      "--bfile",
      shQuote(
        norm_prefix(
          BFILE
        )
      ),

      "--gwas-summary",
      shQuote(
        JOBS$gwas_ma[i]
      ),

      "--beqtl-summary",
      shQuote(
        JOBS$sqtl_prefix[i]
      ),

      "--maf",
      format(
        MAF_THRESHOLD,
        scientific = FALSE
      ),

      "--peqtl-smr",
      format(
        PEQTL_SMR,
        scientific = TRUE
      ),

      "--peqtl-heidi",
      format(
        PEQTL_HEIDI,
        scientific = TRUE
      ),

      "--heidi-mtd",
      as.character(
        HEIDI_METHOD
      ),

      "--heidi-min-m",
      as.character(
        HEIDI_MIN_M
      ),

      "--heidi-max-m",
      as.character(
        HEIDI_MAX_M
      ),

      "--thread-num",
      as.character(
        THREAD_NUM
      ),

      "--out",
      shQuote(
        outprefix
      )
    )

    status <- suppressWarnings(
      system2(
        SMR_EXE,
        args = args,
        stdout = logf,
        stderr = logf
      )
    )

    status <- as.integer(
      status
    )

    state <- if (
      identical(
        status,
        0L
      ) &&
      valid_smr_file(
        smr_file
      )
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  }

  valid_out <- valid_smr_file(
    smr_file
  )

  n_rows <- if (
    valid_out
  ) {
    tryCatch(
      nrow(
        fread(
          smr_file,
          showProgress = FALSE
        )
      ),
      error = function(e) NA_integer_
    )
  } else {
    NA_integer_
  }

  RUN_LIST[[
    length(RUN_LIST) + 1L
  ]] <- data.table(
    job_index = i,
    trait = tr,
    tissue = tissue_name,
    status_code = status,
    run_state = state,
    valid_smr_output = valid_out,
    n_result_rows = n_rows,
    smr_file = smr_file,
    log_file = logf,
    log_tail = if (
      state == "FAIL"
    ) {
      read_tail(
        logf
      )
    } else {
      ""
    }
  )

  cat(
    "State=",
    state,
    " | valid=",
    valid_out,
    " | rows=",
    n_rows,
    "\n",
    sep = ""
  )
}

RUN_QC <- rbindlist(
  RUN_LIST,
  fill = TRUE
)

RUN_QC_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_run_QC.csv"
)

fwrite(
  RUN_QC,
  RUN_QC_FILE
)

FAILED <- RUN_QC[
  run_state == "FAIL" |
  valid_smr_output == FALSE
]

if (nrow(FAILED)) {

  print(
    FAILED[
      ,
      .(
        job_index,
        trait,
        tissue,
        status_code,
        log_file,
        log_tail
      )
    ]
  )

  stop(
    "At least one full-genome sQTL SMR job failed. ",
    "Upload STEP9C2_V3_run_QC.csv. Do not change thresholds."
  )
}

# ============================================================
# 12. READ + STANDARDIZE 32 FULL-GENOME .smr FILES
# ============================================================

RES_LIST <- list()

for (i in seq_len(nrow(JOBS))) {

  f <- paste0(
    JOBS$output_prefix[i],
    ".smr"
  )

  x <- fread(
    f,
    showProgress = FALSE,
    na.strings = c(
      "NA",
      "NaN",
      "nan",
      "-9"
    )
  )

  if (!nrow(x)) {
    next
  }

  nms <- names(x)

  c_probe <- pick_col(
    nms,
    c("ProbeID")
  )

  c_chr <- pick_col(
    nms,
    c(
      "Probe_Chr",
      "ProbeChr"
    ),
    FALSE
  )

  c_gene <- pick_col(
    nms,
    c(
      "Gene",
      "GeneName"
    ),
    FALSE
  )

  c_bp <- pick_col(
    nms,
    c(
      "Probe_bp",
      "ProbeBP"
    ),
    FALSE
  )

  c_snp <- pick_col(
    nms,
    c(
      "topSNP",
      "SNP"
    ),
    FALSE
  )

  c_b <- pick_col(
    nms,
    c("b_SMR"),
    FALSE
  )

  c_se <- pick_col(
    nms,
    c("se_SMR"),
    FALSE
  )

  c_psmr <- pick_col(
    nms,
    c("p_SMR")
  )

  c_pheidi <- pick_col(
    nms,
    c("p_HEIDI"),
    FALSE
  )

  c_nheidi <- pick_col(
    nms,
    c("nsnp_HEIDI"),
    FALSE
  )

  c_pgwas <- pick_col(
    nms,
    c("p_GWAS"),
    FALSE
  )

  c_pqtl <- pick_col(
    nms,
    c(
      "p_eQTL",
      "p_xQTL"
    ),
    FALSE
  )

  out <- data.table(
    trait = JOBS$trait[i],
    tissue = JOBS$tissue[i],

    ProbeID =
      as.character(
        x[[c_probe]]
      ),

    Probe_Chr = if (
      !is.na(c_chr)
    ) {
      as.character(
        x[[c_chr]]
      )
    } else {
      NA_character_
    },

    Gene = if (
      !is.na(c_gene)
    ) {
      as.character(
        x[[c_gene]]
      )
    } else {
      NA_character_
    },

    Probe_bp = if (
      !is.na(c_bp)
    ) {
      as_num(
        x[[c_bp]]
      )
    } else {
      NA_real_
    },

    top_sQTL_SNP = if (
      !is.na(c_snp)
    ) {
      as.character(
        x[[c_snp]]
      )
    } else {
      NA_character_
    },

    b_SMR = if (
      !is.na(c_b)
    ) {
      as_num(
        x[[c_b]]
      )
    } else {
      NA_real_
    },

    se_SMR = if (
      !is.na(c_se)
    ) {
      as_num(
        x[[c_se]]
      )
    } else {
      NA_real_
    },

    p_SMR =
      as_num(
        x[[c_psmr]]
      ),

    p_HEIDI = if (
      !is.na(c_pheidi)
    ) {
      as_num(
        x[[c_pheidi]]
      )
    } else {
      NA_real_
    },

    nsnp_HEIDI = if (
      !is.na(c_nheidi)
    ) {
      as_num(
        x[[c_nheidi]]
      )
    } else {
      NA_real_
    },

    p_GWAS = if (
      !is.na(c_pgwas)
    ) {
      as_num(
        x[[c_pgwas]]
      )
    } else {
      NA_real_
    },

    p_sQTL = if (
      !is.na(c_pqtl)
    ) {
      as_num(
        x[[c_pqtl]]
      )
    } else {
      NA_real_
    },

    n_probes_tissue =
      JOBS$n_probes[i],

    tissue_Bonferroni =
      JOBS$tissue_Bonferroni[i]
  )

  RES_LIST[[
    length(RES_LIST) + 1L
  ]] <- out
}

ALL <- rbindlist(
  RES_LIST,
  fill = TRUE
)

if (!nrow(ALL)) {
  stop(
    "No full-genome sQTL SMR results could be read."
  )
}

ALL[
  ,
  Probe_Chr_num :=
    suppressWarnings(
      as.integer(
        gsub(
          "^chr",
          "",
          as.character(Probe_Chr),
          ignore.case = TRUE
        )
      )
    )
]

# Hard QC: should now include multiple chromosomes.
CHR_PRESENT <- sort(
  unique(
    ALL$Probe_Chr_num[
      is.finite(
        ALL$Probe_Chr_num
      )
    ]
  )
)

if (length(CHR_PRESENT) < 20L) {
  stop(
    "Full-genome result QC failed: only chromosomes ",
    paste(
      CHR_PRESENT,
      collapse = ","
    ),
    " are present."
  )
}

# ============================================================
# 13. MULTIPLE TESTING + HEIDI CLASSIFICATION
# ============================================================

ALL[
  ,
  FDR_tissue :=
    p.adjust(
      p_SMR,
      method = "BH"
    ),
  by = .(
    trait,
    tissue
  )
]

ALL[
  ,
  FDR_trait :=
    p.adjust(
      p_SMR,
      method = "BH"
    ),
  by = trait
]

ALL[
  ,
  FDR_global :=
    p.adjust(
      p_SMR,
      method = "BH"
    )
]

ALL[
  ,
  SMR_Bonferroni_pass :=
    is.finite(p_SMR) &
    p_SMR <
      tissue_Bonferroni
]

ALL[
  ,
  HEIDI_evaluable :=
    is.finite(p_HEIDI) &
    is.finite(nsnp_HEIDI) &
    nsnp_HEIDI >=
      HEIDI_MIN_M
]

ALL[
  ,
  HEIDI_status :=
    fifelse(
      !HEIDI_evaluable,
      "NOT_EVALUABLE",
      fifelse(
        p_HEIDI >
          HEIDI_P_THRESHOLD,
        "PASS_NO_HETEROGENEITY",
        "FAIL_HETEROGENEITY"
      )
    )
]

ALL[
  ,
  Primary_pass :=
    SMR_Bonferroni_pass &
    HEIDI_status ==
      "PASS_NO_HETEROGENEITY"
]

ALL[
  ,
  evidence_class :=
    fifelse(
      Primary_pass,
      "sQTL_SMR_Bonferroni_and_HEIDI_pass",

      fifelse(
        SMR_Bonferroni_pass &
        HEIDI_status ==
          "FAIL_HETEROGENEITY",
        "sQTL_SMR_Bonferroni_but_HEIDI_heterogeneity",

        fifelse(
          SMR_Bonferroni_pass &
          HEIDI_status ==
            "NOT_EVALUABLE",
          "sQTL_SMR_Bonferroni_HEIDI_not_evaluable",

          "No_primary_sQTL_SMR_evidence"
        )
      )
    )
]

MASTER_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_all_32_FULLGENOME_sQTL_SMR_HEIDI_results.csv.gz"
)

fwrite(
  ALL,
  MASTER_FILE
)

# ============================================================
# 14. MAIN RESULT TABLES
# ============================================================

PRIMARY <- ALL[
  Primary_pass == TRUE
][
  order(
    trait,
    p_SMR
  )
]

PRIMARY_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_primary_Bonferroni_HEIDI_pass.csv"
)

fwrite(
  PRIMARY,
  PRIMARY_FILE
)

HEIDI_FAIL <- ALL[
  SMR_Bonferroni_pass == TRUE &
  HEIDI_status ==
    "FAIL_HETEROGENEITY"
][
  order(
    trait,
    p_SMR
  )
]

HEIDI_FAIL_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_Bonferroni_but_HEIDI_heterogeneity.csv"
)

fwrite(
  HEIDI_FAIL,
  HEIDI_FAIL_FILE
)

HEIDI_NE <- ALL[
  SMR_Bonferroni_pass == TRUE &
  HEIDI_status ==
    "NOT_EVALUABLE"
][
  order(
    trait,
    tissue,
    p_SMR
  )
]

HEIDI_NE_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_FULL_sQTL_retest_candidates.csv"
)

fwrite(
  HEIDI_NE,
  HEIDI_NE_FILE
)

# ============================================================
# 15. FROZEN FOUR TARGET REGIONS
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
  )
)

TARGET_LIST <- list()

for (i in seq_len(nrow(REGIONS))) {

  rr <- REGIONS[i]

  z <- ALL[
    Probe_Chr_num ==
      rr$chr &
    is.finite(Probe_bp) &
    Probe_bp >=
      rr$start &
    Probe_bp <=
      rr$stop
  ]

  if (nrow(z)) {

    z[
      ,
      target_region :=
        rr$region
    ]

    TARGET_LIST[[
      length(TARGET_LIST) + 1L
    ]] <- z
  }
}

TARGET <- if (
  length(TARGET_LIST)
) {
  rbindlist(
    TARGET_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

TARGET_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_target_region_sQTL_results.csv"
)

fwrite(
  TARGET,
  TARGET_FILE
)

GRID <- CJ(
  target_region =
    REGIONS$region,
  trait =
    TRAITS,
  unique = TRUE
)

if (nrow(TARGET)) {

  T_SUM <- TARGET[
    ,
    .(
      n_tested_rows = .N,
      n_probes =
        uniqueN(
          ProbeID
        ),

      min_p_SMR =
        safe_min(
          p_SMR
        ),

      n_Bonferroni =
        sum(
          SMR_Bonferroni_pass,
          na.rm = TRUE
        ),

      n_primary_pass =
        sum(
          Primary_pass,
          na.rm = TRUE
        ),

      primary_genes =
        paste(
          sort(
            unique(
              Gene[
                Primary_pass == TRUE
              ]
            )
          ),
          collapse = ";"
        ),

      primary_splicing_events =
        paste(
          sort(
            unique(
              ProbeID[
                Primary_pass == TRUE
              ]
            )
          ),
          collapse = ";"
        )
    ),
    by = .(
      target_region,
      trait
    )
  ]

} else {
  T_SUM <- data.table()
}

TARGET_SUMMARY <- merge(
  GRID,
  T_SUM,
  by = c(
    "target_region",
    "trait"
  ),
  all.x = TRUE
)

TARGET_SUMMARY_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_target_region_summary.csv"
)

fwrite(
  TARGET_SUMMARY,
  TARGET_SUMMARY_FILE
)

# ============================================================
# 16. STRICT CROSS-TRAIT sQTL OVERLAP
#
# Primary overlap criterion:
# SAME tissue + SAME ProbeID is strictest because GTEx LeafCutter
# cluster IDs may differ across tissues.
# We additionally produce a gene-level overlap table separately.
# ============================================================

SHARED_EVENT_LIST <- list()
counter <- 0L

for (i in seq_len(length(TRAITS) - 1L)) {

  for (j in (i + 1L):length(TRAITS)) {

    a <- TRAITS[i]
    b <- TRAITS[j]

    A <- unique(
      PRIMARY[
        trait == a,
        .(
          tissue,
          ProbeID
        )
      ]
    )

    B <- unique(
      PRIMARY[
        trait == b,
        .(
          tissue,
          ProbeID
        )
      ]
    )

    S <- merge(
      A,
      B,
      by = c(
        "tissue",
        "ProbeID"
      )
    )

    if (!nrow(S)) {
      next
    }

    for (k in seq_len(nrow(S))) {

      aa <- PRIMARY[
        trait == a &
        tissue == S$tissue[k] &
        ProbeID == S$ProbeID[k]
      ][
        order(p_SMR)
      ]

      bb <- PRIMARY[
        trait == b &
        tissue == S$tissue[k] &
        ProbeID == S$ProbeID[k]
      ][
        order(p_SMR)
      ]

      counter <- counter + 1L

      SHARED_EVENT_LIST[[
        counter
      ]] <- data.table(
        trait1 = a,
        trait2 = b,
        tissue = S$tissue[k],
        ProbeID = S$ProbeID[k],
        Gene = aa$Gene[1],
        trait1_p_SMR = aa$p_SMR[1],
        trait1_p_HEIDI = aa$p_HEIDI[1],
        trait2_p_SMR = bb$p_SMR[1],
        trait2_p_HEIDI = bb$p_HEIDI[1],
        Probe_Chr = aa$Probe_Chr_num[1],
        Probe_bp = aa$Probe_bp[1]
      )
    }
  }
}

SHARED_EVENT <- if (
  length(SHARED_EVENT_LIST)
) {
  rbindlist(
    SHARED_EVENT_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

SHARED_EVENT_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_strict_cross_trait_shared_splicing_events.csv"
)

fwrite(
  SHARED_EVENT,
  SHARED_EVENT_FILE
)

# Gene-level shared strict primary signals
GENE_PRIMARY <- unique(
  PRIMARY[
    !is.na(Gene) &
    nzchar(trimws(Gene)),
    .(
      trait,
      gene = toupper(
        trimws(Gene)
      )
    )
  ]
)

SHARED_GENE_LIST <- list()
counter <- 0L

for (i in seq_len(length(TRAITS) - 1L)) {

  for (j in (i + 1L):length(TRAITS)) {

    a <- TRAITS[i]
    b <- TRAITS[j]

    shared_genes <- intersect(
      GENE_PRIMARY[
        trait == a,
        gene
      ],
      GENE_PRIMARY[
        trait == b,
        gene
      ]
    )

    if (!length(shared_genes)) {
      next
    }

    counter <- counter + 1L

    SHARED_GENE_LIST[[
      counter
    ]] <- data.table(
      trait1 = a,
      trait2 = b,
      gene = shared_genes
    )
  }
}

SHARED_GENE <- if (
  length(SHARED_GENE_LIST)
) {
  rbindlist(
    SHARED_GENE_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

SHARED_GENE_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_cross_trait_shared_sQTL_genes.csv"
)

fwrite(
  SHARED_GENE,
  SHARED_GENE_FILE
)

# ============================================================
# 17. CROSS-TISSUE RECURRENCE
# ============================================================

if (nrow(PRIMARY)) {

  RECURRENCE <- PRIMARY[
    ,
    .(
      n_tissues =
        uniqueN(
          tissue
        ),

      tissues =
        paste(
          sort(
            unique(
              tissue
            )
          ),
          collapse = ";"
        ),

      min_p_SMR =
        min(
          p_SMR,
          na.rm = TRUE
        ),

      Gene =
        Gene[
          which.min(
            p_SMR
          )
        ],

      Probe_Chr =
        Probe_Chr_num[
          which.min(
            p_SMR
          )
        ],

      Probe_bp =
        Probe_bp[
          which.min(
            p_SMR
          )
        ]
    ),
    by = .(
      trait,
      ProbeID
    )
  ][
    order(
      trait,
      -n_tissues,
      min_p_SMR
    )
  ]

} else {
  RECURRENCE <- data.table()
}

RECURRENCE_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_primary_sQTL_cross_tissue_recurrence.csv"
)

fwrite(
  RECURRENCE,
  RECURRENCE_FILE
)

# ============================================================
# 18. eQTL + sQTL TARGET-REGION INTEGRATION
# ============================================================

INTEGRATION <- data.table()

if (
  file.exists(EQTL_MASTER)
) {

  EQTL <- fread(
    EQTL_MASTER,
    showProgress = FALSE
  )

  if (
    all(
      c(
        "trait",
        "Probe_Chr",
        "Probe_bp",
        "p_SMR",
        "Primary_pass"
      ) %in%
        names(EQTL)
    )
  ) {

    EQTL[
      ,
      Probe_Chr_num :=
        suppressWarnings(
          as.integer(
            gsub(
              "^chr",
              "",
              as.character(
                Probe_Chr
              ),
              ignore.case = TRUE
            )
          )
        )
    ]

    EQTL_TARGET_LIST <- list()

    for (i in seq_len(nrow(REGIONS))) {

      rr <- REGIONS[i]

      z <- EQTL[
        Probe_Chr_num ==
          rr$chr &
        is.finite(Probe_bp) &
        Probe_bp >=
          rr$start &
        Probe_bp <=
          rr$stop
      ]

      if (nrow(z)) {

        z[
          ,
          target_region :=
            rr$region
        ]

        EQTL_TARGET_LIST[[
          length(EQTL_TARGET_LIST) + 1L
        ]] <- z
      }
    }

    EQTL_TARGET <- if (
      length(EQTL_TARGET_LIST)
    ) {
      rbindlist(
        EQTL_TARGET_LIST,
        fill = TRUE
      )
    } else {
      data.table()
    }

    E_SUM <- if (
      nrow(EQTL_TARGET)
    ) {
      EQTL_TARGET[
        ,
        .(
          eQTL_n_rows = .N,
          eQTL_n_primary =
            sum(
              Primary_pass,
              na.rm = TRUE
            ),
          eQTL_min_p_SMR =
            safe_min(
              p_SMR
            )
        ),
        by = .(
          target_region,
          trait
        )
      ]
    } else {
      data.table()
    }

    S_SUM <- if (
      nrow(TARGET)
    ) {
      TARGET[
        ,
        .(
          sQTL_n_rows = .N,
          sQTL_n_primary =
            sum(
              Primary_pass,
              na.rm = TRUE
            ),
          sQTL_min_p_SMR =
            safe_min(
              p_SMR
            )
        ),
        by = .(
          target_region,
          trait
        )
      ]
    } else {
      data.table()
    }

    INTEGRATION <- merge(
      GRID,
      E_SUM,
      by = c(
        "target_region",
        "trait"
      ),
      all.x = TRUE
    )

    INTEGRATION <- merge(
      INTEGRATION,
      S_SUM,
      by = c(
        "target_region",
        "trait"
      ),
      all.x = TRUE
    )

    INTEGRATION[
      ,
      regulatory_pattern :=
        fifelse(
          !is.na(eQTL_n_primary) &
          eQTL_n_primary > 0 &
          !is.na(sQTL_n_primary) &
          sQTL_n_primary > 0,
          "eQTL_and_sQTL",

          fifelse(
            !is.na(eQTL_n_primary) &
            eQTL_n_primary > 0,
            "eQTL_only",

            fifelse(
              !is.na(sQTL_n_primary) &
              sQTL_n_primary > 0,
              "sQTL_only",
              "no_strict_regulatory_signal"
            )
          )
        )
    ]
  }
}

INTEGRATION_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_eQTL_sQTL_target_region_integration.csv"
)

fwrite(
  INTEGRATION,
  INTEGRATION_FILE
)

# ============================================================
# 19. COMPLETION SUMMARY
# ============================================================

COMPLETION <- data.table(
  item = c(
    "source_chr_split_datasets_parsed",
    "selected_tissues_with_22chr",
    "fullgenome_merges_valid",
    "merged_descriptive_cis_all_pass",
    "formal_jobs_total",
    "formal_jobs_valid",
    "chromosomes_present_in_final_SMR",
    "total_fullgenome_sQTL_SMR_rows",
    "primary_sQTL_rows",
    "primary_AF_rows",
    "primary_HFpEF_rows",
    "primary_BMI_rows",
    "primary_OSA_rows",
    "Bonferroni_HEIDI_heterogeneity_rows",
    "Bonferroni_HEIDI_not_evaluable_rows",
    "strict_cross_trait_shared_splicing_events",
    "cross_trait_shared_sQTL_gene_pairs"
  ),

  value = c(
    nrow(SOURCE),
    sum(CHR_QC$PASS),
    sum(MERGE_QC$merged_valid),
    all(
      MERGE_QC$descriptive_cis_status ==
        "PASS"
    ),
    nrow(RUN_QC),
    sum(
      RUN_QC$valid_smr_output
    ),
    paste(
      CHR_PRESENT,
      collapse = ";"
    ),
    nrow(ALL),
    nrow(PRIMARY),
    sum(
      PRIMARY$trait == "AF"
    ),
    sum(
      PRIMARY$trait == "HFpEF"
    ),
    sum(
      PRIMARY$trait == "BMI"
    ),
    sum(
      PRIMARY$trait == "OSA"
    ),
    nrow(HEIDI_FAIL),
    nrow(HEIDI_NE),
    nrow(SHARED_EVENT),
    nrow(SHARED_GENE)
  )
)

COMPLETION_FILE <- file.path(
  TABLE_DIR,
  "STEP9C2_V3_completion_summary.csv"
)

fwrite(
  COMPLETION,
  COMPLETION_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP9C2_V3_sessionInfo.txt"
  )
)

# ============================================================
# 20. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP9C2_V3 COMPLETE — FULL-GENOME GTEx sQTL SMR/HEIDI\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nCompletion summary:\n"
)

print(
  COMPLETION
)

cat(
  "\nTarget-region summary:\n"
)

print(
  TARGET_SUMMARY
)

cat(
  "\nSTRICT shared splicing events across traits:\n"
)

print(
  SHARED_EVENT
)

cat(
  "\nCross-trait shared sQTL genes:\n"
)

print(
  SHARED_GENE
)

cat(
  "\nUPLOAD THESE 10 FILES:\n",
  "1) ", COMPLETION_FILE, "\n",
  "2) ", CHR_QC_FILE, "\n",
  "3) ", MERGE_QC_FILE, "\n",
  "4) ", RUN_QC_FILE, "\n",
  "5) ", SELECTED_FILE, "\n",
  "6) ", PRIMARY_FILE, "\n",
  "7) ", HEIDI_FAIL_FILE, "\n",
  "8) ", HEIDI_NE_FILE, "\n",
  "9) ", TARGET_SUMMARY_FILE, "\n",
  "10) ", SHARED_EVENT_FILE, "\n",
  sep = ""
)

cat(
  "\nAlso upload:\n",
  SHARED_GENE_FILE,
  "\n",
  INTEGRATION_FILE,
  "\n",
  sep = ""
)

cat(
  "====================================================\n"
)
