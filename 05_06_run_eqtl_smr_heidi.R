# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP9B_V2_RUN_32_SMR_HEIDI_GTEX_LITE_FIXED.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP9B_V2 — FORMAL SMR + HEIDI SCREENING ACROSS 32 GTEx ANALYSES
# Project: AF / HFpEF / BMI / OSA
#
# INPUTS CREATED BY STEP9A
# ------------------------
# - 4 GWAS .ma files
# - 8 GTEx v8 cis-eQTL LITE BESD datasets
# - whole-genome 1000G EUR PLINK reference
# - frozen 4-trait x 8-tissue job plan
#
# ANALYSIS DESIGN
# ---------------
# 4 traits x 8 pre-specified tissues = 32 analyses
#
# Primary SMR significance:
#   p_SMR < 0.05 / number_of_probes_in_that_tissue
#
# HEIDI interpretation:
#   p_HEIDI > 0.01 and nsnp_HEIDI >= 3  -> HEIDI_PASS
#   p_HEIDI <= 0.01                     -> HETEROGENEITY
#   p_HEIDI missing or nsnp_HEIDI < 3   -> NOT_EVALUABLE
#
# IMPORTANT:
# HEIDI_PASS means "no evidence of heterogeneity at the selected
# threshold"; it is NOT proof of mediation or causality.
#
# GTEx v8 LITE is a SCREENING resource. If a key SMR signal has
# HEIDI NOT_EVALUABLE because of insufficient SNP coverage, that
# tissue/gene will be flagged for FULL GTEx follow-up rather than
# rescued by relaxing HEIDI thresholds.
#
# MiXeR output is NOT required for this analysis.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PROJECT PATHS
# ============================================================

ROOT <- "D:/A/data/STEP9_SMR"
DATA_ROOT <- "D:/A/data"

SMR_EXE <- file.path(
  ROOT,
  "00_software",
  "smr-1.3.1-win-x86_64",
  "smr-1.3.1-win.exe"
)

TABLE_DIR <- file.path(
  ROOT,
  "05_tables"
)

RESULT_ROOT <- file.path(
  ROOT,
  "04_smr_results"
)

LOG_DIR <- file.path(
  ROOT,
  "07_logs",
  "STEP9B"
)

STEP9B_TABLE_DIR <- file.path(
  TABLE_DIR,
  "STEP9B"
)

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  STEP9B_TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

JOB_PLAN_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_STEP9B_job_plan_32runs.csv"
)

READINESS_FILE <- file.path(
  TABLE_DIR,
  "STEP9A_readiness.csv"
)

TARGET_REGION_FILE <- file.path(
  DATA_ROOT,
  "STEP8_COLOC",
  "03_tables",
  "STEP8A_target_regions.csv"
)

# ============================================================
# 1. FROZEN ANALYSIS PARAMETERS
# ============================================================

RUN_SMR <- TRUE

# If TRUE, a valid existing .smr result is not re-run.
RESUME_EXISTING <- TRUE

THREAD_NUM <- 8L

MAF_THRESHOLD <- 0.01
PEQTL_SMR <- 5e-8
PEQTL_HEIDI <- 1.57e-3

# New HEIDI method.
HEIDI_METHOD <- 1L

HEIDI_MIN_M <- 3L
HEIDI_MAX_M <- 20L

# Primary HEIDI threshold.
HEIDI_P_THRESHOLD <- 0.01

# Do NOT add --max_num_ld here:
# this option was introduced in SMR v1.4.0, whereas the current
# Windows binary used in this project is v1.3.1.

# Pre-specified priority genes for the chr16/FTO-region biology.
PRIORITY_GENES <- c(
  "FTO",
  "IRX3",
  "IRX5",
  "RPGRIP1L"
)

# Regions to extract separately after genome-wide screening.
TARGET_REGION_NAMES <- c(
  "chr16_locus2135_FTO",
  "chr11_recurrent",
  "chr3_BMI_centered",
  "chr2_AF_HFpEF"
)

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
# 3. GENERIC HELPERS
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

# SMR frequently expects a FILE PREFIX rather than an existing file:
#   --bfile         prefix -> prefix.bed/.bim/.fam
#   --beqtl-summary prefix -> prefix.besd/.epi/.esi
#   --out           prefix -> output files created later
#
# normalizePath(prefix, mustWork=TRUE) is therefore WRONG because the
# bare prefix itself usually does not exist. Normalize only its parent
# directory, then append the basename unchanged.
norm_prefix <- function(
  prefix,
  parent_must_exist = TRUE
) {

  prefix <- as.character(prefix)

  parent <- dirname(prefix)
  stem <- basename(prefix)

  parent_norm <- normalizePath(
    parent,
    winslash = "/",
    mustWork = parent_must_exist
  )

  paste0(
    sub(
      "/$",
      "",
      parent_norm
    ),
    "/",
    stem
  )
}

safe_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_.\\-]+",
    "_",
    x
  )
}

as_num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

compact_text <- function(
  x,
  n = 20L
) {

  x <- trimws(
    as.character(x)
  )

  x <- unique(
    x[
      !is.na(x) &
      nzchar(x)
    ]
  )

  if (!length(x)) {
    return("")
  }

  if (length(x) > n) {
    x <- c(
      x[seq_len(n)],
      paste0(
        "... [",
        length(x) - n,
        " more]"
      )
    )
  }

  paste(
    x,
    collapse = " || "
  )
}

read_tail <- function(
  f,
  n = 40L
) {

  if (!file.exists(f)) {
    return("")
  }

  z <- tryCatch(
    readLines(
      f,
      warn = FALSE
    ),
    error = function(e) character(0)
  )

  compact_text(
    tail(
      z,
      n
    ),
    n = n
  )
}

pick_existing_column <- function(
  nms,
  candidates,
  required = TRUE,
  label = "column"
) {

  hit <- candidates[
    candidates %in%
      nms
  ]

  if (length(hit)) {
    return(
      hit[1]
    )
  }

  # Case-insensitive fallback.
  low <- tolower(nms)

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
      "Could not find ",
      label,
      ". Available columns: ",
      paste(
        nms,
        collapse = ", "
      )
    )
  }

  NA_character_
}

# ============================================================
# 4. HARD INPUT QC
# ============================================================

required_files <- c(
  SMR_EXE,
  JOB_PLAN_FILE,
  READINESS_FILE
)

if (any(!file.exists(required_files))) {
  stop(
    "Missing required STEP9B input(s):\n",
    paste(
      required_files[
        !file.exists(required_files)
      ],
      collapse = "\n"
    )
  )
}

READINESS <- fread(
  READINESS_FILE
)

if (
  !"pass" %in%
  names(READINESS)
) {
  stop(
    "STEP9A_readiness.csv lacks 'pass' column."
  )
}

pass_vec <- READINESS$pass

if (!is.logical(pass_vec)) {
  pass_vec <- toupper(
    as.character(
      pass_vec
    )
  ) %in% c(
    "TRUE",
    "T",
    "1",
    "YES"
  )
}

if (!all(pass_vec)) {
  stop(
    "STEP9A readiness is not fully PASS. ",
    "Do not start formal SMR/HEIDI."
  )
}

JOBS <- fread(
  JOB_PLAN_FILE
)

required_job_cols <- c(
  "trait",
  "tissue",
  "gwas_ma",
  "eqtl_prefix",
  "n_probes",
  "per_tissue_probe_bonferroni",
  "reference_mode",
  "whole_genome_bfile",
  "output_prefix"
)

if (!all(
  required_job_cols %in%
  names(JOBS)
)) {
  stop(
    "STEP9A job plan lacks columns: ",
    paste(
      setdiff(
        required_job_cols,
        names(JOBS)
      ),
      collapse = ", "
    )
  )
}

if (nrow(JOBS) != 32L) {
  stop(
    "Expected 32 frozen SMR jobs but found ",
    nrow(JOBS),
    "."
  )
}

if (
  !setequal(
    unique(JOBS$trait),
    c(
      "AF",
      "HFpEF",
      "BMI",
      "OSA"
    )
  )
) {
  stop(
    "Trait set in the job plan is not AF/HFpEF/BMI/OSA."
  )
}

if (uniqueN(JOBS$tissue) != 8L) {
  stop(
    "Expected 8 pre-specified GTEx tissues."
  )
}

# The successful STEP9A detected a whole-genome 1000G EUR prefix.
if (
  any(
    is.na(
      JOBS$whole_genome_bfile
    ) |
    !nzchar(
      JOBS$whole_genome_bfile
    )
  )
) {
  stop(
    "A whole-genome 1000G EUR --bfile prefix is required for STEP9B."
  )
}

# Verify every referenced input exists.
for (i in seq_len(nrow(JOBS))) {

  gwas <- JOBS$gwas_ma[i]
  eqtl <- JOBS$eqtl_prefix[i]
  bfile <- JOBS$whole_genome_bfile[i]

  needed <- c(
    gwas,
    paste0(
      eqtl,
      ".besd"
    ),
    paste0(
      eqtl,
      ".epi"
    ),
    paste0(
      eqtl,
      ".esi"
    ),
    paste0(
      bfile,
      ".bed"
    ),
    paste0(
      bfile,
      ".bim"
    ),
    paste0(
      bfile,
      ".fam"
    )
  )

  if (any(!file.exists(needed))) {
    stop(
      "Missing files for job ",
      JOBS$trait[i],
      " x ",
      JOBS$tissue[i],
      ":\n",
      paste(
        needed[
          !file.exists(needed)
        ],
        collapse = "\n"
      )
    )
  }
}

# ============================================================
# 4B. PREFIX COMPONENT QC
# ============================================================

check_plink_prefix <- function(prefix) {

  needed <- paste0(
    prefix,
    c(
      ".bed",
      ".bim",
      ".fam"
    )
  )

  all(
    file.exists(
      needed
    )
  )
}

check_besd_prefix <- function(prefix) {

  needed <- paste0(
    prefix,
    c(
      ".besd",
      ".epi",
      ".esi"
    )
  )

  all(
    file.exists(
      needed
    )
  )
}

# Re-confirm every unique prefix before launching any of the 32 jobs.
unique_bfiles <- unique(
  as.character(
    JOBS$whole_genome_bfile
  )
)

unique_eqtls <- unique(
  as.character(
    JOBS$eqtl_prefix
  )
)

bad_bfiles <- unique_bfiles[
  !vapply(
    unique_bfiles,
    check_plink_prefix,
    logical(1)
  )
]

bad_eqtls <- unique_eqtls[
  !vapply(
    unique_eqtls,
    check_besd_prefix,
    logical(1)
  )
]

if (length(bad_bfiles)) {
  stop(
    "Invalid PLINK --bfile prefix(es):\n",
    paste(
      bad_bfiles,
      collapse = "\n"
    )
  )
}

if (length(bad_eqtls)) {
  stop(
    "Invalid BESD --beqtl-summary prefix(es):\n",
    paste(
      bad_eqtls,
      collapse = "\n"
    )
  )
}

cat(
  "\nPrefix QC PASS:\n",
  "  PLINK reference prefix(es): ",
  length(unique_bfiles),
  "\n",
  "  GTEx BESD prefix(es): ",
  length(unique_eqtls),
  "\n",
  sep = ""
)

# ============================================================
# 5. SMR COMMAND BUILDER
# ============================================================

build_smr_args <- function(
  gwas,
  eqtl,
  bfile,
  outprefix
) {

  # IMPORTANT:
  # bfile and eqtl are prefixes, not standalone files.
  # Their component files have already been hard-checked above.
  bfile_prefix <- norm_prefix(
    bfile,
    parent_must_exist = TRUE
  )

  eqtl_prefix <- norm_prefix(
    eqtl,
    parent_must_exist = TRUE
  )

  out_prefix <- norm_prefix(
    outprefix,
    parent_must_exist = TRUE
  )

  c(
    "--bfile",
    shQuote(
      bfile_prefix
    ),

    "--gwas-summary",
    shQuote(
      norm_path(
        gwas,
        TRUE
      )
    ),

    "--beqtl-summary",
    shQuote(
      eqtl_prefix
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
      out_prefix
    )
  )
}

# ============================================================
# 6. VALIDATE EXISTING .smr OUTPUT
# ============================================================

smr_file_is_valid <- function(f) {

  if (
    !file.exists(f) ||
    file.info(f)$size <= 0
  ) {
    return(FALSE)
  }

  hdr <- tryCatch(
    fread(
      f,
      nrows = 2L,
      na.strings = c(
        "NA",
        "nan",
        "NaN"
      ),
      showProgress = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(hdr)) {
    return(FALSE)
  }

  nms <- names(hdr)

  has_probe <- any(
    tolower(nms) ==
      "probeid"
  )

  has_smr <- any(
    tolower(nms) ==
      "p_smr"
  )

  has_probe &&
    has_smr
}

# ============================================================
# 7. RUN ALL 32 SMR/HEIDI JOBS
# ============================================================

RUN_QC_LIST <- list()
COMMAND_LIST <- list()

cat(
  "\n====================================================\n",
  "STEP9B_V2 — STARTING 32 FORMAL SMR/HEIDI JOBS\n",
  "====================================================\n",
  sep = ""
)

for (i in seq_len(nrow(JOBS))) {

  tr <- as.character(
    JOBS$trait[i]
  )

  tissue <- as.character(
    JOBS$tissue[i]
  )

  gwas <- as.character(
    JOBS$gwas_ma[i]
  )

  eqtl <- as.character(
    JOBS$eqtl_prefix[i]
  )

  bfile <- as.character(
    JOBS$whole_genome_bfile[i]
  )

  outprefix <- as.character(
    JOBS$output_prefix[i]
  )

  dir.create(
    dirname(
      outprefix
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )

  smr_file <- paste0(
    outprefix,
    ".smr"
  )

  logfile <- file.path(
    LOG_DIR,
    paste0(
      sprintf(
        "%02d",
        i
      ),
      "__",
      safe_filename(tr),
      "__",
      safe_filename(tissue),
      ".log"
    )
  )

  args <- build_smr_args(
    gwas,
    eqtl,
    bfile,
    outprefix
  )

  command_string <- paste(
    shQuote(
      SMR_EXE
    ),
    paste(
      args,
      collapse = " "
    )
  )

  COMMAND_LIST[[
    length(COMMAND_LIST) + 1L
  ]] <- data.table(
    job_index = i,
    trait = tr,
    tissue = tissue,
    command = command_string,
    output_prefix = outprefix
  )

  cat(
    "\n[",
    i,
    "/",
    nrow(JOBS),
    "] ",
    tr,
    " x ",
    tissue,
    "\n",
    sep = ""
  )

  already_valid <- smr_file_is_valid(
    smr_file
  )

  start_time <- Sys.time()

  if (
    RESUME_EXISTING &&
    already_valid
  ) {

    status <- 0L
    run_state <- "SKIPPED_EXISTING_VALID"

    cat(
      "  Existing valid .smr found — skipping.\n"
    )

  } else if (!RUN_SMR) {

    status <- NA_integer_
    run_state <- "DRY_RUN"

  } else {

    # Remove stale partial output before rerun.
    if (file.exists(smr_file)) {
      unlink(
        smr_file
      )
    }

    status <- suppressWarnings(
      system2(
        SMR_EXE,
        args = args,
        stdout = logfile,
        stderr = logfile
      )
    )

    status <- as.integer(
      status
    )

    if (
      identical(
        status,
        0L
      ) &&
      smr_file_is_valid(
        smr_file
      )
    ) {

      run_state <- "PASS"

    } else {

      run_state <- "FAIL"
    }
  }

  end_time <- Sys.time()

  result_valid <- smr_file_is_valid(
    smr_file
  )

  n_result_rows <- NA_integer_

  if (result_valid) {

    n_result_rows <- tryCatch({

      # Count data lines without retaining the whole result.
      con_count <- file(
        smr_file,
        open = "rt"
      )

      on.exit(
        close(con_count),
        add = TRUE
      )

      nlines <- 0L

      repeat {

        z <- readLines(
          con_count,
          n = 100000L,
          warn = FALSE
        )

        if (!length(z)) break

        nlines <- nlines +
          length(z)
      }

      close(
        con_count
      )

      max(
        nlines - 1L,
        0L
      )

    }, error = function(e) {
      NA_integer_
    })
  }

  RUN_QC_LIST[[
    length(RUN_QC_LIST) + 1L
  ]] <- data.table(
    job_index = i,
    trait = tr,
    tissue = tissue,
    status_code = status,
    run_state = run_state,
    result_valid = result_valid,
    n_result_rows = n_result_rows,
    seconds = as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    smr_file = smr_file,
    log_file = logfile,
    log_tail = if (
      run_state == "FAIL"
    ) {
      read_tail(
        logfile,
        50L
      )
    } else {
      ""
    }
  )

  if (run_state == "FAIL") {

    cat(
      "  FAIL — continuing to next job. See log:\n  ",
      logfile,
      "\n",
      sep = ""
    )

  } else {

    cat(
      "  ",
      run_state,
      " | rows=",
      n_result_rows,
      "\n",
      sep = ""
    )
  }
}

RUN_QC <- rbindlist(
  RUN_QC_LIST,
  fill = TRUE
)

COMMANDS <- rbindlist(
  COMMAND_LIST,
  fill = TRUE
)

RUN_QC_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_run_QC.csv"
)

COMMAND_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_exact_commands.csv"
)

fwrite(
  RUN_QC,
  RUN_QC_FILE
)

fwrite(
  COMMANDS,
  COMMAND_FILE
)

# ============================================================
# 8. STOP RESULT INTERPRETATION IF ANY JOB FAILED
#
# We save QC first so failures can be inspected without losing
# completed jobs. Re-running the same script resumes valid outputs.
# ============================================================

FAILED <- RUN_QC[
  run_state ==
    "FAIL" |
  result_valid ==
    FALSE
]

if (nrow(FAILED)) {

  cat(
    "\n====================================================\n",
    "STEP9B RUN PHASE INCOMPLETE\n",
    "====================================================\n",
    "Failed/invalid jobs: ",
    nrow(FAILED),
    "\n",
    sep = ""
  )

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
    "At least one SMR/HEIDI job failed. ",
    "Upload STEP9B_run_QC.csv before changing any scientific parameter. ",
    "Valid completed jobs will be skipped on the next run."
  )
}

# ============================================================
# 9. READ + NORMALIZE ALL 32 .smr FILES
# ============================================================

RESULT_LIST <- list()

for (i in seq_len(nrow(JOBS))) {

  tr <- as.character(
    JOBS$trait[i]
  )

  tissue <- as.character(
    JOBS$tissue[i]
  )

  smr_file <- paste0(
    JOBS$output_prefix[i],
    ".smr"
  )

  x <- fread(
    smr_file,
    na.strings = c(
      "NA",
      "nan",
      "NaN",
      "-9"
    ),
    showProgress = FALSE
  )

  if (!nrow(x)) {
    next
  }

  nms <- names(x)

  col_probe <- pick_existing_column(
    nms,
    c(
      "ProbeID"
    ),
    TRUE,
    "ProbeID"
  )

  col_probe_chr <- pick_existing_column(
    nms,
    c(
      "Probe_Chr",
      "ProbeChr"
    ),
    FALSE,
    "Probe chromosome"
  )

  col_gene <- pick_existing_column(
    nms,
    c(
      "Gene",
      "GeneName"
    ),
    FALSE,
    "Gene"
  )

  col_probe_bp <- pick_existing_column(
    nms,
    c(
      "Probe_bp",
      "ProbeBP"
    ),
    FALSE,
    "Probe position"
  )

  col_snp <- pick_existing_column(
    nms,
    c(
      "SNP"
    ),
    FALSE,
    "top eQTL SNP"
  )

  col_b_smr <- pick_existing_column(
    nms,
    c(
      "b_SMR"
    ),
    FALSE
  )

  col_se_smr <- pick_existing_column(
    nms,
    c(
      "se_SMR"
    ),
    FALSE
  )

  col_p_smr <- pick_existing_column(
    nms,
    c(
      "p_SMR"
    ),
    TRUE,
    "p_SMR"
  )

  col_p_heidi <- pick_existing_column(
    nms,
    c(
      "p_HEIDI"
    ),
    FALSE
  )

  col_nsnp_heidi <- pick_existing_column(
    nms,
    c(
      "nsnp_HEIDI"
    ),
    FALSE
  )

  col_p_gwas <- pick_existing_column(
    nms,
    c(
      "p_GWAS"
    ),
    FALSE
  )

  col_p_eqtl <- pick_existing_column(
    nms,
    c(
      "p_eQTL"
    ),
    FALSE
  )

  out <- data.table(
    trait = tr,
    tissue = tissue,
    ProbeID = as.character(
      x[[col_probe]]
    ),

    Probe_Chr = if (
      !is.na(col_probe_chr)
    ) {
      as.character(
        x[[col_probe_chr]]
      )
    } else {
      NA_character_
    },

    Gene = if (
      !is.na(col_gene)
    ) {
      as.character(
        x[[col_gene]]
      )
    } else {
      NA_character_
    },

    Probe_bp = if (
      !is.na(col_probe_bp)
    ) {
      as_num(
        x[[col_probe_bp]]
      )
    } else {
      NA_real_
    },

    top_eQTL_SNP = if (
      !is.na(col_snp)
    ) {
      as.character(
        x[[col_snp]]
      )
    } else {
      NA_character_
    },

    b_SMR = if (
      !is.na(col_b_smr)
    ) {
      as_num(
        x[[col_b_smr]]
      )
    } else {
      NA_real_
    },

    se_SMR = if (
      !is.na(col_se_smr)
    ) {
      as_num(
        x[[col_se_smr]]
      )
    } else {
      NA_real_
    },

    p_SMR = as_num(
      x[[col_p_smr]]
    ),

    p_HEIDI = if (
      !is.na(col_p_heidi)
    ) {
      as_num(
        x[[col_p_heidi]]
      )
    } else {
      NA_real_
    },

    nsnp_HEIDI = if (
      !is.na(col_nsnp_heidi)
    ) {
      as_num(
        x[[col_nsnp_heidi]]
      )
    } else {
      NA_real_
    },

    p_GWAS = if (
      !is.na(col_p_gwas)
    ) {
      as_num(
        x[[col_p_gwas]]
      )
    } else {
      NA_real_
    },

    p_eQTL = if (
      !is.na(col_p_eqtl)
    ) {
      as_num(
        x[[col_p_eqtl]]
      )
    } else {
      NA_real_
    },

    n_probes_tissue = as.integer(
      JOBS$n_probes[i]
    ),

    tissue_Bonferroni = as.numeric(
      JOBS$per_tissue_probe_bonferroni[i]
    ),

    source_smr_file = norm_path(
      smr_file,
      TRUE
    )
  )

  RESULT_LIST[[
    length(RESULT_LIST) + 1L
  ]] <- out
}

ALL <- rbindlist(
  RESULT_LIST,
  fill = TRUE
)

if (!nrow(ALL)) {
  stop(
    "All 32 jobs completed but no SMR result rows were read."
  )
}

# Remove invalid P rows only for inferential summaries.
ALL[
  !is.finite(p_SMR) |
  p_SMR <= 0 |
  p_SMR > 1,
  p_SMR := NA_real_
]

# ============================================================
# 10. MULTIPLE-TESTING + HEIDI CLASSIFICATION
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
    is.finite(
      p_SMR
    ) &
    p_SMR <
      tissue_Bonferroni
]

ALL[
  ,
  HEIDI_evaluable :=
    is.finite(
      p_HEIDI
    ) &
    is.finite(
      nsnp_HEIDI
    ) &
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
      "SMR_Bonferroni_and_HEIDI_pass",
      fifelse(
        SMR_Bonferroni_pass &
        HEIDI_status ==
          "FAIL_HETEROGENEITY",
        "SMR_Bonferroni_but_HEIDI_heterogeneity",
        fifelse(
          SMR_Bonferroni_pass &
          HEIDI_status ==
            "NOT_EVALUABLE",
          "SMR_Bonferroni_HEIDI_not_evaluable",
          fifelse(
            is.finite(
              FDR_tissue
            ) &
            FDR_tissue <
              0.05,
            "Tissue_FDR_only",
            "No_primary_SMR_evidence"
          )
        )
      )
    )
]

# Clean chromosome field for region matching.
ALL[
  ,
  Probe_Chr_num :=
    suppressWarnings(
      as.integer(
        gsub(
          "^chr",
          "",
          Probe_Chr,
          ignore.case = TRUE
        )
      )
    )
]

# Gene key for cross-tissue recurrence:
# use Gene when available, otherwise ProbeID.
ALL[
  ,
  gene_key :=
    fifelse(
      !is.na(Gene) &
      nzchar(
        trimws(Gene)
      ),
      toupper(
        trimws(Gene)
      ),
      toupper(
        trimws(ProbeID)
      )
    )
]

# ============================================================
# 11. SAVE MASTER RESULT TABLE
# ============================================================

MASTER_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_all_32_SMR_HEIDI_results.csv.gz"
)

fwrite(
  ALL,
  MASTER_FILE
)

# ============================================================
# 12. JOB-LEVEL RESULT QC
# ============================================================

JOB_RESULT_QC <- ALL[
  ,
  .(
    n_result_rows = .N,

    n_SMR_Bonferroni =
      sum(
        SMR_Bonferroni_pass,
        na.rm = TRUE
      ),

    n_HEIDI_evaluable =
      sum(
        HEIDI_evaluable,
        na.rm = TRUE
      ),

    n_primary_pass =
      sum(
        Primary_pass,
        na.rm = TRUE
      ),

    n_Bonferroni_HEIDI_not_evaluable =
      sum(
        SMR_Bonferroni_pass &
        HEIDI_status ==
          "NOT_EVALUABLE",
        na.rm = TRUE
      ),

    n_Bonferroni_HEIDI_fail =
      sum(
        SMR_Bonferroni_pass &
        HEIDI_status ==
          "FAIL_HETEROGENEITY",
        na.rm = TRUE
      ),

    min_p_SMR =
      min(
        p_SMR,
        na.rm = TRUE
      )
  ),
  by = .(
    trait,
    tissue
  )
]

JOB_RESULT_QC[
  !is.finite(
    min_p_SMR
  ),
  min_p_SMR := NA_real_
]

JOB_RESULT_QC_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_job_result_summary.csv"
)

fwrite(
  JOB_RESULT_QC,
  JOB_RESULT_QC_FILE
)

# ============================================================
# 13. PRIMARY + SECONDARY HIT TABLES
# ============================================================

PRIMARY_HITS <- ALL[
  Primary_pass == TRUE
][
  order(
    trait,
    p_SMR
  )
]

PRIMARY_HITS_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_primary_Bonferroni_HEIDI_pass.csv"
)

fwrite(
  PRIMARY_HITS,
  PRIMARY_HITS_FILE
)

BONF_HEIDI_FAIL <- ALL[
  SMR_Bonferroni_pass == TRUE &
  HEIDI_status ==
    "FAIL_HETEROGENEITY"
][
  order(
    trait,
    p_SMR
  )
]

BONF_HEIDI_FAIL_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_Bonferroni_but_HEIDI_heterogeneity.csv"
)

fwrite(
  BONF_HEIDI_FAIL,
  BONF_HEIDI_FAIL_FILE
)

FULL_GTEX_RETEST <- ALL[
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

FULL_GTEX_RETEST[
  ,
  retest_reason :=
    "SMR Bonferroni significant but GTEx LITE HEIDI not evaluable; use corresponding FULL tissue BESD"
]

FULL_GTEX_RETEST_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_FULL_GTEx_retest_candidates.csv"
)

fwrite(
  FULL_GTEX_RETEST,
  FULL_GTEX_RETEST_FILE
)

# ============================================================
# 14. CROSS-TISSUE RECURRENCE
#
# This is descriptive replication across pre-specified tissues.
# It does not convert correlated tissues into independent studies.
# ============================================================

RECURRENCE_SOURCE <- unique(
  PRIMARY_HITS[
    ,
    .(
      trait,
      tissue,
      gene_key,
      Gene,
      ProbeID,
      p_SMR,
      p_HEIDI,
      top_eQTL_SNP
    )
  ]
)

if (nrow(RECURRENCE_SOURCE)) {

  GENE_RECURRENCE <- RECURRENCE_SOURCE[
    ,
    .(
      n_tissues_primary =
        uniqueN(
          tissue
        ),

      tissues_primary =
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

      max_p_HEIDI =
        max(
          p_HEIDI,
          na.rm = TRUE
        ),

      probes =
        paste(
          sort(
            unique(
              ProbeID
            )
          ),
          collapse = ";"
        )
    ),
    by = .(
      trait,
      gene_key
    )
  ][
    order(
      trait,
      -n_tissues_primary,
      min_p_SMR
    )
  ]

} else {

  GENE_RECURRENCE <- data.table(
    trait = character(0),
    gene_key = character(0),
    n_tissues_primary = integer(0),
    tissues_primary = character(0),
    min_p_SMR = numeric(0),
    max_p_HEIDI = numeric(0),
    probes = character(0)
  )
}

GENE_RECURRENCE_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_cross_tissue_primary_recurrence.csv"
)

fwrite(
  GENE_RECURRENCE,
  GENE_RECURRENCE_FILE
)

# ============================================================
# 15. TARGET REGION EXTRACTION
# ============================================================

TARGET_REGION_RESULTS <- data.table()

if (file.exists(TARGET_REGION_FILE)) {

  REGIONS <- fread(
    TARGET_REGION_FILE
  )

  required_region_cols <- c(
    "region",
    "chr",
    "start",
    "stop"
  )

  if (all(
    required_region_cols %in%
    names(REGIONS)
  )) {

    REGIONS[
      ,
      region := as.character(
        region
      )
    ]

    REGIONS <- REGIONS[
      region %in%
        TARGET_REGION_NAMES
    ]

    target_rows <- list()

    for (i in seq_len(nrow(REGIONS))) {

      rr <- REGIONS[i]

      z <- ALL[
        Probe_Chr_num ==
          as.integer(
            rr$chr
          ) &
        is.finite(
          Probe_bp
        ) &
        Probe_bp >=
          as.numeric(
            rr$start
          ) &
        Probe_bp <=
          as.numeric(
            rr$stop
          )
      ]

      if (nrow(z)) {

        z[
          ,
          target_region :=
            as.character(
              rr$region
            )
        ]

        target_rows[[
          length(target_rows) + 1L
        ]] <- z
      }
    }

    if (length(target_rows)) {

      TARGET_REGION_RESULTS <- rbindlist(
        target_rows,
        fill = TRUE
      )[
        order(
          target_region,
          trait,
          p_SMR
        )
      ]
    }
  }
}

TARGET_REGION_FILE_OUT <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_target_region_results_chr16_chr11_chr3_chr2.csv"
)

fwrite(
  TARGET_REGION_RESULTS,
  TARGET_REGION_FILE_OUT
)

# ============================================================
# 16. PRE-SPECIFIED PRIORITY GENE EXTRACTION
# ============================================================

PRIORITY_GENE_RESULTS <- ALL[
  toupper(
    trimws(
      Gene
    )
  ) %in%
    PRIORITY_GENES
][
  order(
    Gene,
    trait,
    p_SMR
  )
]

PRIORITY_GENE_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_priority_genes_FTO_IRX3_IRX5_RPGRIP1L.csv"
)

fwrite(
  PRIORITY_GENE_RESULTS,
  PRIORITY_GENE_FILE
)

# ============================================================
# 17. FULL GTEx TISSUE DOWNLOAD/RETEST PLAN
#
# Summarise only scientifically motivated re-tests:
# SMR Bonferroni significant but LITE HEIDI was not evaluable.
# ============================================================

if (nrow(FULL_GTEX_RETEST)) {

  FULL_TISSUE_PLAN <- FULL_GTEX_RETEST[
    ,
    .(
      n_candidate_rows = .N,
      genes_or_probes =
        paste(
          sort(
            unique(
              gene_key
            )
          ),
          collapse = ";"
        )
    ),
    by = tissue
  ][
    order(
      tissue
    )
  ]

} else {

  FULL_TISSUE_PLAN <- data.table(
    tissue = character(0),
    n_candidate_rows = integer(0),
    genes_or_probes = character(0)
  )
}

FULL_TISSUE_PLAN_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_FULL_GTEx_tissue_download_plan.csv"
)

fwrite(
  FULL_TISSUE_PLAN,
  FULL_TISSUE_PLAN_FILE
)

# ============================================================
# 18. METHODS / PARAMETER PROVENANCE
# ============================================================

METHOD <- data.table(
  parameter = c(
    "SMR_executable",
    "QTL_resource",
    "traits",
    "primary_tissues",
    "n_jobs",
    "LD_reference",
    "MAF_threshold",
    "peqtl_smr",
    "peqtl_heidi",
    "heidi_method",
    "heidi_min_m",
    "heidi_max_m",
    "HEIDI_pass_threshold",
    "primary_SMR_multiple_testing",
    "secondary_multiple_testing",
    "GTEx_LITE_role",
    "MiXeR_dependency"
  ),

  value = c(
    norm_path(
      SMR_EXE,
      TRUE
    ),

    "GTEx v8 cis-eQTL LITE",

    "AF;HFpEF;BMI;OSA",

    paste(
      sort(
        unique(
          JOBS$tissue
        )
      ),
      collapse = ";"
    ),

    nrow(JOBS),

    as.character(
      JOBS$whole_genome_bfile[1]
    ),

    MAF_THRESHOLD,
    PEQTL_SMR,
    PEQTL_HEIDI,
    HEIDI_METHOD,
    HEIDI_MIN_M,
    HEIDI_MAX_M,
    HEIDI_P_THRESHOLD,

    "Per-tissue Bonferroni: 0.05 / number of tested probes",

    "BH-FDR within trait-tissue, across all tissues per trait, and globally",

    "Screening; key SMR hits with non-evaluable HEIDI are re-tested in corresponding FULL GTEx tissue",

    "None"
  )
)

METHOD_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_method_provenance.csv"
)

fwrite(
  METHOD,
  METHOD_FILE
)

# ============================================================
# 19. FINAL READINESS / COMPLETION SUMMARY
# ============================================================

COMPLETION <- data.table(
  item = c(
    "SMR_jobs_total",
    "SMR_jobs_valid",
    "all_jobs_passed",
    "total_SMR_result_rows",
    "primary_Bonferroni_HEIDI_pass_rows",
    "Bonferroni_HEIDI_heterogeneity_rows",
    "Bonferroni_HEIDI_not_evaluable_rows",
    "FULL_GTEx_tissues_flagged_for_retest"
  ),

  value = c(
    nrow(RUN_QC),
    sum(
      RUN_QC$result_valid
    ),
    all(
      RUN_QC$result_valid
    ),
    nrow(ALL),
    nrow(PRIMARY_HITS),
    nrow(BONF_HEIDI_FAIL),
    nrow(FULL_GTEX_RETEST),
    nrow(FULL_TISSUE_PLAN)
  )
)

COMPLETION_FILE <- file.path(
  STEP9B_TABLE_DIR,
  "STEP9B_completion_summary.csv"
)

fwrite(
  COMPLETION,
  COMPLETION_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP9B_sessionInfo.txt"
  )
)

# ============================================================
# 20. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP9B_V2 COMPLETE — FORMAL SMR + HEIDI SCREENING\n",
  "====================================================\n",
  sep = ""
)

print(
  COMPLETION
)

cat(
  "\nPer-job summary:\n"
)

print(
  JOB_RESULT_QC
)

cat(
  "\nPRIMARY INTERPRETATION RULE:\n",
  "  p_SMR < tissue Bonferroni threshold\n",
  "  AND p_HEIDI > ",
  HEIDI_P_THRESHOLD,
  "\n",
  "  AND nsnp_HEIDI >= ",
  HEIDI_MIN_M,
  "\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n",
  "HEIDI PASS = no detected heterogeneity; it does NOT prove ",
  "causal mediation.\n",
  sep = ""
)

cat(
  "\nUPLOAD THESE 8 FILES:\n",
  "1) ", COMPLETION_FILE, "\n",
  "2) ", RUN_QC_FILE, "\n",
  "3) ", JOB_RESULT_QC_FILE, "\n",
  "4) ", PRIMARY_HITS_FILE, "\n",
  "5) ", BONF_HEIDI_FAIL_FILE, "\n",
  "6) ", FULL_GTEX_RETEST_FILE, "\n",
  "7) ", TARGET_REGION_FILE_OUT, "\n",
  "8) ", PRIORITY_GENE_FILE, "\n",
  sep = ""
)

cat(
  "\nIf file 6 contains rows, also upload:\n",
  FULL_TISSUE_PLAN_FILE,
  "\n",
  sep = ""
)

cat(
  "\nYou do NOT need to upload the 32 individual .smr files unless ",
  "a job-level discrepancy needs investigation.\n",
  sep = ""
)

cat(
  "====================================================\n"
)
