# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G2B_AF_SUSIE_LOWRANK_FINEMAP_TO_PIP(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G2B — AF genome-wide SuSiE-RSS fine-mapping for SCAVENGE
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# INPUT (frozen from STEP11G2A):
#   - 515 AF-selected Berisa-Pickrell EUR LDetect blocks (GRCh37/hg19)
#   - Primary AF GWAS: GCST90624412 EUR (source coordinates GRCh38)
#   - 1000 Genomes Phase 3 EUR reference
#
# PRIMARY FINE-MAPPING MODEL:
#   susieR::susie_rss()
#   L = 1
#   uniform prior (prior_weights = NULL)
#   external 1000G EUR reference
#   estimate_residual_variance = FALSE
#   R_mismatch = "none" (when supported)
#
# COMPUTATIONAL STRATEGY:
#   - DO NOT construct giant dense p x p LD matrices.
#   - PLINK exports the matched reference genotypes as A-transpose.
#   - Current susieR accepts reference genotype factor matrix X.
#   - When B < p, susieR uses a low-rank path.
#
# SENSITIVITY (if current susieR supports it):
#   R_finite = TRUE
#   R_mismatch = "eb"
#   This is saved separately and DOES NOT replace the primary SCAVENGE PIP.
#
# ROBUSTNESS:
#   - one-time AF GWAS chromosome splitting (memory-safe Python streaming)
#   - chromosome-by-chromosome processing
#   - block-level checkpoint / resume
#   - failed blocks do not erase completed blocks
#   - every block gets its own harmonized table, PIP, summary, and RDS
#
# IMPORTANT:
#   SCAVENGE PRIMARY INPUT = PIP_primary from the published-style L=1,
#   uniform-prior fine-mapping.
# ==============================================================================

rm(list = ls())
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = max(3600, getOption("timeout"))
)

# ==============================================================================
# 0. FROZEN SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

L_PRIMARY <- 1L
COVERAGE <- 0.95
MIN_ABS_CORR <- 0.5

MIN_MATCHED_SNPS <- 20L

# Keep the project-wide extended MHC exclusion.
MHC_CHR <- 6L
MHC_START <- 25000000L
MHC_END <- 34000000L

# PLINK resource use: deliberately conservative for a desktop machine.
PLINK_THREADS <- 4L
PLINK_MEMORY_MB <- 4096L

# Retain only compact final objects; large temporary .traw files are removed.
KEEP_TRAW <- FALSE

# ==============================================================================
# 1. INPUT / OUTPUT PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

G1_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "01_AF_FINEMAP_PREP"
)

AF_STANDARD_FILE <- file.path(
  G1_ROOT,
  "01_STANDARDIZED_GWAS",
  "AF_GCST90624412_clean_for_finemap.tsv.gz"
)

G2A_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "02_AF_LDETECT_BLOCKS"
)

G2A_READINESS <- file.path(
  G2A_ROOT,
  "01_QC",
  "STEP11G2A_readiness.csv"
)

BLOCK_FILE <- file.path(
  G2A_ROOT,
  "02_SELECTED_BLOCKS",
  "STEP11G2A_AF_selected_LDetect_blocks_hg19.csv"
)

CHAIN_QC_FILE <- file.path(
  G2A_ROOT,
  "01_QC",
  "STEP11G2A_V4_chain_QC.csv"
)

OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "03_AF_SUSIE_FINEMAP"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

AF_SPLIT_DIR <- file.path(
  OUT_ROOT,
  "01_AF_BY_CHR"
)

BLOCK_OUT_DIR <- file.path(
  OUT_ROOT,
  "02_BLOCK_RESULTS"
)

FINAL_DIR <- file.path(
  OUT_ROOT,
  "03_FINAL_PIP"
)

LOG_DIR <- file.path(
  OUT_ROOT,
  "04_LOGS"
)

TMP_DIR <- file.path(
  OUT_ROOT,
  "05_TMP"
)

for (d in c(
  OUT_ROOT,
  QC_DIR,
  AF_SPLIT_DIR,
  BLOCK_OUT_DIR,
  FINAL_DIR,
  LOG_DIR,
  TMP_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

STATUS_FILE <- file.path(
  QC_DIR,
  "STEP11G2B_block_status.csv"
)

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

cran_pkgs <- c(
  "data.table",
  "R.utils",
  "susieR"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

# Ensure susieR supports the low-rank X interface.
if (
  requireNamespace("susieR", quietly = TRUE) &&
  !"X" %in% names(
    formals(
      susieR::susie_rss
    )
  )
) {

  try(
    unloadNamespace("susieR"),
    silent = TRUE
  )

  install.packages(
    "susieR",
    repos = "https://cloud.r-project.org"
  )
}

if (
  !requireNamespace("susieR", quietly = TRUE) ||
  !"X" %in% names(
    formals(
      susieR::susie_rss
    )
  )
) {
  stop(
    paste0(
      "The installed susieR does not support the low-rank X interface.\n",
      "Please update susieR from CRAN and rerun."
    )
  )
}

bioc_pkgs <- c(
  "GenomicRanges",
  "IRanges",
  "rtracklayer"
)

if (any(
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }

  missing_bioc <- bioc_pkgs[
    !vapply(
      bioc_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

library(data.table)

suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(IRanges)
)
suppressPackageStartupMessages(
  library(rtracklayer)
)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

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

allele_key <- function(a1, a2) {

  a1 <- toupper(
    as.character(a1)
  )

  a2 <- toupper(
    as.character(a2)
  )

  lo <- ifelse(
    a1 <= a2,
    a1,
    a2
  )

  hi <- ifelse(
    a1 <= a2,
    a2,
    a1
  )

  paste0(
    lo,
    "/",
    hi
  )
}

valid_base_allele <- function(x) {
  !is.na(x) &
    toupper(x) %chin% c(
      "A",
      "C",
      "G",
      "T"
    )
}

prepare_plain_chain <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Chain file not found: ",
      path
    )
  }

  path <- norm_path(
    path,
    TRUE
  )

  if (!grepl(
    "\\.gz$",
    path,
    ignore.case = TRUE
  )) {
    return(path)
  }

  plain <- sub(
    "\\.gz$",
    "",
    path,
    ignore.case = TRUE
  )

  if (
    !file.exists(plain) ||
    file.info(plain)$size < 10000
  ) {

    cat(
      "Decompressing chain:\n",
      path,
      "\n-> ",
      plain,
      "\n",
      sep = ""
    )

    cin <- gzfile(
      path,
      "rb"
    )

    cout <- file(
      plain,
      "wb"
    )

    on.exit(
      {
        try(close(cin), silent = TRUE)
        try(close(cout), silent = TRUE)
      },
      add = TRUE
    )

    repeat {

      z <- readBin(
        cin,
        what = "raw",
        n = 1024 * 1024
      )

      if (!length(z)) {
        break
      }

      writeBin(
        z,
        cout
      )
    }

    close(cin)
    close(cout)
  }

  norm_path(
    plain,
    TRUE
  )
}

status_update <- function(row_dt) {

  old <- if (file.exists(STATUS_FILE)) {
    fread(
      STATUS_FILE
    )
  } else {
    data.table()
  }

  if (nrow(old)) {
    old <- old[
      selected_block_label !=
        row_dt$selected_block_label
    ]
  }

  out <- rbindlist(
    list(
      old,
      row_dt
    ),
    fill = TRUE
  )

  if (
    "selected_order" %in%
      names(out)
  ) {
    setorder(
      out,
      selected_order
    )
  }

  fwrite(
    out,
    STATUS_FILE
  )

  invisible(out)
}

safe_scalar <- function(
  x,
  default = NA_real_
) {

  if (is.null(x)) {
    return(default)
  }

  y <- suppressWarnings(
    as.numeric(
      unlist(x)
    )
  )

  y <- y[
    is.finite(y)
  ]

  if (!length(y)) {
    return(default)
  }

  y[1]
}

# ==============================================================================
# 4. VERIFY STEP11G2A
# ==============================================================================

required_inputs <- c(
  G2A_READINESS,
  BLOCK_FILE,
  CHAIN_QC_FILE,
  AF_STANDARD_FILE
)

if (!all(
  file.exists(required_inputs)
)) {
  stop(
    paste0(
      "Missing required input(s):\n",
      paste(
        required_inputs[
          !file.exists(required_inputs)
        ],
        collapse = "\n"
      )
    )
  )
}

READY_A <- fread(
  G2A_READINESS
)

if (
  !all(
    READY_A$pass %in%
      c(
        TRUE,
        "TRUE",
        1
      )
  )
) {

  print(
    READY_A
  )

  stop(
    "STEP11G2A is not fully PASS."
  )
}

BLOCKS <- fread(
  BLOCK_FILE
)

required_block_cols <- c(
  "selected_block_label",
  "block_id",
  "selected_order",
  "CHR",
  "start_hg19",
  "end_hg19",
  "n_gws_variants"
)

missing_block_cols <- setdiff(
  required_block_cols,
  names(BLOCKS)
)

if (length(missing_block_cols)) {
  stop(
    "Selected block manifest is missing: ",
    paste(
      missing_block_cols,
      collapse = ", "
    )
  )
}

setorder(
  BLOCKS,
  selected_order
)

cat(
  "\nSTEP11G2A selected blocks: ",
  nrow(BLOCKS),
  "\n",
  sep = ""
)

# ==============================================================================
# 5. IMPORT VALIDATED hg19 -> hg38 CHAIN
# ==============================================================================

CHAIN_QC <- fread(
  CHAIN_QC_FILE
)

CHAIN_IMPORT_FILE <- CHAIN_QC[
  field == "chain_import_file",
  value
]

if (
  length(CHAIN_IMPORT_FILE) != 1L ||
  !file.exists(CHAIN_IMPORT_FILE)
) {

  chain_candidates <- list.files(
    file.path(
      G2A_ROOT,
      "00_REFERENCE"
    ),
    pattern = "^hg19ToHg38[.]over[.]chain([.]gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(chain_candidates)) {
    stop(
      "Could not recover hg19ToHg38 chain from STEP11G2A."
    )
  }

  CHAIN_IMPORT_FILE <- prepare_plain_chain(
    chain_candidates[1]
  )
}

CHAIN_IMPORT_FILE <- prepare_plain_chain(
  CHAIN_IMPORT_FILE
)

CHAIN_19_TO_38 <- rtracklayer::import.chain(
  CHAIN_IMPORT_FILE
)

cat(
  "\nChain loaded:\n",
  CHAIN_IMPORT_FILE,
  "\n",
  sep = ""
)

# ==============================================================================
# 6. AUTO-DETECT COMPLETE 1000G EUR PLINK REFERENCE
# ==============================================================================

find_reference_prefix <- function() {

  chr1_files <- list.files(
    DATA_ROOT,
    pattern = "^1000G[.]EUR[.]QC[.]1[.]bim$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(chr1_files)) {
    stop(
      "Could not locate 1000G.EUR.QC.1.bim under D:/A/data."
    )
  }

  candidates <- unique(
    sub(
      "1[.]bim$",
      "",
      norm_path(
        chr1_files
      ),
      ignore.case = TRUE
    )
  )

  is_complete <- vapply(
    candidates,
    function(pre) {

      all(
        vapply(
          1:22,
          function(chr) {
            all(
              file.exists(
                paste0(
                  pre,
                  chr,
                  c(
                    ".bed",
                    ".bim",
                    ".fam"
                  )
                )
              )
            )
          },
          logical(1)
        )
      )
    },
    logical(1)
  )

  candidates <- candidates[
    is_complete
  ]

  if (!length(candidates)) {
    stop(
      "1000G EUR chr1 candidates were found, but no complete chr1-22 family exists."
    )
  }

  score <- 1000 *
    grepl(
      "STEP10D1",
      candidates,
      ignore.case = TRUE
    ) +
    500 *
    grepl(
      "02_LDSC_reference",
      candidates,
      ignore.case = TRUE
    ) -
    nchar(candidates) / 10000

  candidates[
    which.max(score)
  ]
}

REF_PREFIX <- find_reference_prefix()

cat(
  "\n1000G EUR PLINK prefix:\n",
  REF_PREFIX,
  "\n",
  sep = ""
)

# ==============================================================================
# 7. AUTO-DETECT PLINK 1.9
# ==============================================================================

find_plink <- function() {

  preferred <- c(
    "D:/A/data/MR_最终结果/STEP3_MR/00_reference/plink/plink.exe",
    "D:/A/data/MR_final/STEP3_MR/00_reference/plink/plink.exe",
    "D:/A/plink/plink.exe",
    "D:/A/data/plink/plink.exe"
  )

  preferred <- preferred[
    file.exists(preferred)
  ]

  if (length(preferred)) {
    return(
      norm_path(
        preferred[1],
        TRUE
      )
    )
  }

  hits <- list.files(
    "D:/A",
    pattern = "^plink[.]exe$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(hits)) {
    stop(
      "PLINK 1.9 plink.exe was not found under D:/A."
    )
  }

  hits <- hits[
    !grepl(
      "plink2",
      hits,
      ignore.case = TRUE
    )
  ]

  if (!length(hits)) {
    stop(
      "Only PLINK2 candidates were found; this script expects PLINK 1.9."
    )
  }

  hits[
    order(
      nchar(
        norm_path(
          hits
        )
      )
    )
  ][1]
}

PLINK_EXE <- find_plink()

PLINK_VERSION_LOG <- file.path(
  QC_DIR,
  "STEP11G2B_plink_version.txt"
)

system2(
  PLINK_EXE,
  args = "--version",
  stdout = PLINK_VERSION_LOG,
  stderr = PLINK_VERSION_LOG
)

cat(
  "\nPLINK:\n",
  PLINK_EXE,
  "\n",
  sep = ""
)

# ==============================================================================
# 8. MEMORY-SAFE ONE-TIME SPLIT OF 31.6M AF GWAS BY CHROMOSOME
# ==============================================================================

find_python <- function() {

  candidates <- c(
    "E:/app/python/envs/ldsc39/python.exe",
    "E:/app/python/python.exe",
    Sys.which("python"),
    Sys.which("python3")
  )

  candidates <- unique(
    candidates[
      nzchar(candidates) &
        file.exists(candidates)
    ]
  )

  if (!length(candidates)) {
    stop(
      "No Python executable was found for the one-time memory-safe GWAS split."
    )
  }

  norm_path(
    candidates[1],
    TRUE
  )
}

PYTHON_EXE <- find_python()

SPLIT_MANIFEST <- file.path(
  AF_SPLIT_DIR,
  "AF_chr_split_manifest.csv"
)

chr_files_expected <- file.path(
  AF_SPLIT_DIR,
  sprintf(
    "AF_GCST90624412_chr%02d.tsv.gz",
    1:22
  )
)

split_ready <- (
  file.exists(
    SPLIT_MANIFEST
  ) &&
  all(
    file.exists(
      chr_files_expected
    )
  ) &&
  all(
    file.info(
      chr_files_expected
    )$size > 100
  )
)

if (!split_ready) {

  PY_SPLIT_SCRIPT <- file.path(
    TMP_DIR,
    "split_AF_GWAS_by_chr.py"
  )

  py_lines <- c(
    "import gzip, csv, os, sys",
    "inp = sys.argv[1]",
    "outdir = sys.argv[2]",
    "manifest = sys.argv[3]",
    "os.makedirs(outdir, exist_ok=True)",
    "handles = {}",
    "counts = {str(i): 0 for i in range(1,23)}",
    "try:",
    "    with gzip.open(inp, 'rt', encoding='utf-8', newline='') as fin:",
    "        header = fin.readline()",
    "        if not header:",
    "            raise RuntimeError('Empty standardized AF GWAS')",
    "        cols = header.rstrip('\\n\\r').split('\\t')",
    "        if 'CHR_SOURCE' not in cols:",
    "            raise RuntimeError('CHR_SOURCE not found in standardized AF header')",
    "        chr_idx = cols.index('CHR_SOURCE')",
    "        for c in range(1,23):",
    "            p = os.path.join(outdir, f'AF_GCST90624412_chr{c:02d}.tsv.gz')",
    "            h = gzip.open(p, 'wt', encoding='utf-8', newline='')",
    "            h.write(header)",
    "            handles[str(c)] = h",
    "        for line in fin:",
    "            if not line.strip():",
    "                continue",
    "            parts = line.rstrip('\\n\\r').split('\\t')",
    "            if chr_idx >= len(parts):",
    "                continue",
    "            ch = parts[chr_idx].replace('chr','').replace('CHR','')",
    "            if ch in handles:",
    "                handles[ch].write(line)",
    "                counts[ch] += 1",
    "finally:",
    "    for h in handles.values():",
    "        try: h.close()",
    "        except Exception: pass",
    "with open(manifest, 'w', newline='', encoding='utf-8') as f:",
    "    w = csv.writer(f)",
    "    w.writerow(['CHR','n_rows','file'])",
    "    for c in range(1,23):",
    "        p = os.path.join(outdir, f'AF_GCST90624412_chr{c:02d}.tsv.gz')",
    "        w.writerow([c, counts[str(c)], p])"
  )

  writeLines(
    py_lines,
    PY_SPLIT_SCRIPT
  )

  SPLIT_LOG <- file.path(
    LOG_DIR,
    "STEP11G2B_AF_split_console.log"
  )

  cat(
    "\nSplitting standardized AF GWAS by chromosome (one-time streaming pass)...\n"
  )

  split_status <- suppressWarnings(
    system2(
      PYTHON_EXE,
      args = c(
        shQuote(
          PY_SPLIT_SCRIPT
        ),
        shQuote(
          norm_path(
            AF_STANDARD_FILE,
            TRUE
          )
        ),
        shQuote(
          norm_path(
            AF_SPLIT_DIR,
            FALSE
          )
        ),
        shQuote(
          norm_path(
            SPLIT_MANIFEST,
            FALSE
          )
        )
      ),
      stdout = SPLIT_LOG,
      stderr = SPLIT_LOG
    )
  )

  if (
    !identical(
      as.integer(split_status),
      0L
    ) ||
    !file.exists(
      SPLIT_MANIFEST
    ) ||
    !all(
      file.exists(
        chr_files_expected
      )
    )
  ) {
    stop(
      paste0(
        "AF chromosome split failed.\nInspect:\n",
        SPLIT_LOG
      )
    )
  }
}

SPLIT_QC <- fread(
  SPLIT_MANIFEST
)

if (
  nrow(SPLIT_QC) != 22L ||
  !all(
    1:22 %in%
      SPLIT_QC$CHR
  ) ||
  any(
    SPLIT_QC$n_rows <= 0
  )
) {
  stop(
    "AF chromosome split manifest failed QC."
  )
}

fwrite(
  SPLIT_QC,
  file.path(
    QC_DIR,
    "STEP11G2B_AF_chr_split_QC.csv"
  )
)

# ==============================================================================
# 9. susieR CAPABILITY AUDIT
# ==============================================================================

SUSIE_FORMALS <- names(
  formals(
    susieR::susie_rss
  )
)

HAS_R_FINITE <- "R_finite" %in% SUSIE_FORMALS
HAS_R_MISMATCH <- "R_mismatch" %in% SUSIE_FORMALS

SUSIE_CAP <- data.table(
  field = c(
    "susieR_version",
    "low_rank_X_supported",
    "R_finite_supported",
    "R_mismatch_supported",
    "primary_L",
    "primary_prior",
    "primary_R_mismatch",
    "sensitivity_enabled"
  ),
  value = c(
    as.character(
      packageVersion(
        "susieR"
      )
    ),
    "TRUE",
    as.character(
      HAS_R_FINITE
    ),
    as.character(
      HAS_R_MISMATCH
    ),
    as.character(
      L_PRIMARY
    ),
    "uniform",
    "none",
    as.character(
      HAS_R_FINITE &&
        HAS_R_MISMATCH
    )
  )
)

fwrite(
  SUSIE_CAP,
  file.path(
    QC_DIR,
    "STEP11G2B_susieR_capabilities.csv"
  )
)

# ==============================================================================
# 10. INITIALIZE / RESUME BLOCK STATUS
# ==============================================================================

if (file.exists(STATUS_FILE)) {

  STATUS <- fread(
    STATUS_FILE
  )

} else {

  STATUS <- data.table()
}

# ==============================================================================
# 11. CHROMOSOME-BY-CHROMOSOME FINE-MAPPING
# ==============================================================================

for (chr in 1:22) {

  chr_blocks <- BLOCKS[
    CHR == chr
  ]

  if (!nrow(chr_blocks)) {
    next
  }

  cat(
    "\n============================================================\n",
    "CHR ",
    chr,
    " — ",
    nrow(chr_blocks),
    " selected AF blocks\n",
    "============================================================\n",
    sep = ""
  )

  AF_CHR_FILE <- file.path(
    AF_SPLIT_DIR,
    sprintf(
      "AF_GCST90624412_chr%02d.tsv.gz",
      chr
    )
  )

  if (!file.exists(AF_CHR_FILE)) {
    stop(
      "Missing chromosome-split AF file: ",
      AF_CHR_FILE
    )
  }

  cat(
    "Loading AF chr",
    chr,
    "...\n",
    sep = ""
  )

  AFCHR <- fread(
    AF_CHR_FILE,
    select = c(
      "CHR_SOURCE",
      "BP_SOURCE",
      "A1",
      "A2",
      "BETA",
      "SE",
      "Z",
      "P",
      "EAF",
      "N"
    ),
    showProgress = TRUE
  )

  AFCHR[
    ,
    `:=`(
      A1 = toupper(
        A1
      ),
      A2 = toupper(
        A2
      ),
      allele_key = allele_key(
        A1,
        A2
      )
    )
  ]

  setkey(
    AFCHR,
    BP_SOURCE,
    allele_key
  )

  BIM_FILE <- paste0(
    REF_PREFIX,
    chr,
    ".bim"
  )

  BIM <- fread(
    BIM_FILE,
    header = FALSE,
    col.names = c(
      "CHR_REF",
      "SNP",
      "CM",
      "BP_HG19",
      "REF_A1",
      "REF_A2"
    ),
    showProgress = FALSE
  )

  BIM[
    ,
    `:=`(
      REF_A1 = toupper(
        REF_A1
      ),
      REF_A2 = toupper(
        REF_A2
      ),
      allele_key = allele_key(
        REF_A1,
        REF_A2
      )
    )
  ]

  for (ii in seq_len(nrow(chr_blocks))) {

    b <- chr_blocks[ii]

    LABEL <- b$selected_block_label
    ORDER <- b$selected_order

    BLOCK_DIR <- file.path(
      BLOCK_OUT_DIR,
      LABEL
    )

    dir.create(
      BLOCK_DIR,
      recursive = TRUE,
      showWarnings = FALSE
    )

    PIP_FILE <- file.path(
      BLOCK_DIR,
      paste0(
        LABEL,
        "_PIP.tsv.gz"
      )
    )

    SUMMARY_FILE <- file.path(
      BLOCK_DIR,
      paste0(
        LABEL,
        "_summary.csv"
      )
    )

    FIT_PRIMARY_FILE <- file.path(
      BLOCK_DIR,
      paste0(
        LABEL,
        "_susie_primary.rds"
      )
    )

    # ----------------------------------------------------------
    # RESUME: skip a verified PASS block.
    # ----------------------------------------------------------

    if (
      file.exists(PIP_FILE) &&
      file.exists(SUMMARY_FILE)
    ) {

      prev <- tryCatch(
        fread(
          SUMMARY_FILE
        ),
        error = function(e) NULL
      )

      if (
        !is.null(prev) &&
        nrow(prev) == 1L &&
        identical(
          as.character(
            prev$status
          ),
          "PASS"
        )
      ) {

        cat(
          "[",
          LABEL,
          "] checkpoint PASS — skip\n",
          sep = ""
        )

        next
      }
    }

    cat(
      "\n[",
      LABEL,
      "] chr",
      chr,
      ":",
      b$start_hg19,
      "-",
      b$end_hg19,
      " | GWS=",
      b$n_gws_variants,
      "\n",
      sep = ""
    )

    t0 <- Sys.time()

    result_row <- tryCatch(
      {

        # ======================================================
        # 11A. REFERENCE VARIANTS IN THIS hg19 LDETECT BLOCK
        # ======================================================

        REF_B <- BIM[
          BP_HG19 >=
            b$start_hg19 &
            BP_HG19 <=
            b$end_hg19
        ]

        # Keep MHC excluded even if a selected block straddles
        # the upper MHC boundary.
        if (chr == MHC_CHR) {

          REF_B <- REF_B[
            !(
              BP_HG19 >= MHC_START &
              BP_HG19 <= MHC_END
            )
          ]
        }

        REF_B <- REF_B[
          valid_base_allele(
            REF_A1
          ) &
            valid_base_allele(
              REF_A2
            ) &
            REF_A1 != REF_A2
        ]

        REF_B <- unique(
          REF_B,
          by = "SNP"
        )

        n_ref_block <- nrow(
          REF_B
        )

        if (
          n_ref_block <
            MIN_MATCHED_SNPS
        ) {
          stop(
            "Too few 1000G EUR variants in block: ",
            n_ref_block
          )
        }

        # ======================================================
        # 11B. LIFT REFERENCE VARIANTS hg19 -> GRCh38
        #       (much smaller than lifting the whole 31.6M GWAS)
        # ======================================================

        gr19 <- GRanges(
          seqnames = paste0(
            "chr",
            chr
          ),
          ranges = IRanges(
            start = REF_B$BP_HG19,
            width = 1L
          )
        )

        lifted <- rtracklayer::liftOver(
          gr19,
          CHAIN_19_TO_38
        )

        nmap <- lengths(
          lifted
        )

        unique_idx <- which(
          nmap == 1L
        )

        if (
          length(unique_idx) <
            MIN_MATCHED_SNPS
        ) {
          stop(
            "Too few uniquely lifted 1000G variants: ",
            length(unique_idx)
          )
        }

        ugr <- unlist(
          lifted[
            unique_idx
          ],
          use.names = FALSE
        )

        source_chr <- suppressWarnings(
          as.integer(
            gsub(
              "^chr",
              "",
              as.character(
                seqnames(
                  ugr
                )
              )
            )
          )
        )

        keep_same_chr <- which(
          source_chr == chr
        )

        REF_MAP <- REF_B[
          unique_idx[
            keep_same_chr
          ]
        ]

        REF_MAP[
          ,
          BP_SOURCE :=
            start(
              ugr[
                keep_same_chr
              ]
            )
        ]

        n_ref_unique_lift <- nrow(
          REF_MAP
        )

        # ======================================================
        # 11C. JOIN TO ORIGINAL GRCh38 AF GWAS
        #       exact position + unordered allele pair
        # ======================================================

        MATCH <- merge(
          REF_MAP,
          AFCHR,
          by = c(
            "BP_SOURCE",
            "allele_key"
          ),
          all = FALSE,
          sort = FALSE
        )

        # Defensive: one GWAS row per reference SNP.
        setorder(
          MATCH,
          SNP,
          P
        )

        MATCH <- MATCH[
          ,
          .SD[1],
          by = SNP
        ]

        # Align GWAS effect allele to the BIM A1 allele.
        MATCH[
          ,
          align_to_ref_A1 :=
            fifelse(
              A1 == REF_A1 &
                A2 == REF_A2,
              1,
              fifelse(
                A1 == REF_A2 &
                  A2 == REF_A1,
                -1,
                NA_real_
              )
            )
        ]

        MATCH <- MATCH[
          is.finite(
            align_to_ref_A1
          ) &
            is.finite(
              Z
            ) &
            is.finite(
              BETA
            ) &
            is.finite(
              SE
            ) &
            SE > 0
        ]

        n_matched_preplink <- nrow(
          MATCH
        )

        if (
          n_matched_preplink <
            MIN_MATCHED_SNPS
        ) {
          stop(
            "Too few allele-matched GWAS/reference variants: ",
            n_matched_preplink
          )
        }

        fwrite(
          MATCH,
          file.path(
            BLOCK_DIR,
            paste0(
              LABEL,
              "_harmonized_prePLINK.tsv.gz"
            )
          ),
          sep = "\t",
          compress = "gzip"
        )

        # ======================================================
        # 11D. PLINK EXPORT OF REFERENCE GENOTYPES
        #      A-transpose: variants x ~500 EUR reference samples
        # ======================================================

        SNP_LIST <- file.path(
          BLOCK_DIR,
          paste0(
            LABEL,
            "_matched_snps.txt"
          )
        )

        writeLines(
          MATCH$SNP,
          SNP_LIST
        )

        PLINK_OUT <- file.path(
          TMP_DIR,
          LABEL
        )

        PLINK_LOG <- file.path(
          LOG_DIR,
          paste0(
            LABEL,
            "_plink.log"
          )
        )

        plink_args <- c(
          "--bfile",
          shQuote(
            paste0(
              REF_PREFIX,
              chr
            )
          ),
          "--extract",
          shQuote(
            SNP_LIST
          ),
          "--keep-allele-order",
          "--recode",
          "A-transpose",
          "--write-snplist",
          "--allow-no-sex",
          "--threads",
          as.character(
            PLINK_THREADS
          ),
          "--memory",
          as.character(
            PLINK_MEMORY_MB
          ),
          "--out",
          shQuote(
            PLINK_OUT
          )
        )

        plink_status <- suppressWarnings(
          system2(
            PLINK_EXE,
            args = plink_args,
            stdout = PLINK_LOG,
            stderr = PLINK_LOG
          )
        )

        TRAW_FILE <- paste0(
          PLINK_OUT,
          ".traw"
        )

        if (
          !identical(
            as.integer(
              plink_status
            ),
            0L
          ) ||
          !file.exists(
            TRAW_FILE
          )
        ) {
          stop(
            "PLINK genotype export failed. Inspect: ",
            PLINK_LOG
          )
        }

        TRAW <- fread(
          TRAW_FILE,
          check.names = FALSE,
          showProgress = FALSE
        )

        if (
          nrow(TRAW) <
            MIN_MATCHED_SNPS ||
          ncol(TRAW) <= 6L
        ) {
          stop(
            "Unexpected PLINK .traw dimensions: ",
            nrow(TRAW),
            " x ",
            ncol(TRAW)
          )
        }

        # PLINK .traw first six fields:
        # CHR SNP (C)M POS COUNTED ALT
        TRAW_META <- data.table(
          SNP = as.character(
            TRAW[[2]]
          ),
          COUNTED = toupper(
            as.character(
              TRAW[[5]]
            )
          ),
          ALT = toupper(
            as.character(
              TRAW[[6]]
            )
          ),
          traw_row = seq_len(
            nrow(TRAW)
          )
        )

        H <- merge(
          TRAW_META,
          MATCH,
          by = "SNP",
          all = FALSE,
          sort = FALSE
        )

        setorder(
          H,
          traw_row
        )

        # Align z to the ACTUAL allele counted in PLINK .traw.
        H[
          ,
          align_to_COUNTED :=
            fifelse(
              A1 == COUNTED &
                A2 == ALT,
              1,
              fifelse(
                A1 == ALT &
                  A2 == COUNTED,
                -1,
                NA_real_
              )
            )
        ]

        H <- H[
          is.finite(
            align_to_COUNTED
          )
        ]

        if (
          nrow(H) <
            MIN_MATCHED_SNPS
        ) {
          stop(
            "Too few variants remain after aligning to PLINK COUNTED allele."
          )
        }

        # Reorder TRAW exactly to harmonized H order.
        TRAW_SUB <- TRAW[
          H$traw_row
        ]

        G <- as.matrix(
          TRAW_SUB[
            ,
            7:ncol(TRAW_SUB),
            with = FALSE
          ]
        )

        storage.mode(
          G
        ) <- "double"

        # ======================================================
        # 11E. MEAN IMPUTE THE RARE MISSING REFERENCE GENOTYPES
        # ======================================================

        row_means <- rowMeans(
          G,
          na.rm = TRUE
        )

        bad_mean <- !is.finite(
          row_means
        )

        if (any(bad_mean)) {

          keep <- !bad_mean

          G <- G[
            keep,
            ,
            drop = FALSE
          ]

          H <- H[
            keep
          ]

          row_means <- row_means[
            keep
          ]
        }

        na_idx <- which(
          is.na(G),
          arr.ind = TRUE
        )

        if (nrow(na_idx)) {
          G[
            na_idx
          ] <- row_means[
            na_idx[
              ,
              1
            ]
          ]
        }

        # Remove zero-variance reference variants.
        row_var <- apply(
          G,
          1,
          var
        )

        polymorphic <- is.finite(
          row_var
        ) &
          row_var > 0

        G <- G[
          polymorphic,
          ,
          drop = FALSE
        ]

        H <- H[
          polymorphic
        ]

        n_final <- nrow(
          H
        )

        if (
          n_final <
            MIN_MATCHED_SNPS
        ) {
          stop(
            "Too few polymorphic matched variants remain: ",
            n_final
          )
        }

        # samples x variants factor matrix
        X_REF <- t(
          G
        )

        rm(
          G,
          TRAW,
          TRAW_SUB
        )

        gc(
          verbose = FALSE
        )

        Z_REF <- H$Z *
          H$align_to_COUNTED

        N_EFF <- median(
          H$N[
            is.finite(
              H$N
            ) &
              H$N > 0
          ],
          na.rm = TRUE
        )

        if (!is.finite(N_EFF)) {
          N_EFF <- 1840341
        }

        B_REF <- nrow(
          X_REF
        )

        # ======================================================
        # 11F. PRIMARY SuSiE-RSS
        #      L=1, uniform prior, no functional prior
        # ======================================================

        primary_args <- list(
          z = Z_REF,
          X = X_REF,
          n = N_EFF,
          L = L_PRIMARY,
          prior_weights = NULL,
          estimate_residual_variance = FALSE,
          coverage = COVERAGE,
          min_abs_corr = MIN_ABS_CORR,
          max_iter = 100,
          tol = 1e-4,
          verbose = FALSE
        )

        if (HAS_R_FINITE) {
          primary_args$R_finite <- NULL
        }

        if (HAS_R_MISMATCH) {
          primary_args$R_mismatch <- "none"
        }

        FIT_PRIMARY <- do.call(
          susieR::susie_rss,
          primary_args
        )

        if (
          is.null(
            FIT_PRIMARY$pip
          ) ||
          length(
            FIT_PRIMARY$pip
          ) != n_final
        ) {
          stop(
            "Primary SuSiE did not return one PIP per variant."
          )
        }

        # ======================================================
        # 11G. CURRENT-METHOD SENSITIVITY:
        #      finite-reference + empirical-Bayes LD mismatch
        # ======================================================

        FIT_SENS <- NULL

        if (
          HAS_R_FINITE &&
          HAS_R_MISMATCH
        ) {

          sens_args <- primary_args
          sens_args$R_finite <- TRUE
          sens_args$R_mismatch <- "eb"

          FIT_SENS <- tryCatch(
            do.call(
              susieR::susie_rss,
              sens_args
            ),
            error = function(e) {
              message(
                "[",
                LABEL,
                "] sensitivity fit skipped after error: ",
                conditionMessage(e)
              )
              NULL
            }
          )
        }

        # ======================================================
        # 11H. CREDIBLE SET + PIP EXPORT
        # ======================================================

        cs_primary_idx <- integer(0)

        if (
          !is.null(
            FIT_PRIMARY$sets
          ) &&
          !is.null(
            FIT_PRIMARY$sets$cs
          ) &&
          length(
            FIT_PRIMARY$sets$cs
          )
        ) {
          cs_primary_idx <- sort(
            unique(
              as.integer(
                unlist(
                  FIT_PRIMARY$sets$cs,
                  use.names = FALSE
                )
              )
            )
          )
        }

        PIP_PRIMARY <- as.numeric(
          FIT_PRIMARY$pip
        )

        PIP_SENS <- rep(
          NA_real_,
          n_final
        )

        if (
          !is.null(FIT_SENS) &&
          !is.null(FIT_SENS$pip) &&
          length(FIT_SENS$pip) ==
            n_final
        ) {
          PIP_SENS <- as.numeric(
            FIT_SENS$pip
          )
        }

        OUT <- H[
          ,
          .(
            selected_block_label = LABEL,
            block_id = b$block_id,
            CHR_HG19 = CHR_REF,
            BP_HG19,
            BP_GRCh38 = BP_SOURCE,
            SNP,
            REF_A1,
            REF_A2,
            COUNTED,
            ALT,
            GWAS_A1 = A1,
            GWAS_A2 = A2,
            BETA,
            SE,
            Z_GWAS = Z,
            Z_ref_counted = Z_REF,
            P,
            EAF,
            N,
            PIP_primary = PIP_PRIMARY,
            PIP_sensitivity = PIP_SENS,
            in_primary_CS95 =
              seq_len(.N) %in%
              cs_primary_idx
          )
        ]

        setorder(
          OUT,
          -PIP_primary
        )

        fwrite(
          OUT,
          PIP_FILE,
          sep = "\t",
          compress = "gzip"
        )

        saveRDS(
          FIT_PRIMARY,
          FIT_PRIMARY_FILE,
          compress = TRUE
        )

        FIT_SENS_FILE <- NA_character_

        if (!is.null(FIT_SENS)) {

          FIT_SENS_FILE <- file.path(
            BLOCK_DIR,
            paste0(
              LABEL,
              "_susie_sensitivity.rds"
            )
          )

          saveRDS(
            FIT_SENS,
            FIT_SENS_FILE,
            compress = TRUE
          )
        }

        # ======================================================
        # 11I. BLOCK QC
        # ======================================================

        primary_converged <- if (
          !is.null(
            FIT_PRIMARY$converged
          )
        ) {
          isTRUE(
            FIT_PRIMARY$converged
          )
        } else {
          TRUE
        }

        sens_converged <- if (
          is.null(FIT_SENS)
        ) {
          NA
        } else if (
          !is.null(
            FIT_SENS$converged
          )
        ) {
          isTRUE(
            FIT_SENS$converged
          )
        } else {
          TRUE
        }

        pip_cor <- if (
          all(
            is.finite(
              PIP_SENS
            )
          ) &&
          length(
            unique(
              PIP_SENS
            )
          ) > 1L
        ) {
          suppressWarnings(
            cor(
              PIP_PRIMARY,
              PIP_SENS,
              method = "spearman"
            )
          )
        } else {
          NA_real_
        }

        diag_sens <- if (
          !is.null(FIT_SENS)
        ) {
          FIT_SENS$R_finite_diagnostics
        } else {
          NULL
        }

        r_over_B <- safe_scalar(
          if (!is.null(diag_sens)) {
            diag_sens$r_over_B
          } else {
            NULL
          }
        )

        lambda_bias <- safe_scalar(
          if (!is.null(diag_sens)) {
            diag_sens$lambda_bias
          } else {
            NULL
          }
        )

        B_corrected <- safe_scalar(
          if (!is.null(diag_sens)) {
            diag_sens$B_corrected
          } else {
            NULL
          }
        )

        summary_dt <- data.table(
          selected_block_label = LABEL,
          block_id = b$block_id,
          selected_order = ORDER,
          CHR = chr,
          start_hg19 = b$start_hg19,
          end_hg19 = b$end_hg19,
          status = "PASS",
          error_message = NA_character_,
          n_GWS = b$n_gws_variants,
          n_ref_block = n_ref_block,
          n_ref_unique_liftover =
            n_ref_unique_lift,
          n_matched_prePLINK =
            n_matched_preplink,
          n_final_variants =
            n_final,
          n_reference_samples =
            B_REF,
          n_eff =
            N_EFF,
          reference_match_fraction =
            n_matched_preplink /
            max(
              1,
              n_ref_unique_lift
            ),
          primary_converged =
            primary_converged,
          sum_PIP_primary =
            sum(
              PIP_PRIMARY,
              na.rm = TRUE
            ),
          max_PIP_primary =
            max(
              PIP_PRIMARY,
              na.rm = TRUE
            ),
          n_PIP_primary_ge_0_1 =
            sum(
              PIP_PRIMARY >=
                0.1,
              na.rm = TRUE
            ),
          n_PIP_primary_ge_0_01 =
            sum(
              PIP_PRIMARY >=
                0.01,
              na.rm = TRUE
            ),
          primary_CS95_size =
            length(
              cs_primary_idx
            ),
          sensitivity_run =
            !is.null(
              FIT_SENS
            ),
          sensitivity_converged =
            sens_converged,
          primary_vs_sensitivity_PIP_spearman =
            pip_cor,
          sensitivity_r_over_B =
            r_over_B,
          sensitivity_lambda_bias =
            lambda_bias,
          sensitivity_B_corrected =
            B_corrected,
          elapsed_seconds =
            as.numeric(
              difftime(
                Sys.time(),
                t0,
                units = "secs"
              )
            ),
          pip_file =
            norm_path(
              PIP_FILE,
              TRUE
            ),
          primary_fit_file =
            norm_path(
              FIT_PRIMARY_FILE,
              TRUE
            )
        )

        fwrite(
          summary_dt,
          SUMMARY_FILE
        )

        if (!KEEP_TRAW) {

          unlink(
            c(
              paste0(
                PLINK_OUT,
                ".traw"
              ),
              paste0(
                PLINK_OUT,
                ".snplist"
              ),
              paste0(
                PLINK_OUT,
                ".nosex"
              ),
              paste0(
                PLINK_OUT,
                ".log"
              )
            )[
              file.exists(
                c(
                  paste0(
                    PLINK_OUT,
                    ".traw"
                  ),
                  paste0(
                    PLINK_OUT,
                    ".snplist"
                  ),
                  paste0(
                    PLINK_OUT,
                    ".nosex"
                  ),
                  paste0(
                    PLINK_OUT,
                    ".log"
                  )
                )
              )
            ]
          )
        }

        rm(
          REF_B,
          REF_MAP,
          MATCH,
          TRAW_META,
          H,
          X_REF,
          FIT_PRIMARY,
          FIT_SENS,
          OUT
        )

        gc(
          verbose = FALSE
        )

        summary_dt
      },

      error = function(e) {

        msg <- conditionMessage(
          e
        )

        fail_dt <- data.table(
          selected_block_label = LABEL,
          block_id = b$block_id,
          selected_order = ORDER,
          CHR = chr,
          start_hg19 = b$start_hg19,
          end_hg19 = b$end_hg19,
          status = "FAIL",
          error_message = msg,
          n_GWS = b$n_gws_variants,
          elapsed_seconds =
            as.numeric(
              difftime(
                Sys.time(),
                t0,
                units = "secs"
              )
            )
        )

        fwrite(
          fail_dt,
          SUMMARY_FILE
        )

        message(
          "[",
          LABEL,
          "] FAIL: ",
          msg
        )

        fail_dt
      }
    )

    STATUS <- status_update(
      result_row
    )

    cat(
      "[",
      LABEL,
      "] ",
      result_row$status,
      if (
        identical(
          result_row$status,
          "PASS"
        )
      ) {
        paste0(
          " | n=",
          result_row$n_final_variants,
          " | max PIP=",
          signif(
            result_row$max_PIP_primary,
            4
          )
        )
      } else {
        paste0(
          " | ",
          result_row$error_message
        )
      },
      "\n",
      sep = ""
    )
  }

  rm(
    AFCHR,
    BIM
  )

  gc(
    verbose = FALSE
  )
}

# ==============================================================================
# 12. FINAL STATUS AUDIT
# ==============================================================================

STATUS <- fread(
  STATUS_FILE
)

STATUS <- merge(
  BLOCKS[
    ,
    .(
      selected_block_label,
      expected_order =
        selected_order
    )
  ],
  STATUS,
  by = "selected_block_label",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  STATUS,
  expected_order
)

fwrite(
  STATUS,
  file.path(
    QC_DIR,
    "STEP11G2B_block_status_FINAL.csv"
  )
)

PASS_BLOCKS <- STATUS[
  status == "PASS"
]

FAIL_BLOCKS <- STATUS[
  is.na(status) |
    status != "PASS"
]

fwrite(
  FAIL_BLOCKS,
  file.path(
    QC_DIR,
    "STEP11G2B_failed_blocks.csv"
  )
)

cat(
  "\n============================================================\n",
  "BLOCK STATUS\n",
  "============================================================\n",
  "Expected blocks: ",
  nrow(BLOCKS),
  "\nPASS: ",
  nrow(PASS_BLOCKS),
  "\nFAIL/incomplete: ",
  nrow(FAIL_BLOCKS),
  "\n",
  sep = ""
)

# ==============================================================================
# 13. COMBINE ALL PRIMARY PIP RESULTS
# ==============================================================================

PIP_LIST <- list()

for (lab in PASS_BLOCKS$selected_block_label) {

  f <- file.path(
    BLOCK_OUT_DIR,
    lab,
    paste0(
      lab,
      "_PIP.tsv.gz"
    )
  )

  if (file.exists(f)) {
    PIP_LIST[[
      length(PIP_LIST) + 1L
    ]] <- fread(
      f,
      showProgress = FALSE
    )
  }
}

ALL_PIP <- if (length(PIP_LIST)) {
  rbindlist(
    PIP_LIST,
    fill = TRUE
  )
} else {
  data.table()
}

if (nrow(ALL_PIP)) {

  setorder(
    ALL_PIP,
    CHR_HG19,
    BP_HG19,
    -PIP_primary
  )

  ALL_PIP_FILE <- file.path(
    FINAL_DIR,
    "AF_SCAVENGE_PIP_PRIMARY_ALL.tsv.gz"
  )

  fwrite(
    ALL_PIP,
    ALL_PIP_FILE,
    sep = "\t",
    compress = "gzip"
  )

  PIP001 <- ALL_PIP[
    PIP_primary >=
      0.001
  ]

  fwrite(
    PIP001,
    file.path(
      FINAL_DIR,
      "AF_SCAVENGE_PIP_PRIMARY_GE_0.001.tsv.gz"
    ),
    sep = "\t",
    compress = "gzip"
  )

  BED <- ALL_PIP[
    PIP_primary > 0,
    .(
      CHR = paste0(
        "chr",
        CHR_HG19
      ),
      START0 = pmax(
        0L,
        BP_HG19 - 1L
      ),
      END = BP_HG19,
      SNP,
      PIP = PIP_primary,
      selected_block_label
    )
  ]

  fwrite(
    BED,
    file.path(
      FINAL_DIR,
      "AF_SCAVENGE_PIP_PRIMARY_hg19.bed"
    ),
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )

  # Sensitivity PIP is saved separately and never silently substituted.
  if (
    "PIP_sensitivity" %in%
      names(ALL_PIP) &&
    any(
      is.finite(
        ALL_PIP$PIP_sensitivity
      )
    )
  ) {

    fwrite(
      ALL_PIP[
        is.finite(
          PIP_sensitivity
        )
      ],
      file.path(
        FINAL_DIR,
        "AF_SCAVENGE_PIP_SENSITIVITY_FINITE_REF_EB.tsv.gz"
      ),
      sep = "\t",
      compress = "gzip"
    )
  }

} else {

  ALL_PIP_FILE <- NA_character_
}

# ==============================================================================
# 14. FINAL QC TABLES
# ==============================================================================

BLOCK_QC <- STATUS[
  status == "PASS",
  .(
    n_pass_blocks = .N,
    median_final_variants =
      median(
        n_final_variants,
        na.rm = TRUE
      ),
    max_final_variants =
      max(
        n_final_variants,
        na.rm = TRUE
      ),
    median_reference_match_fraction =
      median(
        reference_match_fraction,
        na.rm = TRUE
      ),
    median_max_PIP =
      median(
        max_PIP_primary,
        na.rm = TRUE
      ),
    n_primary_nonconverged =
      sum(
        !primary_converged,
        na.rm = TRUE
      ),
    median_primary_vs_sensitivity_spearman =
      median(
        primary_vs_sensitivity_PIP_spearman,
        na.rm = TRUE
      )
  )
]

fwrite(
  BLOCK_QC,
  file.path(
    QC_DIR,
    "STEP11G2B_finemap_summary_QC.csv"
  )
)

N_DUP_SNP <- if (nrow(ALL_PIP)) {
  sum(
    duplicated(
      ALL_PIP$SNP
    )
  )
} else {
  NA_integer_
}

READINESS <- data.table(
  check = c(
    "STEP11G2A_all_pass",
    "AF_chr1_22_split_complete",
    "hg19ToHg38_chain_loaded",
    "1000G_EUR_chr1_22_complete",
    "PLINK_detected",
    "susieR_low_rank_X_supported",
    "all_selected_blocks_attempted",
    "all_selected_blocks_PASS",
    "all_primary_SuSiE_converged",
    "combined_primary_PIP_created",
    "no_duplicate_SNP_across_blocks",
    "STEP11G3_SCAVENGE_ready"
  ),
  pass = c(
    all(
      READY_A$pass %in%
        c(
          TRUE,
          "TRUE",
          1
        )
    ),
    nrow(SPLIT_QC) == 22L &&
      all(
        SPLIT_QC$n_rows > 0
      ),
    length(
      CHAIN_19_TO_38
    ) > 0,
    all(
      vapply(
        1:22,
        function(chr) {
          all(
            file.exists(
              paste0(
                REF_PREFIX,
                chr,
                c(
                  ".bed",
                  ".bim",
                  ".fam"
                )
              )
            )
          )
        },
        logical(1)
      )
    ),
    file.exists(
      PLINK_EXE
    ),
    "X" %in%
      names(
        formals(
          susieR::susie_rss
        )
      ),
    nrow(STATUS) ==
      nrow(BLOCKS) &&
      all(
        !is.na(
          STATUS$status
        )
      ),
    nrow(FAIL_BLOCKS) == 0L,
    nrow(PASS_BLOCKS) ==
      nrow(BLOCKS) &&
      all(
        PASS_BLOCKS$primary_converged %in%
          c(
            TRUE,
            "TRUE",
            1
          )
      ),
    !is.na(
      ALL_PIP_FILE
    ) &&
      file.exists(
        ALL_PIP_FILE
      ) &&
      nrow(
        ALL_PIP
      ) > 0L,
    isTRUE(
      N_DUP_SNP == 0L
    ),
    nrow(FAIL_BLOCKS) == 0L &&
      !is.na(
        ALL_PIP_FILE
      ) &&
      file.exists(
        ALL_PIP_FILE
      )
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G2B_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "primary_AF_GWAS",
    "GWAS_source_build",
    "fine_mapping_build",
    "fine_mapping_partition",
    "n_selected_blocks",
    "LD_reference",
    "reference_genotype_representation",
    "primary_model",
    "primary_L",
    "primary_prior",
    "primary_residual_variance",
    "primary_R_mismatch",
    "sensitivity_model",
    "SCAVENGE_primary_input"
  ),
  value = c(
    "GCST90624412 EUR",
    "GRCh38",
    "GRCh37/hg19",
    "Berisa-Pickrell EUR LDetect blocks",
    as.character(
      nrow(BLOCKS)
    ),
    "1000 Genomes Phase 3 EUR",
    "PLINK A-transpose reference genotype matrix X; susieR low-rank RSS path",
    "susieR::susie_rss",
    "L=1",
    "uniform prior_weights=NULL",
    "estimate_residual_variance=FALSE",
    "none",
    if (
      HAS_R_FINITE &&
      HAS_R_MISMATCH
    ) {
      "R_finite=TRUE + R_mismatch=eb (sensitivity only)"
    } else {
      "not available in installed susieR"
    },
    "AF_SCAVENGE_PIP_PRIMARY_ALL.tsv.gz"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G2B_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G2B_sessionInfo.txt"
  )
)

# ==============================================================================
# 15. FINAL SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G2B FINISHED\n",
  "============================================================\n\n",
  sep = ""
)

print(
  BLOCK_QC
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

if (nrow(FAIL_BLOCKS)) {

  cat(
    "\nIMPORTANT: some blocks failed. Do NOT start SCAVENGE yet.\n",
    "Upload:\n",
    file.path(
      QC_DIR,
      "STEP11G2B_failed_blocks.csv"
    ),
    "\n",
    file.path(
      QC_DIR,
      "STEP11G2B_block_status_FINAL.csv"
    ),
    "\n",
    sep = ""
  )

} else {

  cat(
    "\nALL BLOCKS PASS.\n",
    "Primary SCAVENGE PIP:\n",
    ALL_PIP_FILE,
    "\n\nNext stage: STEP11G3 PIP -> GSE238242 scATAC -> gchromVAR / SCAVENGE.\n",
    sep = ""
  )
}
