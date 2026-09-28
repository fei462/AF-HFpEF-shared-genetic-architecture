# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G3A_V2_FIX_WEIGHT_ASSAY_AF_PIP_TO_GCHROMVAR_SEEDS(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G3A V2 — AF fine-mapped PIP -> GSE238242 single-cell gchromVAR
#              -> initial trait Z score + SCAVENGE seed cells
#
# Project: AF–HFpEF–BMI–OSA shared genetics
# Dataset: GSE238242 left atrial appendage snATAC
#
# FROZEN DESIGN
#   - 7 donors
#   - 11,986 author-annotated cells
#   - 12 author cell types
#   - 212,084 GRCh38 ATAC peaks
#   - NO reclustering
#   - NO reannotation
#   - NO biological filtering by cell type / rhythm
#   - all cells enter the primary SCAVENGE analysis
#
# INPUT
#   STEP11G2B:
#     AF_SCAVENGE_PIP_PRIMARY_ALL.tsv.gz
#     primary SuSiE PIP from 515 EUR LDetect blocks
#
# PRIMARY SCAVENGE VARIANT INPUT
#   PIP_primary >= 0.001
#
# RATIONALE
#   This threshold preserves ~98.5% of the primary PIP mass in the frozen
#   STEP11G2B result and retains all 515 fine-mapped blocks.
#
# OFFICIAL SCAVENGE / gchromVAR LOGIC
#   1) peak-by-cell scATAC count matrix
#   2) fine-mapped posterior probabilities
#   3) add GC bias
#   4) 200 matched background peak iterations
#   5) weighted gchromVAR bias-corrected Z per cell
#   6) initial one-tailed P < 0.05 cells
#   7) if >5% eligible, keep top 5% Z as seed cells
#   8) scale factor = mean Z among top 1% cells
#
# THIS SCRIPT STOPS BEFORE NETWORK PROPAGATION.
# STEP11G3B will perform:
#   TF-IDF -> 30-d LSI -> Harmony donor correction -> mutual kNN (k=30)
#   -> random walk with restart gamma=0.05 -> TRS -> permutation inference.
#
# MEMORY STRATEGY
#   - stream each donor's dense TSV in small peak batches
#   - convert every batch immediately to sparse dgCMatrix
#   - save each donor sparse matrix as a checkpoint
#   - combine only after all 7 donor matrices pass QC
#
# IMPORTANT
#   Re-running this script reuses successful checkpoints.
# ==============================================================================

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = max(3600, getOption("timeout"))
)

set.seed(9527)

# ==============================================================================
# 0. FROZEN SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

EXPECTED_N_CELLS <- 11986L
EXPECTED_N_CELLTYPES <- 12L
EXPECTED_N_PEAKS <- 212084L
EXPECTED_N_BLOCKS <- 515L

PIP_MIN <- 0.001

BACKGROUND_ITERATIONS <- 200L

SEED_PERCENT <- 0.05
SCALE_PERCENT <- 0.01

# 500 was already validated in the earlier GSE238242 streaming workflow.
BATCH_N_PEAKS <- 500L

EXPECTED_DONORS <- c(
  "CF69",
  "CF77",
  "CF89",
  "CF91",
  "CF93",
  "CF97",
  "CF102"
)

# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

QC_ROOT <- file.path(
  ROOT,
  "00_QC"
)

EXTRACT_DIR <- file.path(
  ROOT,
  "01_EXTRACTED"
)

OBJECT_DIR <- file.path(
  ROOT,
  "02_OBJECTS"
)

META_RDS <- file.path(
  OBJECT_DIR,
  "GSE238242_author_metadata_LOCKED.rds"
)

PEAK_RDS <- file.path(
  OBJECT_DIR,
  "GSE238242_ATAC_peak_universe_GRCh38.rds"
)

PEAK_COORD_FILE <- file.path(
  OBJECT_DIR,
  "GSE238242_ATAC_peak_coordinates_GRCh38.csv"
)

FILE_MANIFEST <- file.path(
  QC_ROOT,
  "STEP11A_extracted_file_manifest.csv"
)

G2B_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "03_AF_SUSIE_FINEMAP"
)

G2B_READINESS <- file.path(
  G2B_ROOT,
  "00_QC",
  "STEP11G2B_readiness.csv"
)

G2B_STATUS <- file.path(
  G2B_ROOT,
  "00_QC",
  "STEP11G2B_block_status_FINAL.csv"
)

PIP_FILE <- file.path(
  G2B_ROOT,
  "03_FINAL_PIP",
  "AF_SCAVENGE_PIP_PRIMARY_ALL.tsv.gz"
)

OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "04_GCHROMVAR"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

MATRIX_DIR <- file.path(
  OUT_ROOT,
  "01_SCATAC_MATRIX"
)

DONOR_DIR <- file.path(
  MATRIX_DIR,
  "DONOR_CHECKPOINTS"
)

TRAIT_DIR <- file.path(
  OUT_ROOT,
  "02_AF_TRAIT_INPUT"
)

GCV_DIR <- file.path(
  OUT_ROOT,
  "03_GCHROMVAR"
)

FIG_DIR <- file.path(
  OUT_ROOT,
  "04_FIGURES"
)

LOG_DIR <- file.path(
  OUT_ROOT,
  "05_LOGS"
)

for (d in c(
  OUT_ROOT,
  QC_DIR,
  MATRIX_DIR,
  DONOR_DIR,
  TRAIT_DIR,
  GCV_DIR,
  FIG_DIR,
  LOG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

RSE_FILE <- file.path(
  MATRIX_DIR,
  "GSE238242_SCAVENGE_peak_by_cell_GRCh38.rds"
)

CELL_MANIFEST_FILE <- file.path(
  MATRIX_DIR,
  "GSE238242_SCAVENGE_cell_manifest.csv"
)

TRAIT_BED <- file.path(
  TRAIT_DIR,
  "AF_SCAVENGE_PIP_PRIMARY_GE0.001_GRCh38.bed"
)

TRAIT_TABLE <- file.path(
  TRAIT_DIR,
  "AF_SCAVENGE_PIP_PRIMARY_GE0.001_GRCh38.csv.gz"
)

SE_GC_FILE <- file.path(
  GCV_DIR,
  "GSE238242_SCAVENGE_SE_with_GCbias.rds"
)

TRAIT_IMPORT_FILE <- file.path(
  GCV_DIR,
  "AF_gchromVAR_trait_weights.rds"
)

BACKGROUND_FILE <- file.path(
  GCV_DIR,
  "GSE238242_gchromVAR_background200.rds"
)

DEV_FILE <- file.path(
  GCV_DIR,
  "AF_gchromVAR_weighted_deviations.rds"
)

CELL_SCORE_FILE <- file.path(
  GCV_DIR,
  "STEP11G3A_AF_cell_gchromVAR_seed_scores.csv.gz"
)

# ==============================================================================
# 2. PACKAGE INSTALLATION / LOADING
# ==============================================================================

install_cran_if_missing <- function(pkg) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "CRAN package installation failed: ",
      pkg
    )
  }
}

for (p in c(
  "data.table",
  "Matrix",
  "ggplot2",
  "remotes"
)) {
  install_cran_if_missing(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_pkgs <- c(
  "SummarizedExperiment",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "GenomeInfoDb",
  "BiocParallel",
  "chromVAR",
  "BSgenome.Hsapiens.UCSC.hg38"
)

missing_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc)) {

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

still_missing_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(still_missing_bioc)) {
  stop(
    paste0(
      "Bioconductor package installation failed:\n",
      paste(
        still_missing_bioc,
        collapse = "\n"
      )
    )
  )
}

# gchromVAR remains a GitHub package.
if (!requireNamespace("gchromVAR", quietly = TRUE)) {

  cat(
    "\ngchromVAR is not installed. Installing the official package from GitHub...\n"
  )

  try(
    remotes::install_github(
      "caleblareau/gchromVAR",
      dependencies = TRUE,
      upgrade = "never",
      build_vignettes = FALSE
    ),
    silent = FALSE
  )
}

if (!requireNamespace("gchromVAR", quietly = TRUE)) {
  stop(
    paste0(
      "gchromVAR could not be installed automatically.\n",
      "Install the official package with:\n",
      "remotes::install_github('caleblareau/gchromVAR')\n",
      "Then rerun this same STEP11G3A script."
    )
  )
}

library(data.table)
library(Matrix)
library(ggplot2)

suppressPackageStartupMessages(
  library(SummarizedExperiment)
)
suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(IRanges)
)
suppressPackageStartupMessages(
  library(S4Vectors)
)
suppressPackageStartupMessages(
  library(chromVAR)
)

# Windows-safe deterministic execution.
BiocParallel::register(
  BiocParallel::SerialParam()
)

HG38_GENOME <-
  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

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

clean_quotes <- function(x) {
  gsub(
    "\"",
    "",
    x,
    fixed = TRUE
  )
}

read_header_barcodes <- function(path) {

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  line <- readLines(
    con,
    n = 1L,
    warn = FALSE
  )

  if (!length(line)) {
    stop(
      "Empty ATAC file: ",
      path
    )
  }

  clean_quotes(
    strsplit(
      line,
      "\t",
      fixed = TRUE
    )[[1]]
  )
}

get_one_file <- function(
  manifest,
  donor,
  modality = "ATAC"
) {

  d <- donor
  m <- modality

  z <- manifest[
    donor == d &
      modality == m,
    full_path
  ]

  if (length(z) != 1L) {
    stop(
      "Expected exactly one ",
      modality,
      " file for ",
      donor,
      "; found ",
      length(z)
    )
  }

  if (!file.exists(z)) {
    stop(
      "Manifest file does not exist:\n",
      z
    )
  }

  norm_path(
    z,
    TRUE
  )
}

rbind_sparse_list <- function(x) {

  if (!length(x)) {
    stop(
      "Sparse block list is empty."
    )
  }

  out <- tryCatch(
    do.call(
      rbind,
      x
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {

    out <- Reduce(
      function(a, b) {
        rbind(
          a,
          b
        )
      },
      x
    )
  }

  methods::as(
    out,
    "dgCMatrix"
  )
}

cbind_sparse_list <- function(x) {

  if (!length(x)) {
    stop(
      "Sparse donor list is empty."
    )
  }

  out <- tryCatch(
    do.call(
      cbind,
      x
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {

    out <- Reduce(
      function(a, b) {
        cbind(
          a,
          b
        )
      },
      x
    )
  }

  methods::as(
    out,
    "dgCMatrix"
  )
}

# Exact SCAVENGE seed rule reproduced from seedindex():
# one-tailed P<0.05 first; if >5% then rank(-Z)<=5%N.
make_seed_index <- function(
  z_score,
  percent_cut = 0.05
) {

  if (
    percent_cut <= 0 ||
    percent_cut >= 1
  ) {
    stop(
      "percent_cut must be between 0 and 1."
    )
  }

  finite <- is.finite(
    z_score
  )

  p_one <- rep(
    NA_real_,
    length(z_score)
  )

  p_one[finite] <- pnorm(
    z_score[finite],
    lower.tail = FALSE
  )

  initial <- finite &
    p_one <= 0.05

  s_percent <- sum(
    initial
  ) *
    100 /
    length(z_score)

  if (
    s_percent >
      100 * percent_cut
  ) {

    r <- rank(
      -z_score,
      na.last = "keep"
    )

    seed <- finite &
      r <=
      (
        percent_cut *
          length(z_score)
      )

  } else {

    seed <- initial
  }

  seed[
    is.na(seed)
  ] <- FALSE

  list(
    seed = seed,
    initial = initial,
    p_one = p_one,
    initial_percent = s_percent
  )
}

# Exact SCAVENGE cal_scalefactor logic.
make_scale_factor <- function(
  z_score,
  percent_cut = 0.01
) {

  finite <- is.finite(
    z_score
  )

  z <- z_score[
    finite
  ]

  r <- rank(
    -z
  )

  idx <- r <=
    (
      percent_cut *
        length(z)
    )

  mean(
    z[idx]
  )
}

# ==============================================================================
# 4. HARD INPUT CHECK + STEP11G2B READINESS
# ==============================================================================

required_inputs <- c(
  metadata = META_RDS,
  peak_universe = PEAK_RDS,
  peak_coordinates = PEAK_COORD_FILE,
  file_manifest = FILE_MANIFEST,
  G2B_readiness = G2B_READINESS,
  G2B_status = G2B_STATUS,
  primary_PIP = PIP_FILE
)

missing_inputs <- required_inputs[
  !file.exists(
    required_inputs
  )
]

if (length(missing_inputs)) {
  stop(
    paste0(
      "Missing required STEP11G3A input(s):\n",
      paste(
        names(missing_inputs),
        missing_inputs,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

READY_G2B <- fread(
  G2B_READINESS
)

if (
  nrow(READY_G2B) != 12L ||
  !all(
    READY_G2B$pass %in%
      c(
        TRUE,
        "TRUE",
        1
      )
  )
) {
  print(
    READY_G2B
  )

  stop(
    "STEP11G2B readiness is not fully PASS."
  )
}

STATUS_G2B <- fread(
  G2B_STATUS
)

if (
  nrow(STATUS_G2B) !=
    EXPECTED_N_BLOCKS ||
  any(
    STATUS_G2B$status != "PASS"
  ) ||
  any(
    !STATUS_G2B$primary_converged
  )
) {
  stop(
    "STEP11G2B 515-block status is not fully PASS/converged."
  )
}

# Preserve low-reference-match blocks for PRIMARY analysis,
# but freeze their IDs for a future sensitivity analysis.
LOW_MATCH_BLOCKS <- STATUS_G2B[
  reference_match_fraction <
    0.90,
  .(
    selected_block_label,
    block_id,
    CHR,
    reference_match_fraction,
    n_final_variants
  )
]

fwrite(
  LOW_MATCH_BLOCKS,
  file.path(
    QC_DIR,
    "STEP11G3A_low_reference_match_blocks_for_sensitivity.csv"
  )
)

# ==============================================================================
# 5. LOAD FROZEN METADATA / PEAK UNIVERSE / MANIFEST
# ==============================================================================

META <- as.data.table(
  readRDS(
    META_RDS
  )
)

required_meta <- c(
  "barcode",
  "cell_type",
  "donor",
  "rhythm",
  "sex"
)

if (!all(
  required_meta %in%
    names(META)
)) {
  stop(
    paste0(
      "Locked metadata missing fields: ",
      paste(
        setdiff(
          required_meta,
          names(META)
        ),
        collapse = ", "
      )
    )
  )
}

META[
  ,
  `:=`(
    barcode =
      as.character(barcode),
    cell_type =
      as.character(cell_type),
    donor =
      as.character(donor),
    rhythm =
      as.character(rhythm),
    sex =
      as.character(sex)
  )
]

if (
  nrow(META) !=
    EXPECTED_N_CELLS
) {
  stop(
    "Expected 11,986 frozen cells; observed ",
    nrow(META)
  )
}

if (
  uniqueN(
    META$cell_type
  ) !=
    EXPECTED_N_CELLTYPES
) {
  stop(
    "Expected 12 frozen author cell types; observed ",
    uniqueN(
      META$cell_type
    )
  )
}

if (
  !setequal(
    unique(
      META$donor
    ),
    EXPECTED_DONORS
  )
) {
  stop(
    "Frozen donor set differs from expected 7 donors."
  )
}

PEAKS <- readRDS(
  PEAK_RDS
)

PEAKS <- as.character(
  PEAKS
)

if (
  length(PEAKS) !=
    EXPECTED_N_PEAKS
) {
  stop(
    "Expected 212,084 frozen ATAC peaks; observed ",
    length(PEAKS)
  )
}

if (anyDuplicated(PEAKS)) {
  stop(
    "Frozen ATAC peak universe contains duplicated IDs."
  )
}

PEAK_COORDS <- fread(
  PEAK_COORD_FILE
)

if (
  nrow(PEAK_COORDS) !=
    EXPECTED_N_PEAKS
) {
  stop(
    "Peak coordinate table row count mismatch."
  )
}

required_peak_cols <- c(
  "peak_id",
  "chr",
  "start",
  "end"
)

if (!all(
  required_peak_cols %in%
    names(PEAK_COORDS)
)) {
  stop(
    "Frozen peak coordinate table lacks peak_id/chr/start/end."
  )
}

if (!identical(
  as.character(
    PEAK_COORDS$peak_id
  ),
  PEAKS
)) {
  stop(
    "Frozen peak coordinate order does not match peak universe RDS."
  )
}

if (
  anyNA(
    PEAK_COORDS$start
  ) ||
  anyNA(
    PEAK_COORDS$end
  ) ||
  any(
    PEAK_COORDS$start < 1L
  ) ||
  any(
    PEAK_COORDS$end <
      PEAK_COORDS$start
  )
) {
  stop(
    "Invalid GRCh38 peak coordinates for RangedSummarizedExperiment."
  )
}

MANIFEST <- fread(
  FILE_MANIFEST
)

required_manifest_cols <- c(
  "donor",
  "modality",
  "full_path"
)

if (!all(
  required_manifest_cols %in%
    names(MANIFEST)
)) {
  stop(
    "STEP11A file manifest lacks donor/modality/full_path."
  )
}

# ==============================================================================
# 6. BUILD / REUSE DONOR-LEVEL SPARSE scATAC MATRICES
# ==============================================================================

DONOR_QC_LIST <- list()

for (d in EXPECTED_DONORS) {

  DONOR_MATRIX_FILE <- file.path(
    DONOR_DIR,
    paste0(
      "GSE238242_",
      d,
      "_ATAC_sparse_GRCh38.rds"
    )
  )

  DONOR_META_FILE <- file.path(
    DONOR_DIR,
    paste0(
      "GSE238242_",
      d,
      "_cell_manifest.csv"
    )
  )

  # ----------------------------------------------------------------
  # Reuse a complete donor checkpoint.
  # ----------------------------------------------------------------

  checkpoint_ok <- FALSE

  if (
    file.exists(
      DONOR_MATRIX_FILE
    ) &&
    file.exists(
      DONOR_META_FILE
    )
  ) {

    dm <- tryCatch(
      readRDS(
        DONOR_MATRIX_FILE
      ),
      error = function(e) NULL
    )

    dmeta <- tryCatch(
      fread(
        DONOR_META_FILE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(dm) &&
      !is.null(dmeta) &&
      nrow(dm) ==
        EXPECTED_N_PEAKS &&
      ncol(dm) ==
        nrow(dmeta) &&
      identical(
        rownames(dm),
        PEAKS
      ) &&
      identical(
        colnames(dm),
        dmeta$cell_id
      )
    ) {
      checkpoint_ok <- TRUE
    }

    rm(
      dm,
      dmeta
    )

    invisible(
      gc(
        verbose = FALSE
      )
    )
  }

  if (!checkpoint_ok) {

    cat(
      "\n============================================================\n",
      "BUILDING SPARSE ATAC — ",
      d,
      "\n",
      "============================================================\n",
      sep = ""
    )

    ATAC_FILE <- get_one_file(
      MANIFEST,
      d,
      "ATAC"
    )

    BARCODES <- read_header_barcodes(
      ATAC_FILE
    )

    META_D <- META[
      donor == d
    ]

    if (
      !setequal(
        BARCODES,
        META_D$barcode
      )
    ) {
      stop(
        "ATAC / metadata barcode set mismatch for donor ",
        d
      )
    }

    if (
      anyDuplicated(BARCODES) ||
      anyDuplicated(
        META_D$barcode
      )
    ) {
      stop(
        "Duplicated donor-level barcode detected for ",
        d
      )
    }

    idx_meta <- match(
      BARCODES,
      META_D$barcode
    )

    if (anyNA(idx_meta)) {
      stop(
        "Metadata alignment failed for donor ",
        d
      )
    }

    META_ORDERED <- copy(
      META_D[
        idx_meta
      ]
    )

    META_ORDERED[
      ,
      cell_id :=
        paste(
          donor,
          barcode,
          sep = "#"
        )
    ]

    if (anyDuplicated(
      META_ORDERED$cell_id
    )) {
      stop(
        "Duplicated global cell_id within donor ",
        d
      )
    }

    con <- gzfile(
      ATAC_FILE,
      open = "rt"
    )

    # Skip header.
    invisible(
      readLines(
        con,
        n = 1L,
        warn = FALSE
      )
    )

    SPARSE_BLOCKS <- list()

    row_start <- 1L
    n_rows_seen <- 0L
    batch_id <- 0L
    nnz_total <- 0

    repeat {

      lines <- readLines(
        con,
        n = BATCH_N_PEAKS,
        warn = FALSE
      )

      if (!length(lines)) {
        break
      }

      batch_id <- batch_id + 1L

      DT <- fread(
        text = paste(
          lines,
          collapse = "\n"
        ),
        sep = "\t",
        header = FALSE,
        quote = "\"",
        data.table = TRUE,
        showProgress = FALSE
      )

      n_batch <- nrow(
        DT
      )

      row_end <- row_start +
        n_batch -
        1L

      if (
        row_end >
          EXPECTED_N_PEAKS
      ) {

        close(
          con
        )

        stop(
          "ATAC file contains more peaks than frozen universe for ",
          d
        )
      }

      peak_batch <- as.character(
        DT[[1]]
      )

      if (!identical(
        peak_batch,
        PEAKS[
          row_start:row_end
        ]
      )) {

        close(
          con
        )

        stop(
          "ATAC peak order mismatch for ",
          d,
          " at rows ",
          row_start,
          "-",
          row_end
        )
      }

      X <- as.matrix(
        DT[
          ,
          -1,
          with = FALSE
        ]
      )

      storage.mode(
        X
      ) <- "double"

      if (
        ncol(X) !=
          length(
            BARCODES
          )
      ) {

        close(
          con
        )

        stop(
          "ATAC cell count mismatch for ",
          d,
          " at batch ",
          batch_id
        )
      }

      SX <- Matrix::Matrix(
        X,
        sparse = TRUE
      )

      SX <- methods::as(
        SX,
        "dgCMatrix"
      )

      nnz_total <- nnz_total +
        length(
          SX@x
        )

      SPARSE_BLOCKS[[
        length(
          SPARSE_BLOCKS
        ) + 1L
      ]] <- SX

      rm(
        DT,
        X,
        SX
      )

      row_start <- row_end +
        1L

      n_rows_seen <-
        n_rows_seen +
        n_batch

      if (
        batch_id %% 25L ==
          0L
      ) {

        cat(
          "  ",
          d,
          ": processed ",
          format(
            n_rows_seen,
            big.mark = ","
          ),
          " / ",
          format(
            EXPECTED_N_PEAKS,
            big.mark = ","
          ),
          " peaks\n",
          sep = ""
        )

        invisible(
          gc(
            verbose = FALSE
          )
        )
      }
    }

    close(
      con
    )

    if (
      n_rows_seen !=
        EXPECTED_N_PEAKS
    ) {
      stop(
        "ATAC peak count mismatch for ",
        d,
        ": ",
        n_rows_seen,
        " vs ",
        EXPECTED_N_PEAKS
      )
    }

    DONOR_MAT <- rbind_sparse_list(
      SPARSE_BLOCKS
    )

    rm(
      SPARSE_BLOCKS
    )

    rownames(
      DONOR_MAT
    ) <- PEAKS

    colnames(
      DONOR_MAT
    ) <- META_ORDERED$cell_id

    if (
      nrow(DONOR_MAT) !=
        EXPECTED_N_PEAKS ||
      ncol(DONOR_MAT) !=
        nrow(
          META_ORDERED
        )
    ) {
      stop(
        "Sparse donor matrix dimension mismatch for ",
        d
      )
    }

    saveRDS(
      DONOR_MAT,
      DONOR_MATRIX_FILE,
      compress = FALSE
    )

    fwrite(
      META_ORDERED,
      DONOR_META_FILE
    )

    cat(
      "Saved donor checkpoint: ",
      d,
      " | ",
      nrow(DONOR_MAT),
      " peaks x ",
      ncol(DONOR_MAT),
      " cells | nnz=",
      format(
        length(
          DONOR_MAT@x
        ),
        big.mark = ","
      ),
      "\n",
      sep = ""
    )

    rm(
      DONOR_MAT,
      META_ORDERED
    )

    invisible(
      gc(
        verbose = FALSE
      )
    )
  }

  # Audit the final checkpoint.
  DM <- readRDS(
    DONOR_MATRIX_FILE
  )

  DMC <- fread(
    DONOR_META_FILE
  )

  DONOR_QC_LIST[[
    length(
      DONOR_QC_LIST
    ) + 1L
  ]] <- data.table(
    donor = d,
    n_peaks = nrow(DM),
    n_cells = ncol(DM),
    nnz = length(DM@x),
    peak_order_match =
      identical(
        rownames(DM),
        PEAKS
      ),
    cell_id_order_match =
      identical(
        colnames(DM),
        DMC$cell_id
      ),
    matrix_file =
      norm_path(
        DONOR_MATRIX_FILE,
        TRUE
      )
  )

  rm(
    DM,
    DMC
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )
}

DONOR_QC <- rbindlist(
  DONOR_QC_LIST
)

fwrite(
  DONOR_QC,
  file.path(
    QC_DIR,
    "STEP11G3A_donor_sparse_matrix_QC.csv"
  )
)

if (
  !all(
    DONOR_QC$n_peaks ==
      EXPECTED_N_PEAKS
  ) ||
  sum(
    DONOR_QC$n_cells
  ) !=
    EXPECTED_N_CELLS ||
  !all(
    DONOR_QC$peak_order_match
  ) ||
  !all(
    DONOR_QC$cell_id_order_match
  )
) {
  stop(
    "Donor sparse matrix QC failed."
  )
}

# ==============================================================================
# 7. COMBINE / REUSE FULL PEAK-BY-CELL RANGEDSUMMARIZEDEXPERIMENT
# ==============================================================================

RSE_VALID <- FALSE

if (file.exists(
  RSE_FILE
)) {

  SE <- tryCatch(
    readRDS(
      RSE_FILE
    ),
    error = function(e) NULL
  )

  if (
    !is.null(SE) &&
    nrow(SE) ==
      EXPECTED_N_PEAKS &&
    ncol(SE) ==
      EXPECTED_N_CELLS
  ) {
    RSE_VALID <- TRUE
  }
}

if (!RSE_VALID) {

  cat(
    "\nCombining 7 donor sparse matrices...\n"
  )

  DONOR_MATS <- list()
  CELL_META_LIST <- list()

  for (d in EXPECTED_DONORS) {

    dm_file <- file.path(
      DONOR_DIR,
      paste0(
        "GSE238242_",
        d,
        "_ATAC_sparse_GRCh38.rds"
      )
    )

    dc_file <- file.path(
      DONOR_DIR,
      paste0(
        "GSE238242_",
        d,
        "_cell_manifest.csv"
      )
    )

    DONOR_MATS[[
      d
    ]] <- readRDS(
      dm_file
    )

    CELL_META_LIST[[
      d
    ]] <- fread(
      dc_file
    )
  }

  COUNTS <- cbind_sparse_list(
    DONOR_MATS
  )

  CELL_META <- rbindlist(
    CELL_META_LIST,
    use.names = TRUE,
    fill = TRUE
  )

  if (
    nrow(COUNTS) !=
      EXPECTED_N_PEAKS ||
    ncol(COUNTS) !=
      EXPECTED_N_CELLS ||
    nrow(CELL_META) !=
      EXPECTED_N_CELLS
  ) {
    stop(
      "Combined peak-by-cell dimensions do not match frozen design."
    )
  }

  if (!identical(
    colnames(COUNTS),
    CELL_META$cell_id
  )) {
    stop(
      "Combined count matrix and cell metadata order mismatch."
    )
  }

  if (!identical(
    rownames(COUNTS),
    PEAKS
  )) {
    stop(
      "Combined count matrix peak order mismatch."
    )
  }

  DEPTH <- Matrix::colSums(
    COUNTS
  )

  N_ACCESSIBLE <- diff(
    COUNTS@p
  )

  CELL_META[
    ,
    `:=`(
      total_ATAC_counts =
        as.numeric(DEPTH),
      n_accessible_peaks =
        as.integer(
          N_ACCESSIBLE
        )
    )
  ]

  if (
    any(
      !is.finite(
        CELL_META$total_ATAC_counts
      )
    ) ||
    any(
      CELL_META$total_ATAC_counts <=
        0
    )
  ) {
    stop(
      "At least one frozen cell has zero/invalid ATAC depth. ",
      "Do not silently filter cells."
    )
  }

  PEAK_GR <- GRanges(
    seqnames =
      as.character(
        PEAK_COORDS$chr
      ),
    ranges = IRanges(
      start =
        as.integer(
          PEAK_COORDS$start
        ),
      end =
        as.integer(
          PEAK_COORDS$end
        )
    )
  )

  names(
    PEAK_GR
  ) <- PEAKS

  CD <- S4Vectors::DataFrame(
    CELL_META
  )

  rownames(
    CD
  ) <- CELL_META$cell_id

  SE <- SummarizedExperiment(
    assays = list(
      counts = COUNTS
    ),
    rowRanges = PEAK_GR,
    colData = CD
  )

  if (
    nrow(SE) !=
      EXPECTED_N_PEAKS ||
    ncol(SE) !=
      EXPECTED_N_CELLS
  ) {
    stop(
      "Final RangedSummarizedExperiment dimensions are incorrect."
    )
  }

  saveRDS(
    SE,
    RSE_FILE,
    compress = FALSE
  )

  fwrite(
    CELL_META,
    CELL_MANIFEST_FILE
  )

  rm(
    DONOR_MATS,
    CELL_META_LIST,
    COUNTS,
    CD
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )

} else {

  CELL_META <- as.data.table(
    as.data.frame(
      SummarizedExperiment::colData(
        SE
      )
    )
  )

  if (
    !"cell_id" %in%
      names(
        CELL_META
      )
  ) {
    CELL_META[
      ,
      cell_id :=
        colnames(
          SE
        )
    ]
  }
}

# Final matrix QC.
COUNTS <- SummarizedExperiment::assay(
  SE,
  "counts"
)

if (!inherits(
  COUNTS,
  "sparseMatrix"
)) {
  stop(
    "SCAVENGE count matrix is not sparse."
  )
}

TOTAL_NNZ <- length(
  methods::as(
    COUNTS,
    "dgCMatrix"
  )@x
)

MATRIX_QC <- data.table(
  metric = c(
    "n_peaks",
    "n_cells",
    "n_cell_types",
    "n_donors",
    "n_nonzero_matrix_entries",
    "min_cell_depth",
    "median_cell_depth",
    "max_cell_depth"
  ),
  value = c(
    nrow(SE),
    ncol(SE),
    uniqueN(
      as.character(
        SummarizedExperiment::colData(
          SE
        )$cell_type
      )
    ),
    uniqueN(
      as.character(
        SummarizedExperiment::colData(
          SE
        )$donor
      )
    ),
    TOTAL_NNZ,
    min(
      SummarizedExperiment::colData(
        SE
      )$total_ATAC_counts
    ),
    median(
      SummarizedExperiment::colData(
        SE
      )$total_ATAC_counts
    ),
    max(
      SummarizedExperiment::colData(
        SE
      )$total_ATAC_counts
    )
  )
)

fwrite(
  MATRIX_QC,
  file.path(
    QC_DIR,
    "STEP11G3A_full_scATAC_matrix_QC.csv"
  )
)

# ==============================================================================
# 8. PREPARE PRIMARY AF FINE-MAPPED PIP >= 0.001 IN GRCh38
# ==============================================================================

cat(
  "\nReading primary AF fine-mapping PIP mass...\n"
)

PIP_MASS_DT <- fread(
  PIP_FILE,
  select = "PIP_primary",
  showProgress = TRUE
)

TOTAL_PIP_MASS <- sum(
  PIP_MASS_DT$PIP_primary,
  na.rm = TRUE
)

rm(
  PIP_MASS_DT
)

invisible(
  gc(
    verbose = FALSE
  )
)

PIP <- fread(
  PIP_FILE,
  select = c(
    "selected_block_label",
    "CHR_HG19",
    "BP_GRCh38",
    "SNP",
    "PIP_primary"
  ),
  showProgress = TRUE
)

PIP <- PIP[
  is.finite(
    PIP_primary
  ) &
    PIP_primary >=
      PIP_MIN &
    is.finite(
      BP_GRCh38
    ) &
    BP_GRCh38 >
      0 &
    CHR_HG19 %between%
      c(
        1L,
        22L
      )
]

PIP[
  ,
  `:=`(
    CHR_GRCh38 =
      as.integer(
        CHR_HG19
      ),
    chr =
      paste0(
        "chr",
        as.integer(
          CHR_HG19
        )
      ),
    start0 =
      as.integer(
        BP_GRCh38
      ) -
      1L,
    end =
      as.integer(
        BP_GRCh38
      )
  )
]

setorder(
  PIP,
  CHR_GRCh38,
  BP_GRCh38,
  -PIP_primary
)

RETAINED_PIP_MASS <- sum(
  PIP$PIP_primary
)

RETAINED_MASS_FRACTION <-
  RETAINED_PIP_MASS /
  TOTAL_PIP_MASS

N_BLOCKS_PIP <- uniqueN(
  PIP$selected_block_label
)

if (
  N_BLOCKS_PIP !=
    EXPECTED_N_BLOCKS
) {
  stop(
    "PIP>=0.001 does not retain all 515 fine-mapped blocks; observed ",
    N_BLOCKS_PIP
  )
}

if (
  !is.finite(
    RETAINED_MASS_FRACTION
  ) ||
  RETAINED_MASS_FRACTION <
    0.95
) {
  stop(
    "PIP>=0.001 retains <95% of total fine-mapping PIP mass."
  )
}

fwrite(
  PIP[
    ,
    .(
      selected_block_label,
      CHR_GRCh38,
      BP_GRCh38,
      SNP,
      PIP_primary
    )
  ],
  TRAIT_TABLE,
  compress = "gzip"
)

# Official gchromVAR-style BED:
# column1 chr
# column2 0-based start
# column3 1-based/half-open end
# column4 locus/block label
# column5 posterior probability
fwrite(
  PIP[
    ,
    .(
      chr,
      start0,
      end,
      selected_block_label,
      PIP_primary
    )
  ],
  TRAIT_BED,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE
)

PIP_QC <- data.table(
  metric = c(
    "PIP_threshold",
    "n_variants_retained",
    "n_blocks_retained",
    "total_primary_PIP_mass",
    "retained_primary_PIP_mass",
    "retained_PIP_mass_fraction"
  ),
  value = c(
    PIP_MIN,
    nrow(PIP),
    N_BLOCKS_PIP,
    TOTAL_PIP_MASS,
    RETAINED_PIP_MASS,
    RETAINED_MASS_FRACTION
  )
)

fwrite(
  PIP_QC,
  file.path(
    QC_DIR,
    "STEP11G3A_AF_PIP_input_QC.csv"
  )
)

# ==============================================================================
# 9. PIP -> ATAC PEAK OVERLAP: MANUAL AUDIT + OFFICIAL importBedScore
# ==============================================================================

PEAK_GR <- SummarizedExperiment::rowRanges(
  SE
)

PIP_GR <- GRanges(
  seqnames =
    PIP$chr,
  ranges = IRanges(
    start =
      as.integer(
        PIP$BP_GRCh38
      ),
    width = 1L
  )
)

hits <- findOverlaps(
  PEAK_GR,
  PIP_GR,
  ignore.strand = TRUE
)

N_VARIANTS_OVERLAP <- uniqueN(
  subjectHits(
    hits
  )
)

N_PEAKS_WEIGHTED <- uniqueN(
  queryHits(
    hits
  )
)

MANUAL_WEIGHT <- numeric(
  nrow(SE)
)

if (length(hits)) {

  tmp_weight <- data.table(
    peak_index =
      queryHits(
        hits
      ),
    weight =
      PIP$PIP_primary[
        subjectHits(
          hits
        )
      ]
  )[
    ,
    .(
      weight =
        sum(
          weight
        )
    ),
    by = peak_index
  ]

  MANUAL_WEIGHT[
    tmp_weight$peak_index
  ] <- tmp_weight$weight
}

if (
  sum(
    MANUAL_WEIGHT
  ) <=
    0
) {
  stop(
    "None of the retained AF fine-mapped PIP overlaps the GSE238242 peak universe."
  )
}

TRAIT_IMPORT <- NULL
IMPORT_METHOD <- NA_character_
IMPORT_CONCORDANCE <- NA_real_
IMPORT_MAX_ABS_DIFF <- NA_real_

official_import <- tryCatch(
  gchromVAR::importBedScore(
    ranges = PEAK_GR,
    files = TRAIT_BED,
    colidx = 5,
    FUN = sum,
    default.val = 0
  ),
  error = function(e) {
    message(
      "Official gchromVAR::importBedScore failed: ",
      conditionMessage(e)
    )
    NULL
  }
)

if (!is.null(
  official_import
)) {

  official_w <- as.numeric(
    SummarizedExperiment::assay(
      official_import
    )[
      ,
      1
    ]
  )

  if (
    length(
      official_w
    ) ==
      length(
        MANUAL_WEIGHT
      )
  ) {

    IMPORT_MAX_ABS_DIFF <- max(
      abs(
        official_w -
          MANUAL_WEIGHT
      ),
      na.rm = TRUE
    )

    if (
      stats::sd(
        official_w
      ) >
        0 &&
      stats::sd(
        MANUAL_WEIGHT
      ) >
        0
    ) {

      IMPORT_CONCORDANCE <- suppressWarnings(
        cor(
          official_w,
          MANUAL_WEIGHT,
          method = "pearson"
        )
      )
    }

    # Use the official implementation when it agrees with direct
    # 1-bp point-overlap arithmetic.
    if (
      is.finite(
        IMPORT_MAX_ABS_DIFF
      ) &&
      IMPORT_MAX_ABS_DIFF <
        1e-8
    ) {

      TRAIT_IMPORT <- official_import
      IMPORT_METHOD <-
        "official_gchromVAR_importBedScore"

    }
  }
}

# Conservative fallback:
# direct GRanges overlap with exact GRCh38 SNP coordinate.
if (is.null(
  TRAIT_IMPORT
)) {

  WEIGHT_MAT <- Matrix::Matrix(
    MANUAL_WEIGHT,
    ncol = 1L,
    sparse = TRUE
  )

  colnames(
    WEIGHT_MAT
  ) <- "AF"

  TRAIT_IMPORT <- SummarizedExperiment(
    assays = list(
      counts = WEIGHT_MAT
    ),
    rowRanges = PEAK_GR
  )

  IMPORT_METHOD <-
    "direct_GRanges_PIP_sum_fallback"

  warning(
    paste0(
      "Using direct GRanges overlap weights because the installed ",
      "gchromVAR importBedScore result was unavailable or differed from ",
      "the direct coordinate audit. This is recorded in provenance."
    )
  )
}

saveRDS(
  TRAIT_IMPORT,
  TRAIT_IMPORT_FILE,
  compress = TRUE
)

OVERLAP_QC <- data.table(
  metric = c(
    "n_PIP_variants",
    "n_PIP_variants_overlapping_at_least_one_peak",
    "variant_overlap_fraction",
    "n_ATAC_peaks_with_nonzero_AF_PIP",
    "sum_PIP_weight_over_peaks",
    "weight_import_method",
    "official_vs_manual_pearson",
    "official_vs_manual_max_abs_diff"
  ),
  value = c(
    nrow(PIP),
    N_VARIANTS_OVERLAP,
    N_VARIANTS_OVERLAP /
      nrow(PIP),
    N_PEAKS_WEIGHTED,
    sum(
      MANUAL_WEIGHT
    ),
    IMPORT_METHOD,
    IMPORT_CONCORDANCE,
    IMPORT_MAX_ABS_DIFF
  )
)

fwrite(
  OVERLAP_QC,
  file.path(
    QC_DIR,
    "STEP11G3A_PIP_peak_overlap_QC.csv"
  )
)

# ==============================================================================
# 10. ADD GC BIAS — CHECKPOINT
# ==============================================================================

if (file.exists(
  SE_GC_FILE
)) {

  SE_GC <- readRDS(
    SE_GC_FILE
  )

  gc_valid <- (
    nrow(SE_GC) ==
      EXPECTED_N_PEAKS &&
      ncol(SE_GC) ==
        EXPECTED_N_CELLS &&
      "bias" %in%
        names(
          S4Vectors::mcols(
            SummarizedExperiment::rowRanges(
              SE_GC
            )
          )
        )
  )

} else {

  gc_valid <- FALSE
}

if (!gc_valid) {

  cat(
    "\nAdding GC bias using BSgenome.Hsapiens.UCSC.hg38...\n"
  )

  SE_GC <- chromVAR::addGCBias(
    object = SE,
    genome = HG38_GENOME
  )

  saveRDS(
    SE_GC,
    SE_GC_FILE,
    compress = FALSE
  )
}

# ==============================================================================
# 11. 200 MATCHED BACKGROUND PEAK ITERATIONS — CHECKPOINT
# ==============================================================================

if (file.exists(
  BACKGROUND_FILE
)) {

  BG <- readRDS(
    BACKGROUND_FILE
  )

  bg_valid <- (
    is.matrix(BG) ||
      inherits(
        BG,
        "Matrix"
      )
  ) &&
    nrow(BG) ==
      EXPECTED_N_PEAKS &&
    ncol(BG) ==
      BACKGROUND_ITERATIONS

} else {

  bg_valid <- FALSE
}

if (!bg_valid) {

  cat(
    "\nComputing ",
    BACKGROUND_ITERATIONS,
    " chromVAR matched-background iterations...\n",
    sep = ""
  )

  set.seed(
    9527
  )

  BG <- chromVAR::getBackgroundPeaks(
    object = SE_GC,
    niterations =
      BACKGROUND_ITERATIONS
  )

  saveRDS(
    BG,
    BACKGROUND_FILE,
    compress = TRUE
  )
}

# ==============================================================================
# 12. WEIGHTED gchromVAR DEVIATIONS — CHECKPOINT
# ==============================================================================

if (file.exists(
  DEV_FILE
)) {

  DEV <- readRDS(
    DEV_FILE
  )

  dev_valid <- (
    "z" %in%
      SummarizedExperiment::assayNames(
        DEV
      )
  )

} else {

  dev_valid <- FALSE
}

if (!dev_valid) {

  cat(
    "\nRunning weighted gchromVAR for AF fine-mapped PIP...\n"
  )

  # IMPORTANT compatibility fix:
  # gchromVAR::computeWeightedDeviations() accepts a numeric/Matrix weight
  # matrix directly.  Passing a SummarizedExperiment makes chromVAR call
  # annotationMatches(), which rejects a fallback assay named "counts".
  # Therefore always extract the first trait-weight assay and pass the
  # peak x trait Matrix directly.
  if (inherits(TRAIT_IMPORT, "SummarizedExperiment")) {

    TRAIT_WEIGHT_MAT <- SummarizedExperiment::assay(
      TRAIT_IMPORT,
      1
    )

  } else if (
    inherits(TRAIT_IMPORT, "Matrix") ||
    is.matrix(TRAIT_IMPORT)
  ) {

    TRAIT_WEIGHT_MAT <- TRAIT_IMPORT

  } else {

    stop(
      "Unsupported TRAIT_IMPORT class: ",
      paste(class(TRAIT_IMPORT), collapse = ", ")
    )
  }

  TRAIT_WEIGHT_MAT <- methods::as(
    TRAIT_WEIGHT_MAT,
    "dgCMatrix"
  )

  if (
    nrow(TRAIT_WEIGHT_MAT) != nrow(SE_GC) ||
    ncol(TRAIT_WEIGHT_MAT) != 1L
  ) {
    stop(
      paste0(
        "Trait-weight matrix dimension mismatch: ",
        nrow(TRAIT_WEIGHT_MAT), " x ", ncol(TRAIT_WEIGHT_MAT),
        "; expected ", nrow(SE_GC), " x 1."
      )
    )
  }

  if (sum(TRAIT_WEIGHT_MAT) <= 0) {
    stop(
      "Trait-weight matrix contains no positive AF fine-mapping weight."
    )
  }

  colnames(TRAIT_WEIGHT_MAT) <- "AF"

  cat(
    "Trait weight matrix ready: ",
    nrow(TRAIT_WEIGHT_MAT),
    " peaks x ",
    ncol(TRAIT_WEIGHT_MAT),
    " trait; nonzero peaks = ",
    Matrix::nnzero(TRAIT_WEIGHT_MAT),
    "\n",
    sep = ""
  )

  DEV <- gchromVAR::computeWeightedDeviations(
    object = SE_GC,
    weights = TRAIT_WEIGHT_MAT,
    background_peaks = BG
  )

  if (!"z" %in%
      SummarizedExperiment::assayNames(
        DEV
      )) {
    stop(
      "gchromVAR output does not contain the expected 'z' assay."
    )
  }

  saveRDS(
    DEV,
    DEV_FILE,
    compress = TRUE
  )
}

# ==============================================================================
# 13. EXTRACT CELL-LEVEL gchromVAR Z SCORE
# ==============================================================================

ZMAT <- SummarizedExperiment::assay(
  DEV,
  "z"
)

if (
  ncol(ZMAT) ==
    EXPECTED_N_CELLS
) {

  if (nrow(ZMAT) != 1L) {
    stop(
      "Expected one AF trait row in gchromVAR Z assay; observed ",
      nrow(ZMAT)
    )
  }

  Z <- as.numeric(
    ZMAT[
      1,
      ]
  )

} else if (
  nrow(ZMAT) ==
    EXPECTED_N_CELLS
) {

  if (ncol(ZMAT) != 1L) {
    stop(
      "Expected one AF trait column in gchromVAR Z assay."
    )
  }

  Z <- as.numeric(
    ZMAT[
      ,
      1
    ]
  )

} else {

  stop(
    paste0(
      "Unexpected gchromVAR Z dimensions: ",
      paste(
        dim(ZMAT),
        collapse = " x "
      )
    )
  )
}

N_FINITE_Z <- sum(
  is.finite(
    Z
  )
)

FINITE_Z_FRACTION <-
  N_FINITE_Z /
  EXPECTED_N_CELLS

if (
  FINITE_Z_FRACTION <
    0.99
) {
  stop(
    "Fewer than 99% of frozen cells received a finite gchromVAR Z score."
  )
}

# ==============================================================================
# 14. OFFICIAL SCAVENGE SEED RULE + SCALE FACTOR
# ==============================================================================

SEED_RES <- make_seed_index(
  Z,
  percent_cut =
    SEED_PERCENT
)

SCALE_FACTOR <- make_scale_factor(
  Z,
  percent_cut =
    SCALE_PERCENT
)

SEED <- SEED_RES$seed

if (
  sum(SEED) <= 0L
) {
  stop(
    "No SCAVENGE seed cells were selected."
  )
}

if (
  sum(SEED) >
    ceiling(
      SEED_PERCENT *
        EXPECTED_N_CELLS
    )
) {
  stop(
    "Seed-cell fraction unexpectedly exceeds the frozen 5% cap."
  )
}

# ==============================================================================
# 15. CELL-LEVEL RESULT TABLE
# ==============================================================================

CD <- as.data.table(
  as.data.frame(
    SummarizedExperiment::colData(
      SE_GC
    )
  )
)

if (
  nrow(CD) !=
    EXPECTED_N_CELLS
) {
  stop(
    "colData cell count mismatch."
  )
}

if (!"cell_id" %in%
    names(CD)) {
  CD[
    ,
    cell_id :=
      colnames(
        SE_GC
      )
  ]
}

CELL_SCORE <- copy(
  CD
)

CELL_SCORE[
  ,
  `:=`(
    AF_gchromVAR_Z = Z,
    AF_gchromVAR_one_tailed_P =
      SEED_RES$p_one,
    AF_initial_P_lt_0_05 =
      SEED_RES$initial,
    AF_SCAVENGE_seed =
      SEED
  )
]

CELL_SCORE[
  ,
  AF_gchromVAR_rank :=
    frank(
      -AF_gchromVAR_Z,
      ties.method = "average",
      na.last = "keep"
    )
]

CELL_SCORE[
  ,
  AF_gchromVAR_percentile :=
    1 -
    (
      AF_gchromVAR_rank -
        1
    ) /
    (
      .N -
        1
    )
]

setorder(
  CELL_SCORE,
  AF_gchromVAR_rank
)

fwrite(
  CELL_SCORE,
  CELL_SCORE_FILE,
  compress = "gzip"
)

# ==============================================================================
# 16. CELL-TYPE + DONOR SUMMARY
# ==============================================================================

CT_SUMMARY <- CELL_SCORE[
  ,
  .(
    n_cells = .N,
    n_seed =
      sum(
        AF_SCAVENGE_seed
      ),
    seed_fraction =
      mean(
        AF_SCAVENGE_seed
      ),
    mean_Z =
      mean(
        AF_gchromVAR_Z,
        na.rm = TRUE
      ),
    median_Z =
      median(
        AF_gchromVAR_Z,
        na.rm = TRUE
      ),
    q75_Z =
      quantile(
        AF_gchromVAR_Z,
        0.75,
        na.rm = TRUE
      ),
    q90_Z =
      quantile(
        AF_gchromVAR_Z,
        0.90,
        na.rm = TRUE
      ),
    q95_Z =
      quantile(
        AF_gchromVAR_Z,
        0.95,
        na.rm = TRUE
      )
  ),
  by = cell_type
]

setorder(
  CT_SUMMARY,
  -median_Z
)

fwrite(
  CT_SUMMARY,
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_gchromVAR_celltype_summary.csv"
  )
)

DONOR_CT_SUMMARY <- CELL_SCORE[
  ,
  .(
    n_cells = .N,
    n_seed =
      sum(
        AF_SCAVENGE_seed
      ),
    seed_fraction =
      mean(
        AF_SCAVENGE_seed
      ),
    mean_Z =
      mean(
        AF_gchromVAR_Z,
        na.rm = TRUE
      ),
    median_Z =
      median(
        AF_gchromVAR_Z,
        na.rm = TRUE
      )
  ),
  by = .(
    donor,
    rhythm,
    cell_type
  )
]

fwrite(
  DONOR_CT_SUMMARY,
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_gchromVAR_donor_celltype_summary.csv"
  )
)

SEED_SUMMARY <- data.table(
  metric = c(
    "n_cells_total",
    "n_finite_Z",
    "finite_Z_fraction",
    "n_initial_one_tailed_P_lt_0.05",
    "initial_enriched_fraction",
    "seed_percent_cap",
    "n_final_seed_cells",
    "final_seed_fraction",
    "scale_percent",
    "scale_factor"
  ),
  value = c(
    EXPECTED_N_CELLS,
    N_FINITE_Z,
    FINITE_Z_FRACTION,
    sum(
      SEED_RES$initial
    ),
    mean(
      SEED_RES$initial
    ),
    SEED_PERCENT,
    sum(
      SEED
    ),
    mean(
      SEED
    ),
    SCALE_PERCENT,
    SCALE_FACTOR
  )
)

fwrite(
  SEED_SUMMARY,
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_SCAVENGE_seed_summary.csv"
  )
)

# ==============================================================================
# 17. PUBLICATION-QUALITY EXPLORATORY FIGURES
#     Descriptive only; formal inference follows network propagation/permutation.
# ==============================================================================

CELL_SCORE[
  ,
  cell_type_plot :=
    factor(
      cell_type,
      levels =
        CT_SUMMARY$cell_type
    )
]

# Muted journal-style palette, one stable color per author cell type.
CT_LEVELS <- levels(
  CELL_SCORE$cell_type_plot
)

BASE_PALETTE <- c(
  "#3C5488",
  "#E64B35",
  "#00A087",
  "#4DBBD5",
  "#F39B7F",
  "#8491B4",
  "#91D1C2",
  "#DC0000",
  "#7E6148",
  "#B09C85",
  "#6A6599",
  "#5F9E6E"
)

CT_COLORS <- setNames(
  BASE_PALETTE[
    seq_along(
      CT_LEVELS
    )
  ],
  CT_LEVELS
)

p_z <- ggplot(
  CELL_SCORE,
  aes(
    x = cell_type_plot,
    y = AF_gchromVAR_Z,
    fill = cell_type_plot
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    linetype = 2,
    colour = "grey55"
  ) +
  geom_boxplot(
    width = 0.68,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  scale_fill_manual(
    values = CT_COLORS
  ) +
  labs(
    x = NULL,
    y = "AF gchromVAR Z score"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    plot.margin = margin(
      8,
      10,
      8,
      8
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3A_AF_gchromVAR_Z_by_celltype.tiff"
  ),
  p_z,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3A_AF_gchromVAR_Z_by_celltype.pdf"
  ),
  p_z,
  width = 7.2,
  height = 4.8,
  units = "in"
)

p_seed <- ggplot(
  CT_SUMMARY,
  aes(
    x = factor(
      cell_type,
      levels =
        cell_type
    ),
    y = seed_fraction,
    fill = cell_type
  )
) +
  geom_col(
    width = 0.72,
    linewidth = 0.35,
    colour = "black"
  ) +
  scale_fill_manual(
    values =
      CT_COLORS[
        CT_SUMMARY$cell_type
      ]
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(
          100 * x,
          1
        ),
        "%"
      )
    },
    expand = expansion(
      mult = c(
        0,
        0.06
      )
    )
  ) +
  labs(
    x = NULL,
    y = "AF seed-cell fraction"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    plot.margin = margin(
      8,
      10,
      8,
      8
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3A_AF_seed_fraction_by_celltype.tiff"
  ),
  p_seed,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3A_AF_seed_fraction_by_celltype.pdf"
  ),
  p_seed,
  width = 7.2,
  height = 4.8,
  units = "in"
)

# ==============================================================================
# 18. QC / READINESS
# ==============================================================================

WEIGHT_VEC <- as.numeric(
  SummarizedExperiment::assay(
    TRAIT_IMPORT
  )[
    ,
    1
  ]
)

N_WEIGHTED_PEAKS <- sum(
  WEIGHT_VEC >
    0,
  na.rm = TRUE
)

READINESS <- data.table(
  check = c(
    "STEP11G2B_all_pass",
    "515_finemapped_blocks_present",
    "7_donors_present",
    "11986_frozen_cells_present",
    "12_author_celltypes_present",
    "212084_ATAC_peaks_present",
    "all_donor_sparse_matrix_QC_pass",
    "full_peak_by_cell_matrix_sparse",
    "all_cells_have_positive_ATAC_depth",
    "PIP_ge_0.001_retains_all_515_blocks",
    "PIP_ge_0.001_retains_ge_95pct_total_mass",
    "AF_PIP_overlaps_ATAC_peaks",
    "gchromVAR_weights_nonzero",
    "200_background_iterations_created",
    "gchromVAR_Z_finite_ge_99pct",
    "SCAVENGE_seed_cells_created",
    "seed_fraction_le_5pct",
    "scale_factor_finite_positive",
    "STEP11G3B_network_ready"
  ),
  pass = c(
    all(
      READY_G2B$pass %in%
        c(
          TRUE,
          "TRUE",
          1
        )
    ),
    nrow(
      STATUS_G2B
    ) ==
      EXPECTED_N_BLOCKS &&
      all(
        STATUS_G2B$status ==
          "PASS"
      ),
    uniqueN(
      META$donor
    ) ==
      7L,
    nrow(META) ==
      EXPECTED_N_CELLS,
    uniqueN(
      META$cell_type
    ) ==
      EXPECTED_N_CELLTYPES,
    length(PEAKS) ==
      EXPECTED_N_PEAKS,
    all(
      DONOR_QC$n_peaks ==
        EXPECTED_N_PEAKS
    ) &&
      sum(
        DONOR_QC$n_cells
      ) ==
        EXPECTED_N_CELLS &&
      all(
        DONOR_QC$peak_order_match
      ) &&
      all(
        DONOR_QC$cell_id_order_match
      ),
    inherits(
      COUNTS,
      "sparseMatrix"
    ),
    all(
      SummarizedExperiment::colData(
        SE_GC
      )$total_ATAC_counts >
        0
    ),
    N_BLOCKS_PIP ==
      EXPECTED_N_BLOCKS,
    RETAINED_MASS_FRACTION >=
      0.95,
    N_VARIANTS_OVERLAP >
      0L &&
      N_PEAKS_WEIGHTED >
        0L,
    N_WEIGHTED_PEAKS >
      0L &&
      sum(
        WEIGHT_VEC,
        na.rm = TRUE
      ) >
        0,
    nrow(BG) ==
      EXPECTED_N_PEAKS &&
      ncol(BG) ==
        BACKGROUND_ITERATIONS,
    FINITE_Z_FRACTION >=
      0.99,
    sum(SEED) >
      0L,
    mean(SEED) <=
      SEED_PERCENT +
        0.001,
    is.finite(
      SCALE_FACTOR
    ) &&
      SCALE_FACTOR >
        0,
    FINITE_Z_FRACTION >=
      0.99 &&
      sum(SEED) >
        0L
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G3A_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "dataset",
    "tissue",
    "genome_build",
    "cells",
    "author_celltypes",
    "reclustering",
    "reannotation",
    "primary_AF_finemapping",
    "PIP_threshold",
    "PIP_blocks_retained",
    "PIP_mass_fraction_retained",
    "PIP_peak_weight_import",
    "gchromVAR_background_iterations",
    "gchromVAR_seed",
    "seed_rule",
    "scale_factor_rule",
    "network_propagation_status"
  ),
  value = c(
    "GSE238242",
    "Left atrial appendage",
    "GRCh38",
    as.character(
      EXPECTED_N_CELLS
    ),
    as.character(
      EXPECTED_N_CELLTYPES
    ),
    "No",
    "No; author cell_type labels frozen",
    "STEP11G2B SuSiE-RSS L=1 uniform-prior primary PIP",
    as.character(
      PIP_MIN
    ),
    as.character(
      N_BLOCKS_PIP
    ),
    as.character(
      RETAINED_MASS_FRACTION
    ),
    IMPORT_METHOD,
    as.character(
      BACKGROUND_ITERATIONS
    ),
    "bias-corrected gchromVAR Z",
    "one-tailed P<0.05; if >5% eligible retain top 5% by Z",
    "mean Z among top 1% cells",
    "Not run; STEP11G3B after G3A QC"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G3A_method_provenance.csv"
  )
)

PACKAGE_VERSIONS <- data.table(
  package = c(
    "R",
    "Matrix",
    "SummarizedExperiment",
    "GenomicRanges",
    "chromVAR",
    "gchromVAR",
    "BSgenome.Hsapiens.UCSC.hg38"
  ),
  version = c(
    as.character(
      getRversion()
    ),
    as.character(
      packageVersion(
        "Matrix"
      )
    ),
    as.character(
      packageVersion(
        "SummarizedExperiment"
      )
    ),
    as.character(
      packageVersion(
        "GenomicRanges"
      )
    ),
    as.character(
      packageVersion(
        "chromVAR"
      )
    ),
    as.character(
      packageVersion(
        "gchromVAR"
      )
    ),
    as.character(
      packageVersion(
        "BSgenome.Hsapiens.UCSC.hg38"
      )
    )
  )
)

fwrite(
  PACKAGE_VERSIONS,
  file.path(
    QC_DIR,
    "STEP11G3A_package_versions.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G3A_sessionInfo.txt"
  )
)

# ==============================================================================
# 19. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G3A COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Peak-by-cell matrix: ",
  format(
    nrow(SE_GC),
    big.mark = ","
  ),
  " peaks x ",
  format(
    ncol(SE_GC),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "PIP >= ",
  PIP_MIN,
  ": ",
  format(
    nrow(PIP),
    big.mark = ","
  ),
  " variants across ",
  N_BLOCKS_PIP,
  " blocks\n",
  "Retained PIP mass: ",
  sprintf(
    "%.2f%%",
    100 *
      RETAINED_MASS_FRACTION
  ),
  "\n",
  sep = ""
)

cat(
  "PIP variants overlapping ATAC peaks: ",
  format(
    N_VARIANTS_OVERLAP,
    big.mark = ","
  ),
  "\n",
  "ATAC peaks receiving AF PIP weight: ",
  format(
    N_WEIGHTED_PEAKS,
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Finite gchromVAR Z: ",
  N_FINITE_Z,
  " / ",
  EXPECTED_N_CELLS,
  " (",
  sprintf(
    "%.2f%%",
    100 *
      FINITE_Z_FRACTION
  ),
  ")\n",
  sep = ""
)

cat(
  "Initial one-tailed P<0.05 cells: ",
  sum(
    SEED_RES$initial
  ),
  "\n",
  "Final SCAVENGE seed cells: ",
  sum(
    SEED
  ),
  " (",
  sprintf(
    "%.2f%%",
    100 *
      mean(
        SEED
      )
  ),
  ")\n",
  "Scale factor (top 1% mean Z): ",
  signif(
    SCALE_FACTOR,
    5
  ),
  "\n\n",
  sep = ""
)

cat(
  "Cell-type descriptive summary:\n"
)

print(
  CT_SUMMARY
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

if (!all(
  READINESS$pass
)) {

  cat(
    "\nSTEP11G3A has one or more failed readiness checks.\n",
    "Do NOT start network propagation yet.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP11G3A PASSED.\n",
    "Next: STEP11G3B TF-IDF -> LSI -> Harmony -> m-kNN -> ",
    "RWR -> TRS -> permutation inference.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD THIS FOLDER AFTER COMPLETION:\n",
  QC_DIR,
  "\n\nAnd these small result files:\n",
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_gchromVAR_celltype_summary.csv"
  ),
  "\n",
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_gchromVAR_donor_celltype_summary.csv"
  ),
  "\n",
  file.path(
    GCV_DIR,
    "STEP11G3A_AF_SCAVENGE_seed_summary.csv"
  ),
  "\n",
  CELL_SCORE_FILE,
  "\n",
  sep = ""
)
