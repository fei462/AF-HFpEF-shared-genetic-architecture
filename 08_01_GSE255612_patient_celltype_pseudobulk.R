# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP12A_GSE255612_AUDIT_AND_PATIENT_CELLTYPE_PSEUDOBULK(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP12A — GSE255612 AF snRNA-seq
#            LOCAL AUDIT + METADATA LOCK + PATIENT×CELL-TYPE PSEUDOBULK
#
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# NEXT BOARD AFTER GSE238242 S-LDSC / SCAVENGE
#
# DATASET
#   GSE255612 / SCP2489
#   Human left atrium snRNA-seq
#   Final published dataset:
#     179,697 nuclei
#     34 individuals
#     18 AF
#     16 controls
#     15 author-defined cell types
#
# PURPOSE OF STEP12A
#   1) Verify the four already-downloaded GEO files.
#   2) Read the MatrixMarket counts without constructing a Seurat object.
#   3) Align genes, barcodes and author metadata exactly.
#   4) Auto-identify and LOCK:
#        - cell barcode
#        - donor/sample
#        - AF/control status
#        - author cell type
#   5) Preserve the author annotations; NO reclustering / reannotation.
#   6) Build donor × cell-type pseudobulk raw counts.
#   7) Export QC/readiness for STEP12B AF-vs-control differential expression.
#
# IMPORTANT
#   - Statistical unit in STEP12B will be the PATIENT/DONOR, not the nucleus.
#   - This script performs NO differential expression.
#   - No candidate genes are selected from GSE255612 in this step.
#   - We first lock the independent disease-tissue dataset, then test the
#     candidates frozen upstream from genetics/S-LDSC/SCAVENGE.
#
# OUTPUT
#   D:/A/data/STEP12_GSE255612_AF_snRNA/
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))

# ==============================================================================
# 0. FROZEN EXPECTATIONS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

EXPECTED_CELLS <- 179697L
EXPECTED_GENES <- 36600L
EXPECTED_DONORS <- 34L
EXPECTED_AF <- 18L
EXPECTED_CTRL <- 16L
EXPECTED_CELLTYPES <- 15L

# For downstream pseudobulk DE:
# donor × cell-type groups with fewer than this many nuclei are retained in the
# raw pseudobulk object, but flagged as underpowered and excluded from STEP12B
# primary testing unless explicitly justified.
MIN_CELLS_PER_DONOR_CELLTYPE <- 20L

# ==============================================================================
# 1. INPUT FILES — already downloaded by the user
# ==============================================================================

MTX_FILE <- file.path(
  DATA_ROOT,
  "GSE255612_AF_snRNA_Matrix_V1.mtx.gz"
)

GENE_FILE <- file.path(
  DATA_ROOT,
  "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz"
)

BARCODE_FILE <- file.path(
  DATA_ROOT,
  "GSE255612_AF_snRNA_Processed_Expression_Matrix_barcodes_V1.tsv.gz"
)

META_FILE <- file.path(
  DATA_ROOT,
  "GSE255612_AF_snRNA_MetaData.txt.gz"
)

# ==============================================================================
# 2. OUTPUT DIRECTORIES
# ==============================================================================

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
OBJECT_DIR <- file.path(OUT_ROOT, "01_OBJECTS")
PB_DIR <- file.path(OUT_ROOT, "02_PSEUDOBULK")
TABLE_DIR <- file.path(OUT_ROOT, "03_TABLES")
FIG_DIR <- file.path(OUT_ROOT, "04_FIGURES")

for (d in c(
  OUT_ROOT,
  QC_DIR,
  OBJECT_DIR,
  PB_DIR,
  TABLE_DIR,
  FIG_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

META_LOCK_FILE <- file.path(
  OBJECT_DIR,
  "GSE255612_author_metadata_LOCKED.csv.gz"
)

COUNTS_FILE <- file.path(
  OBJECT_DIR,
  "GSE255612_counts_sparse_genes_by_cells.rds"
)

GENE_LOCK_FILE <- file.path(
  OBJECT_DIR,
  "GSE255612_gene_manifest_LOCKED.csv"
)

PB_COUNTS_FILE <- file.path(
  PB_DIR,
  "GSE255612_patient_by_celltype_pseudobulk_counts.rds"
)

PB_META_FILE <- file.path(
  PB_DIR,
  "GSE255612_patient_by_celltype_pseudobulk_metadata.csv"
)

# ==============================================================================
# 3. PACKAGES
# ==============================================================================

cran_pkgs <- c(
  "data.table",
  "Matrix",
  "ggplot2"
)

for (p in cran_pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {

    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }

  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Could not install/load package: ", p)
  }
}

library(data.table)
library(Matrix)
library(ggplot2)

# ==============================================================================
# 4. HELPERS
# ==============================================================================

norm_text <- function(x) {

  x <- as.character(x)

  x <- trimws(
    gsub(
      '"',
      "",
      x,
      fixed = TRUE
    )
  )

  x
}

safe_names <- function(x) {

  y <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )

  y <- gsub(
    "^_+|_+$",
    "",
    y
  )

  make.unique(y)
}

read_one_col <- function(path) {

  x <- fread(
    path,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )

  if (!ncol(x)) {
    stop("Empty one-column file: ", path)
  }

  norm_text(x[[1]])
}

# ------------------------------------------------------------------
# Score metadata columns by name/value structure.
# We auto-detect only if there is a clear winner.
# ------------------------------------------------------------------

column_profile <- function(DT) {

  rbindlist(
    lapply(
      names(DT),
      function(nm) {

        x <- norm_text(DT[[nm]])

        x_nonempty <- x[
          !is.na(x) &
            nzchar(x)
        ]

        values_lower <- tolower(
          unique(
            x_nonempty
          )
        )

        data.table(
          column = nm,
          n_unique =
            uniqueN(
              x_nonempty
            ),
          example =
            paste(
              head(
                unique(
                  x_nonempty
                ),
                8
              ),
              collapse = " | "
            ),
          has_AF =
            any(
              values_lower %in%
                c(
                  "af",
                  "atrial fibrillation",
                  "atrial_fibrillation",
                  "case",
                  "cases"
                )
            ) ||
            any(
              grepl(
                "atrial.?fibrillation|^af$",
                values_lower
              )
            ),
          has_CTRL =
            any(
              values_lower %in%
                c(
                  "ctrl",
                  "control",
                  "controls",
                  "normal",
                  "non-af",
                  "non_af",
                  "sinus rhythm",
                  "sr"
                )
            ) ||
            any(
              grepl(
                "control|non.?af|normal|sinus",
                values_lower
              )
            ),
          has_CM =
            any(
              grepl(
                "cardiac muscle|cardiomyocyte|cardiomyocytes|^cm$",
                values_lower
              )
            ),
          has_macrophage =
            any(
              grepl(
                "macrophage",
                values_lower
              )
            )
        )
      }
    ),
    fill = TRUE
  )
}

pick_named_column <- function(
  profiles,
  patterns,
  structural_filter = rep(TRUE, nrow(profiles))
) {

  name_score <- vapply(
    profiles$column,
    function(nm) {

      z <- tolower(nm)

      max(
        vapply(
          seq_along(patterns),
          function(i) {

            if (
              grepl(
                patterns[i],
                z
              )
            ) {
              length(patterns) -
                i +
                1
            } else {
              0
            }
          },
          numeric(1)
        )
      )
    },
    numeric(1)
  )

  idx <- which(
    structural_filter &
      name_score >
        0
  )

  if (!length(idx)) {
    return(NA_character_)
  }

  profiles$column[
    idx[
      which.max(
        name_score[idx]
      )
    ]
  ]
}

normalize_condition <- function(x) {

  y <- tolower(
    norm_text(x)
  )

  out <- rep(
    NA_character_,
    length(y)
  )

  out[
    y %in%
      c(
        "af",
        "atrial fibrillation",
        "atrial_fibrillation",
        "case",
        "cases"
      ) |
      grepl(
        "atrial.?fibrillation|^af$",
        y
      )
  ] <- "AF"

  out[
    y %in%
      c(
        "ctrl",
        "control",
        "controls",
        "normal",
        "non-af",
        "non_af",
        "sinus rhythm",
        "sr"
      ) |
      grepl(
        "control|non.?af|normal|sinus",
        y
      )
  ] <- "CTRL"

  out
}

# ==============================================================================
# 5. HARD FILE CHECK
# ==============================================================================

required_files <- c(
  MTX_FILE,
  GENE_FILE,
  BARCODE_FILE,
  META_FILE
)

if (!all(
  file.exists(
    required_files
  )
)) {

  stop(
    paste0(
      "Missing GSE255612 input file(s):\n",
      paste(
        required_files[
          !file.exists(
            required_files
          )
        ],
        collapse = "\n"
      )
    )
  )
}

FILE_QC <- data.table(
  file = basename(
    required_files
  ),
  path = required_files,
  size_GB =
    file.info(
      required_files
    )$size /
    1024^3
)

fwrite(
  FILE_QC,
  file.path(
    QC_DIR,
    "STEP12A_input_file_QC.csv"
  )
)

# ==============================================================================
# 6. READ GENES / BARCODES / METADATA
# ==============================================================================

cat("\nReading genes...\n")

GENE_RAW <- fread(
  GENE_FILE,
  header = FALSE,
  showProgress = FALSE
)

cat("Reading barcodes...\n")

BARCODES <- read_one_col(
  BARCODE_FILE
)

cat("Reading metadata...\n")

META <- fread(
  META_FILE,
  header = TRUE,
  data.table = TRUE,
  showProgress = TRUE
)

names(META) <- safe_names(
  names(META)
)

for (nm in names(META)) {

  if (
    is.character(
      META[[nm]]
    ) ||
    is.factor(
      META[[nm]]
    )
  ) {
    META[
      ,
      (nm) :=
        norm_text(
          get(nm)
        )
    ]
  }
}

# ==============================================================================
# 7. GENE MANIFEST
# ==============================================================================

if (ncol(GENE_RAW) == 1L) {

  GENE_MANIFEST <- data.table(
    gene_index = seq_len(
      nrow(
        GENE_RAW
      )
    ),
    gene_id =
      norm_text(
        GENE_RAW[[1]]
      ),
    gene_symbol =
      norm_text(
        GENE_RAW[[1]]
      )
  )

} else {

  g1 <- norm_text(
    GENE_RAW[[1]]
  )

  g2 <- norm_text(
    GENE_RAW[[2]]
  )

  # Ensembl IDs typically begin ENSG.
  if (
    mean(
      grepl(
        "^ENSG",
        g1
      )
    ) >
      mean(
        grepl(
          "^ENSG",
          g2
        )
      )
  ) {

    GENE_MANIFEST <- data.table(
      gene_index =
        seq_len(
          nrow(
            GENE_RAW
          )
        ),
      gene_id = g1,
      gene_symbol = g2
    )

  } else {

    GENE_MANIFEST <- data.table(
      gene_index =
        seq_len(
          nrow(
            GENE_RAW
          )
        ),
      gene_id = g2,
      gene_symbol = g1
    )
  }
}

GENE_MANIFEST[
  ,
  gene_id_clean :=
    sub(
      "\\.[0-9]+$",
      "",
      gene_id
    )
]

GENE_MANIFEST[
  ,
  gene_key :=
    ifelse(
      nzchar(
        gene_id
      ),
      gene_id,
      gene_symbol
    )
]

if (anyDuplicated(
  GENE_MANIFEST$gene_key
)) {

  GENE_MANIFEST[
    ,
    gene_key :=
      make.unique(
        gene_key
      )
  ]
}

fwrite(
  GENE_MANIFEST,
  GENE_LOCK_FILE
)

# ==============================================================================
# 8. METADATA COLUMN AUDIT + AUTO-DETECTION
# ==============================================================================

PROF <- column_profile(
  META
)

fwrite(
  PROF,
  file.path(
    QC_DIR,
    "STEP12A_metadata_column_profile.csv"
  )
)

# ------------------------------------------------------------------
# 8A. Barcode column:
# choose the metadata column with the highest exact overlap with GEO barcodes.
# ------------------------------------------------------------------

barcode_overlap <- vapply(
  names(META),
  function(nm) {

    x <- norm_text(
      META[[nm]]
    )

    sum(
      unique(
        x
      ) %in%
        BARCODES
    )
  },
  numeric(1)
)

BARCODE_COL <- names(
  barcode_overlap
)[
  which.max(
    barcode_overlap
  )
]

MAX_BARCODE_OVERLAP <- max(
  barcode_overlap
)

if (
  MAX_BARCODE_OVERLAP <
    0.90 *
      length(
        BARCODES
      )
) {

  # Some metadata files encode the cell ID in the first column with a
  # sample prefix.  Freeze an audit rather than guessing silently.
  fwrite(
    data.table(
      column =
        names(
          barcode_overlap
        ),
      barcode_overlap =
        as.numeric(
          barcode_overlap
        )
    )[
      order(
        -barcode_overlap
      )
    ],
    file.path(
      QC_DIR,
      "STEP12A_barcode_column_audit.csv"
    )
  )

  stop(
    paste0(
      "Could not unambiguously align metadata to GEO barcodes.\n",
      "Best metadata column = ",
      BARCODE_COL,
      "; exact overlap = ",
      MAX_BARCODE_OVERLAP,
      " / ",
      length(BARCODES),
      ".\nUpload STEP12A_barcode_column_audit.csv and ",
      "STEP12A_metadata_column_profile.csv."
    )
  )
}

# ------------------------------------------------------------------
# 8B. Cell type
# ------------------------------------------------------------------

CELLTYPE_COL <- pick_named_column(
  PROF,
  patterns = c(
    "^cell_type$",
    "^celltype$",
    "cell_type",
    "celltype",
    "annotation",
    "cell.*type",
    "cluster"
  ),
  structural_filter =
    PROF$n_unique >=
      5 &
    PROF$n_unique <=
      50
)

if (is.na(
  CELLTYPE_COL
)) {

  idx <- which(
    PROF$n_unique >=
      10 &
    PROF$n_unique <=
      25 &
    (
      PROF$has_CM |
        PROF$has_macrophage
    )
  )

  if (length(idx) == 1L) {
    CELLTYPE_COL <- PROF$column[idx]
  }
}

# ------------------------------------------------------------------
# 8C. Condition
# ------------------------------------------------------------------

CONDITION_COL <- pick_named_column(
  PROF,
  patterns = c(
    "^condition$",
    "^diagnosis$",
    "disease",
    "condition",
    "diagnosis",
    "phenotype",
    "status",
    "^af$"
  ),
  structural_filter =
    PROF$n_unique >=
      2 &
    PROF$n_unique <=
      10 &
    PROF$has_AF &
    PROF$has_CTRL
)

if (is.na(
  CONDITION_COL
)) {

  idx <- which(
    PROF$n_unique >=
      2 &
    PROF$n_unique <=
      10 &
    PROF$has_AF &
    PROF$has_CTRL
  )

  if (length(idx) == 1L) {
    CONDITION_COL <- PROF$column[idx]
  }
}

# ------------------------------------------------------------------
# 8D. Donor/sample
# Expected 34 biological samples.
# ------------------------------------------------------------------

DONOR_COL <- pick_named_column(
  PROF,
  patterns = c(
    "^donor$",
    "^sample$",
    "sample_id",
    "donor",
    "patient",
    "subject",
    "individual",
    "sample"
  ),
  structural_filter =
    PROF$n_unique >=
      30 &
    PROF$n_unique <=
      40
)

if (is.na(
  DONOR_COL
)) {

  idx <- which(
    PROF$n_unique >=
      30 &
    PROF$n_unique <=
      40
  )

  # Prefer the candidate with exactly 34 unique values.
  if (length(idx)) {

    d <- abs(
      PROF$n_unique[idx] -
        EXPECTED_DONORS
    )

    best <- idx[
      which(
        d ==
          min(d)
      )
    ]

    if (length(best) == 1L) {
      DONOR_COL <- PROF$column[best]
    }
  }
}

DETECTION <- data.table(
  role = c(
    "barcode",
    "donor",
    "condition",
    "cell_type"
  ),
  column = c(
    BARCODE_COL,
    DONOR_COL,
    CONDITION_COL,
    CELLTYPE_COL
  )
)

fwrite(
  DETECTION,
  file.path(
    QC_DIR,
    "STEP12A_detected_metadata_columns.csv"
  )
)

if (
  any(
    is.na(
      DETECTION$column
    )
  )
) {

  print(
    DETECTION
  )

  stop(
    paste0(
      "One or more author metadata fields could not be identified ",
      "unambiguously.\nUpload:\n",
      file.path(
        QC_DIR,
        "STEP12A_metadata_column_profile.csv"
      ),
      "\n",
      file.path(
        QC_DIR,
        "STEP12A_detected_metadata_columns.csv"
      )
    )
  )
}

# ==============================================================================
# 9. ALIGN AUTHOR METADATA TO BARCODE FILE
# ==============================================================================

META[
  ,
  barcode :=
    norm_text(
      get(
        BARCODE_COL
      )
    )
]

META[
  ,
  donor :=
    norm_text(
      get(
        DONOR_COL
      )
    )
]

META[
  ,
  condition_original :=
    norm_text(
      get(
        CONDITION_COL
      )
    )
]

META[
  ,
  condition :=
    normalize_condition(
      condition_original
    )
]

META[
  ,
  cell_type :=
    norm_text(
      get(
        CELLTYPE_COL
      )
    )
]

if (
  anyDuplicated(
    META$barcode
  )
) {
  stop("Author metadata contains duplicated cell barcodes.")
}

idx_meta <- match(
  BARCODES,
  META$barcode
)

if (anyNA(
  idx_meta
)) {
  stop(
    "Some GEO barcodes cannot be matched to author metadata."
  )
}

META_LOCK <- META[
  idx_meta
]

if (!identical(
  META_LOCK$barcode,
  BARCODES
)) {
  stop("Metadata/barcode order locking failed.")
}

if (anyNA(
  META_LOCK$condition
)) {

  bad <- unique(
    META_LOCK[
      is.na(
        condition
      ),
      condition_original
    ]
  )

  stop(
    paste0(
      "Condition normalization left unknown labels: ",
      paste(
        bad,
        collapse = ", "
      )
    )
  )
}

# A donor must map to exactly one AF/CTRL status.
DONOR_COND <- unique(
  META_LOCK[
    ,
    .(
      donor,
      condition
    )
  ]
)

donor_condition_n <- DONOR_COND[
  ,
  .N,
  by = donor
]

if (any(
  donor_condition_n$N !=
    1L
)) {
  stop(
    "At least one donor maps to multiple AF/control labels."
  )
}

fwrite(
  META_LOCK,
  META_LOCK_FILE,
  compress = "gzip"
)

# ==============================================================================
# 10. LOCK DATASET COUNTS / SAMPLE STRUCTURE BEFORE MATRIX READ
# ==============================================================================

N_CELLS_META <- nrow(
  META_LOCK
)

N_DONORS <- uniqueN(
  META_LOCK$donor
)

N_CT <- uniqueN(
  META_LOCK$cell_type
)

DONOR_STATUS <- unique(
  META_LOCK[
    ,
    .(
      donor,
      condition
    )
  ]
)

N_AF <- sum(
  DONOR_STATUS$condition ==
    "AF"
)

N_CTRL <- sum(
  DONOR_STATUS$condition ==
    "CTRL"
)

STRUCTURE_QC <- data.table(
  metric = c(
    "n_barcodes",
    "n_metadata_rows",
    "n_genes_manifest",
    "n_donors",
    "n_AF_donors",
    "n_CTRL_donors",
    "n_author_celltypes",
    "barcode_metadata_exact_order_match"
  ),
  value = c(
    length(BARCODES),
    N_CELLS_META,
    nrow(GENE_MANIFEST),
    N_DONORS,
    N_AF,
    N_CTRL,
    N_CT,
    identical(
      META_LOCK$barcode,
      BARCODES
    )
  )
)

fwrite(
  STRUCTURE_QC,
  file.path(
    QC_DIR,
    "STEP12A_dataset_structure_QC.csv"
  )
)

# Freeze author cell type names before any downstream analysis.
CELLTYPE_COUNTS <- META_LOCK[
  ,
  .(
    n_cells = .N,
    n_donors =
      uniqueN(
        donor
      ),
    n_AF_cells =
      sum(
        condition ==
          "AF"
      ),
    n_CTRL_cells =
      sum(
        condition ==
          "CTRL"
      )
  ),
  by = cell_type
][
  order(
    -n_cells
  )
]

fwrite(
  CELLTYPE_COUNTS,
  file.path(
    TABLE_DIR,
    "GSE255612_author_celltype_counts.csv"
  )
)

# ==============================================================================
# 11. READ MATRIXMARKET SPARSE COUNTS
# ==============================================================================

if (file.exists(
  COUNTS_FILE
)) {

  cat(
    "\nReusing locked sparse count checkpoint...\n"
  )

  COUNTS <- readRDS(
    COUNTS_FILE
  )

} else {

  cat(
    "\nReading 956-MB compressed MatrixMarket counts.\n",
    "This is the only large read in STEP12A...\n",
    sep = ""
  )

  con <- gzfile(
    MTX_FILE,
    open = "rt"
  )

  COUNTS <- Matrix::readMM(
    con
  )

  close(
    con
  )

  COUNTS <- methods::as(
    COUNTS,
    "dgCMatrix"
  )

  # Detect and correct orientation only from exact dimensions.
  if (
    nrow(COUNTS) ==
      nrow(GENE_MANIFEST) &&
    ncol(COUNTS) ==
      length(BARCODES)
  ) {

    MATRIX_ORIENTATION <-
      "genes_by_cells"

  } else if (
    nrow(COUNTS) ==
      length(BARCODES) &&
    ncol(COUNTS) ==
      nrow(GENE_MANIFEST)
  ) {

    COUNTS <- Matrix::t(
      COUNTS
    )

    COUNTS <- methods::as(
      COUNTS,
      "dgCMatrix"
    )

    MATRIX_ORIENTATION <-
      "cells_by_genes_transposed_to_genes_by_cells"

  } else {

    stop(
      paste0(
        "Matrix dimensions do not match genes/barcodes.\n",
        "Matrix = ",
        nrow(COUNTS),
        " x ",
        ncol(COUNTS),
        "\nGenes = ",
        nrow(GENE_MANIFEST),
        "\nBarcodes = ",
        length(BARCODES)
      )
    )
  }

  rownames(
    COUNTS
  ) <- GENE_MANIFEST$gene_key

  colnames(
    COUNTS
  ) <- BARCODES

  saveRDS(
    COUNTS,
    COUNTS_FILE,
    compress = FALSE
  )
}

if (
  nrow(COUNTS) !=
    nrow(GENE_MANIFEST) ||
  ncol(COUNTS) !=
    length(BARCODES)
) {
  stop("Locked sparse counts dimension mismatch.")
}

if (!identical(
  colnames(COUNTS),
  META_LOCK$barcode
)) {
  stop("Sparse count matrix and locked metadata cell order differ.")
}

# ==============================================================================
# 12. CELL-LEVEL BASIC COUNT QC
# ==============================================================================

CELL_UMI <- Matrix::colSums(
  COUNTS
)

CELL_NGENE <- Matrix::colSums(
  COUNTS >
    0
)

CELL_QC <- data.table(
  metric = c(
    "n_cells",
    "n_genes",
    "matrix_nnz",
    "min_UMI",
    "median_UMI",
    "max_UMI",
    "min_detected_genes",
    "median_detected_genes",
    "max_detected_genes"
  ),
  value = c(
    ncol(COUNTS),
    nrow(COUNTS),
    length(COUNTS@x),
    min(CELL_UMI),
    median(CELL_UMI),
    max(CELL_UMI),
    min(CELL_NGENE),
    median(CELL_NGENE),
    max(CELL_NGENE)
  )
)

fwrite(
  CELL_QC,
  file.path(
    QC_DIR,
    "STEP12A_sparse_count_QC.csv"
  )
)

# Do NOT re-filter author nuclei here.
# The publication already performed QC and this board is independent validation.

# ==============================================================================
# 13. BUILD DONOR × CELL-TYPE PSEUDOBULK
# ==============================================================================

cat(
  "\nBuilding patient × cell-type pseudobulk counts...\n"
)

PB_CELL_META <- META_LOCK[
  ,
  .(
    barcode,
    donor,
    condition,
    cell_type
  )
]

PB_CELL_META[
  ,
  pb_id :=
    paste(
      donor,
      cell_type,
      sep = "__"
    )
]

PB_LEVELS <- unique(
  PB_CELL_META$pb_id
)

PB_FACTOR <- factor(
  PB_CELL_META$pb_id,
  levels = PB_LEVELS
)

# Sparse cell x pseudobulk membership.
DESIGN_PB <- Matrix::sparse.model.matrix(
  ~ 0 + PB_FACTOR
)

colnames(
  DESIGN_PB
) <- sub(
  "^PB_FACTOR",
  "",
  colnames(
    DESIGN_PB
  )
)

# genes x cells  %*%  cells x patient-celltype
PB_COUNTS <- COUNTS %*%
  DESIGN_PB

PB_COUNTS <- methods::as(
  PB_COUNTS,
  "dgCMatrix"
)

PB_SAMPLE_META <- PB_CELL_META[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    pb_id,
    donor,
    condition,
    cell_type
  )
]

idx_pb <- match(
  colnames(
    PB_COUNTS
  ),
  PB_SAMPLE_META$pb_id
)

if (anyNA(
  idx_pb
)) {
  stop(
    "Pseudobulk matrix/sample metadata alignment failed."
  )
}

PB_SAMPLE_META <- PB_SAMPLE_META[
  idx_pb
]

PB_SAMPLE_META[
  ,
  library_size :=
    as.numeric(
      Matrix::colSums(
        PB_COUNTS
      )
    )
]

PB_SAMPLE_META[
  ,
  genes_detected :=
    as.integer(
      Matrix::colSums(
        PB_COUNTS >
          0
      )
    )
]

PB_SAMPLE_META[
  ,
  eligible_primary_DE :=
    n_cells >=
      MIN_CELLS_PER_DONOR_CELLTYPE &
    library_size >
      0
]

saveRDS(
  PB_COUNTS,
  PB_COUNTS_FILE,
  compress = TRUE
)

fwrite(
  PB_SAMPLE_META,
  PB_META_FILE
)

# ==============================================================================
# 14. PSEUDOBULK READINESS PER CELL TYPE
# ==============================================================================

CT_PB_READY <- PB_SAMPLE_META[
  ,
  .(
    n_donors_all =
      uniqueN(
        donor
      ),
    n_AF_donors_all =
      uniqueN(
        donor[
          condition ==
            "AF"
        ]
      ),
    n_CTRL_donors_all =
      uniqueN(
        donor[
          condition ==
            "CTRL"
        ]
      ),
    n_eligible_donors =
      uniqueN(
        donor[
          eligible_primary_DE
        ]
      ),
    n_eligible_AF =
      uniqueN(
        donor[
          eligible_primary_DE &
            condition ==
              "AF"
        ]
      ),
    n_eligible_CTRL =
      uniqueN(
        donor[
          eligible_primary_DE &
            condition ==
              "CTRL"
        ]
      ),
    median_cells_per_donor =
      as.numeric(
        median(
          n_cells
        )
      ),
    min_cells_per_donor =
      min(
        n_cells
      )
  ),
  by = cell_type
]

# Primary DE readiness criterion:
# at least 5 AF + 5 CTRL biological donors with >=20 nuclei each.
CT_PB_READY[
  ,
  ready_for_STEP12B :=
    n_eligible_AF >=
      5L &
    n_eligible_CTRL >=
      5L
]

setorder(
  CT_PB_READY,
  -ready_for_STEP12B,
  -n_eligible_donors
)

fwrite(
  CT_PB_READY,
  file.path(
    QC_DIR,
    "STEP12A_celltype_pseudobulk_readiness.csv"
  )
)

# Identify author cardiomyocyte label for the PRIMARY validation cell type.
CM_LABELS <- CT_PB_READY[
  grepl(
    "cardiac muscle|cardiomyocyte|cardiomyocytes|^cm$",
    cell_type,
    ignore.case = TRUE
  ),
  cell_type
]

MAC_LABELS <- CT_PB_READY[
  grepl(
    "macrophage",
    cell_type,
    ignore.case = TRUE
  ),
  cell_type
]

TARGET_CT <- data.table(
  target = c(
    "Primary_genetics_supported_celltype",
    "Published_AF_DE_context_celltype"
  ),
  expected_biology = c(
    "Cardiomyocyte",
    "Macrophage"
  ),
  detected_author_label = c(
    if (
      length(
        CM_LABELS
      ) ==
        1L
    ) {
      CM_LABELS
    } else {
      paste(
        CM_LABELS,
        collapse = " | "
      )
    },
    if (
      length(
        MAC_LABELS
      ) ==
        1L
    ) {
      MAC_LABELS
    } else {
      paste(
        MAC_LABELS,
        collapse = " | "
      )
    }
  )
)

fwrite(
  TARGET_CT,
  file.path(
    QC_DIR,
    "STEP12A_target_celltype_labels.csv"
  )
)

# ==============================================================================
# 15. DESCRIPTIVE FIGURE — nuclei per donor/cell type
# ==============================================================================

PLOT_CT <- CELLTYPE_COUNTS[
  order(
    n_cells
  )
]

PLOT_CT[
  ,
  cell_type_plot :=
    factor(
      cell_type,
      levels =
        cell_type
    )
]

p <- ggplot(
  PLOT_CT,
  aes(
    x = n_cells,
    y = cell_type_plot
  )
) +
  geom_col(
    width = 0.72
  ) +
  labs(
    x = "Number of nuclei",
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.55
    )
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP12A_GSE255612_author_celltype_counts.tiff"
  ),
  p,
  width = 6.2,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================================================================
# 16. READINESS
# ==============================================================================

CM_READY <- (
  length(
    CM_LABELS
  ) ==
    1L &&
  CT_PB_READY[
    cell_type ==
      CM_LABELS,
    ready_for_STEP12B
  ][1] %in%
    TRUE
)

READINESS <- data.table(
  check = c(
    "four_GEO_files_present",
    "179697_barcodes_present",
    "179697_metadata_rows_present",
    "36600_genes_present",
    "metadata_barcode_alignment_exact",
    "34_donors_present",
    "18_AF_donors_present",
    "16_CTRL_donors_present",
    "15_author_celltypes_present",
    "sparse_matrix_dimensions_correct",
    "pseudobulk_matrix_created",
    "pseudobulk_metadata_created",
    "cardiomyocyte_author_label_unique",
    "cardiomyocyte_ready_for_patient_level_DE",
    "STEP12A_complete"
  ),
  pass = c(
    all(
      file.exists(
        required_files
      )
    ),
    length(
      BARCODES
    ) ==
      EXPECTED_CELLS,
    nrow(
      META_LOCK
    ) ==
      EXPECTED_CELLS,
    nrow(
      GENE_MANIFEST
    ) ==
      EXPECTED_GENES,
    identical(
      META_LOCK$barcode,
      BARCODES
    ),
    N_DONORS ==
      EXPECTED_DONORS,
    N_AF ==
      EXPECTED_AF,
    N_CTRL ==
      EXPECTED_CTRL,
    N_CT ==
      EXPECTED_CELLTYPES,
    nrow(
      COUNTS
    ) ==
      EXPECTED_GENES &&
      ncol(
        COUNTS
      ) ==
        EXPECTED_CELLS,
    nrow(
      PB_COUNTS
    ) ==
      EXPECTED_GENES &&
      ncol(
        PB_COUNTS
      ) ==
        nrow(
          PB_SAMPLE_META
        ),
    file.exists(
      PB_META_FILE
    ),
    length(
      CM_LABELS
    ) ==
      1L,
    CM_READY,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP12A_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "dataset",
    "tissue",
    "assay",
    "published_final_nuclei",
    "published_final_AF",
    "published_final_control",
    "published_cell_types",
    "annotation_policy",
    "statistical_unit_for_next_step",
    "pseudobulk_definition",
    "minimum_nuclei_per_donor_celltype_primary_DE",
    "primary_validation_celltype",
    "STEP12A_differential_expression"
  ),
  value = c(
    "GSE255612 / SCP2489",
    "Human left atrium",
    "snRNA-seq",
    as.character(
      EXPECTED_CELLS
    ),
    as.character(
      EXPECTED_AF
    ),
    as.character(
      EXPECTED_CTRL
    ),
    as.character(
      EXPECTED_CELLTYPES
    ),
    "Use author metadata; no reclustering or reannotation",
    "Biological donor/patient",
    "Raw UMI counts summed within donor x author cell type",
    as.character(
      MIN_CELLS_PER_DONOR_CELLTYPE
    ),
    "Cardiomyocyte, motivated upstream by AF S-LDSC/SCAVENGE",
    "Not performed; STEP12B after QC"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP12A_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP12A_sessionInfo.txt"
  )
)

# ==============================================================================
# 17. FINAL SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP12A GSE255612 COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

print(
  STRUCTURE_QC
)

cat(
  "\nDetected metadata fields:\n"
)
print(
  DETECTION
)

cat(
  "\nAuthor cell-type pseudobulk readiness:\n"
)
print(
  CT_PB_READY
)

cat(
  "\nPrimary/control target labels:\n"
)
print(
  TARGET_CT
)

cat(
  "\nReadiness:\n"
)
print(
  READINESS
)

if (all(
  READINESS$pass
)) {

  cat(
    "\nSTEP12A PASSED.\n",
    "Next: STEP12B patient-level AF vs CTRL pseudobulk DE, ",
    "with cardiomyocyte as the prespecified primary cell type.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP12A has failed readiness checks.\n",
    "Do not run differential expression yet.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  PB_META_FILE,
  "\n",
  sep = ""
)
