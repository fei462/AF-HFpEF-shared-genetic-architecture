# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G3B_V3_OFFICIAL_EQUIVALENT_SCAVENGE_1000PERM(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G3B V2 — FULL SCAVENGE NETWORK PROPAGATION + 1000 PERMUTATIONS
# Project: AF–HFpEF–BMI–OSA shared genetics
# Dataset: GSE238242 left atrial appendage snATAC
#
# PREREQUISITE:
#   STEP11G3A must be fully PASS.
#
# FROZEN PRIMARY SETTINGS (SCAVENGE workflow):
#   peak-by-cell matrix: all 11,986 frozen cells
#   TF-IDF: binary=TRUE, TF=TRUE, log_TF=TRUE
#   LSI dimensions: 30
#   Harmony correction: donor ONLY
#       (rhythm is NOT corrected because it may contain biological signal)
#   mutual kNN: k = 30
#   random walk restart gamma = 0.05
#   RWR stationary cutoff = 1e-5
#   TRS: cap NP at 95th percentile -> min-max -> x scale factor
#   permutation: 1,000 degree-matched seed permutations
#   cell significance: empirical P <= 0.05
#
# IMPORTANT:
#   - NO reclustering
#   - NO reannotation
#   - NO pre-filtering to CM
#   - NO correction for rhythm/AF status
#   - all author-defined cell types are retained
#
# ROBUSTNESS:
#   - V2 fixes degree-stratum list indexing: ALL_BY_DEGREE[[deg_name]]
#   - checkpoint raw LSI
#   - checkpoint Harmony LSI
#   - checkpoint mutual-kNN graph
#   - checkpoint real RWR/TRS
#   - permutation is run in batches and each batch is checkpointed
#   - rerunning the script resumes completed work
#
# PERMUTATION IMPLEMENTATION:
#   SCAVENGE's published null samples seed cells matched on graph degree,
#   then reruns the same random walk. This script preserves that exact null
#   logic but propagates many permutations simultaneously in matrix batches,
#   making it Windows-safe and restartable.
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

LSI_DIMS <- 30L
K_MKNN <- 30L

RW_GAMMA <- 0.05
RW_STATIONARY_CUTOFF <- 1e-5
RW_MAX_ITER <- 10000L

TRS_CAP_QUANTILE <- 0.95

PERMUTATION_TIMES <- 1000L
PERM_BATCH_SIZE <- 100L
TRUE_CELL_ALPHA <- 0.05

HARMONY_VARIABLE <- "donor"
HARMONY_NCORES <- 1L

# UMAP is visualization only.
RUN_UMAP <- TRUE
UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.30

# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

G3A_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "04_GCHROMVAR"
)

G3A_QC <- file.path(
  G3A_ROOT,
  "00_QC"
)

G3A_READINESS <- file.path(
  G3A_QC,
  "STEP11G3A_readiness.csv"
)

RSE_FILE <- file.path(
  G3A_ROOT,
  "01_SCATAC_MATRIX",
  "GSE238242_SCAVENGE_peak_by_cell_GRCh38.rds"
)

CELL_SCORE_FILE <- file.path(
  G3A_ROOT,
  "03_GCHROMVAR",
  "STEP11G3A_AF_cell_gchromVAR_seed_scores.csv.gz"
)

SEED_SUMMARY_FILE <- file.path(
  G3A_ROOT,
  "03_GCHROMVAR",
  "STEP11G3A_AF_SCAVENGE_seed_summary.csv"
)

OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "05_SCAVENGE_NETWORK_V3_OFFICIAL_EQUIV"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

EMBED_DIR <- file.path(
  OUT_ROOT,
  "01_EMBEDDING"
)

NETWORK_DIR <- file.path(
  OUT_ROOT,
  "02_NETWORK"
)

RWR_DIR <- file.path(
  OUT_ROOT,
  "03_RWR_TRS"
)

PERM_DIR <- file.path(
  OUT_ROOT,
  "04_PERMUTATION"
)

FINAL_DIR <- file.path(
  OUT_ROOT,
  "05_FINAL"
)

FIG_DIR <- file.path(
  OUT_ROOT,
  "06_FIGURES"
)

for (d in c(
  OUT_ROOT,
  QC_DIR,
  EMBED_DIR,
  NETWORK_DIR,
  RWR_DIR,
  PERM_DIR,
  FINAL_DIR,
  FIG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LSI_FILE <- file.path(
  EMBED_DIR,
  "STEP11G3B_LSI30_raw.rds"
)

HARMONY_FILE <- file.path(
  EMBED_DIR,
  "STEP11G3B_LSI30_Harmony_donor.rds"
)

HARMONY_QC_FILE <- file.path(
  QC_DIR,
  "STEP11G3B_Harmony_donor_QC.csv"
)

UMAP_FILE <- file.path(
  EMBED_DIR,
  "STEP11G3B_Harmony_UMAP.rds"
)

GRAPH_FILE <- file.path(
  NETWORK_DIR,
  "STEP11G3B_mutualKNN_k30.rds"
)

GRAPH_QC_FILE <- file.path(
  QC_DIR,
  "STEP11G3B_graph_QC.csv"
)

REAL_RWR_FILE <- file.path(
  RWR_DIR,
  "STEP11G3B_real_RWR_TRS.rds"
)

CELL_RESULT_FILE <- file.path(
  FINAL_DIR,
  "STEP11G3B_AF_SCAVENGE_cell_results.csv.gz"
)

CELLTYPE_RESULT_FILE <- file.path(
  FINAL_DIR,
  "STEP11G3B_AF_SCAVENGE_celltype_summary.csv"
)

CELLTYPE_ENRICH_FILE <- file.path(
  FINAL_DIR,
  "STEP11G3B_AF_SCAVENGE_celltype_enrichment_Fisher.csv"
)

DONOR_CELLTYPE_FILE <- file.path(
  FINAL_DIR,
  "STEP11G3B_AF_SCAVENGE_donor_celltype_summary.csv"
)

RHYTHM_CELLTYPE_FILE <- file.path(
  FINAL_DIR,
  "STEP11G3B_AF_SCAVENGE_rhythm_celltype_DESCRIPTIVE.csv"
)

READINESS_FILE <- file.path(
  QC_DIR,
  "STEP11G3B_readiness.csv"
)

# ==============================================================================
# 2. PACKAGES
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
      "Failed to install/load CRAN package: ",
      pkg
    )
  }
}

for (p in c(
  "data.table",
  "Matrix",
  "irlba",
  "RANN",
  "igraph",
  "harmony",
  "ggplot2"
)) {
  install_cran_if_missing(p)
}

if (RUN_UMAP) {
  install_cran_if_missing("uwot")
}

if (!requireNamespace(
  "SummarizedExperiment",
  quietly = TRUE
)) {

  if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
  )) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }

  BiocManager::install(
    "SummarizedExperiment",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace(
  "SummarizedExperiment",
  quietly = TRUE
)) {
  stop(
    "SummarizedExperiment is required."
  )
}

library(data.table)
library(Matrix)
library(ggplot2)

# Optional official package:
# If already installed, we use its TF-IDF / LSI / mKNN / TRS utilities.
# We do NOT force a GitHub installation because this workstation has shown
# certificate/network restrictions; exact-compatible local fallbacks are below.
HAS_SCAVENGE <- requireNamespace(
  "SCAVENGE",
  quietly = TRUE
)

# ==============================================================================
# 3. SCAVENGE OFFICIAL-SOURCE-EQUIVALENT CORE FUNCTIONS
# ==============================================================================
#
# These local functions reproduce the current public SCAVENGE source logic for:
#   tfidf()
#   do_lsi()
#   getmutualknn()
#   capOutlierQuantile()
#   max_min_scale()
#
# This avoids dependence on GitHub installation while preserving the published
# package behavior.  In particular:
#   - TF-IDF uses log1p(), NOT log()
#   - getmutualknn() keeps reciprocal neighbours and, when none exist,
#     connects the cell to its nearest neighbour.
# ==============================================================================

sc_tfidf_source_equiv <- function(
  bmat,
  mat_binary = TRUE,
  TF = TRUE,
  log_TF = TRUE,
  scale_factor = 100000
) {

  bmat <- methods::as(
    bmat,
    "dgCMatrix"
  )

  if (mat_binary) {
    bmat@x[
      bmat@x >= 1
    ] <- 1
    message("[info] binarize matrix")
  }

  if (TF) {

    cs <- Matrix::colSums(
      bmat
    )

    if (any(
      !is.finite(cs) |
        cs <= 0
    )) {
      stop(
        "TF-IDF input contains zero/invalid cell depth."
      )
    }

    tf <- Matrix::t(
      Matrix::t(
        bmat
      ) /
        cs
    )

  } else {

    tf <- bmat
  }

  message("[info] calculate tf")

  if (log_TF) {

    if (TF) {

      tf@x <- log1p(
        tf@x *
          scale_factor
      )

    } else {

      tf@x <- log1p(
        tf@x
      )
    }
  }

  message("[info] calculate idf")

  idf <- log(
    1 +
      ncol(bmat) /
      Matrix::rowSums(
        bmat
      )
  )

  fast_tfidf <- function(
    tf,
    idf
  ) {

    tf <- Matrix::t(
      tf
    )

    tf@x <- tf@x *
      rep.int(
        idf,
        diff(
          tf@p
        )
      )

    tf <- Matrix::t(
      tf
    )

    tf
  }

  message("[info] fast log tf-idf")

  tf_idf_counts <- fast_tfidf(
    tf,
    idf
  )

  rownames(
    tf_idf_counts
  ) <- rownames(
    bmat
  )

  colnames(
    tf_idf_counts
  ) <- colnames(
    bmat
  )

  methods::as(
    tf_idf_counts,
    "dgCMatrix"
  )
}

sc_do_lsi_source_equiv <- function(
  mat,
  dims = 30
) {

  message(
    "SVD analysis of TF-IDF matrix"
  )

  pca_results <- irlba::irlba(
    Matrix::t(
      mat
    ),
    nv = dims,
    fastpath = FALSE
  )

  PCA_result <- pca_results$u %*%
    diag(
      pca_results$d
    )

  rownames(
    PCA_result
  ) <- colnames(
    mat
  )

  colnames(
    PCA_result
  ) <- paste0(
    "LSI_",
    seq_len(
      dims
    )
  )

  PCA_result
}

sc_getmutualknn_source_equiv <- function(
  lsimat,
  num_k = 30
) {

  stopifnot(
    "num_k must be numeric" =
      is.numeric(
        num_k
      )
  )

  if (
    nrow(lsimat) <=
      num_k
  ) {
    stop(
      "Number of cells must exceed k."
    )
  }

  mymat <- as.matrix(
    lsimat
  )

  K <- num_k

  message(
    "[info] fast knn"
  )

  # Current SCAVENGE source uses RANN::nn2(mymat, k=K).
  # The first index is the cell itself.
  knn_info <- RANN::nn2(
    mymat,
    k = K
  )

  knn <- knn_info$nn.idx

  edge_list <- vector(
    "list",
    nrow(
      knn
    )
  )

  for (
    i in seq_len(
      nrow(
        knn
      )
    )
  ) {

    candidates <- knn[
      i,
      -1,
      drop = TRUE
    ]

    reciprocal <- vapply(
      candidates,
      function(j) {
        any(
          knn[
            j,
            ,
            drop = TRUE
          ] ==
            i
        )
      },
      logical(1)
    )

    if (sum(
      reciprocal
    ) >
      0L) {

      neighbours <- candidates[
        reciprocal
      ]

    } else {

      # Exact SCAVENGE fallback semantics:
      # if there are no reciprocal neighbours, keep the nearest neighbour.
      neighbours <- knn[
        i,
        2
      ]
    }

    edge_list[[
      i
    ]] <- cbind(
      rep.int(
        i,
        length(
          neighbours
        )
      ),
      neighbours
    )
  }

  el <- do.call(
    rbind,
    edge_list
  )

  n <- nrow(
    knn
  )

  adj <- Matrix::sparseMatrix(
    i = el[
      ,
      1
    ],
    j = el[
      ,
      2
    ],
    x = 1,
    dims = c(
      n,
      n
    )
  )

  # SCAVENGE source:
  # mutualknn <- 1*((adj + t(adj)) > 0)
  mutualknn <- (
    (
      adj +
        Matrix::t(
          adj
        )
    ) >
      0
  ) *
    1

  Matrix::diag(
    mutualknn
  ) <- 0

  mutualknn <- methods::as(
    mutualknn,
    "dgCMatrix"
  )

  colnames(
    mutualknn
  ) <- rownames(
    mutualknn
  ) <- rownames(
    lsimat
  )

  mutualknn
}

sc_cap_outlier_source_equiv <- function(
  x,
  q_ceiling = 0.95
) {

  q <- stats::quantile(
    x,
    q_ceiling
  )

  x[
    x >
      q
  ] <- q

  x
}

sc_max_min_source_equiv <- function(
  x
) {

  (
    x -
      min(
        x
      )
  ) /
    (
      max(
        x
      ) -
        min(
          x
        )
    )
}

# ==============================================================================
# 4. RWR ENGINE
# ==============================================================================

# Build the same column-normalized transition matrix used by SCAVENGE.
make_transition <- function(
  adj
) {

  adj <- methods::as(
    adj,
    "dgCMatrix"
  )

  cs <- Matrix::colSums(
    adj
  )

  if (any(
    !is.finite(cs) |
      cs <= 0
  )) {
    stop(
      "Transition graph has zero/invalid degree columns."
    )
  }

  Matrix::t(
    Matrix::t(
      adj
    ) /
      cs
  )
}

# Single/multi-seed RWR.
# P0 is cells x walks; every column must sum to 1.
rwr_matrix <- function(
  W,
  P0,
  gamma = 0.05,
  stationary_cutoff = 1e-5,
  max_iter = 10000L,
  verbose = FALSE
) {

  P <- as.matrix(
    P0
  )

  if (
    nrow(P) !=
      nrow(W)
  ) {
    stop(
      "P0 / graph dimension mismatch."
    )
  }

  cs <- colSums(
    P
  )

  if (any(
    !is.finite(cs) |
      cs <= 0
  )) {
    stop(
      "Every P0 column must contain positive seed mass."
    )
  }

  P <- sweep(
    P,
    2,
    cs,
    "/"
  )

  for (iter in seq_len(
    max_iter
  )) {

    P_new <- (
      1 -
        gamma
    ) *
      (
        W %*%
          P
      ) +
      gamma *
      P0

    delta_each <- colSums(
      abs(
        as.matrix(
          P_new
        ) -
          P
      )
    )

    delta <- max(
      delta_each
    )

    P <- as.matrix(
      P_new
    )

    if (
      is.finite(delta) &&
      delta <=
        stationary_cutoff
    ) {

      if (verbose) {
        message(
          "Stationary step: ",
          iter
        )
        message(
          "Stationary Delta: ",
          signif(
            delta,
            6
          )
        )
      }

      return(
        list(
          score = P,
          iterations = iter,
          delta = delta
        )
      )
    }
  }

  stop(
    "RWR failed to reach stationary cutoff within ",
    max_iter,
    " iterations."
  )
}

# ==============================================================================
# 5. HARD INPUT / STEP11G3A CHECK
# ==============================================================================

required_inputs <- c(
  G3A_readiness =
    G3A_READINESS,
  peak_by_cell_RSE =
    RSE_FILE,
  cell_gchromVAR =
    CELL_SCORE_FILE,
  seed_summary =
    SEED_SUMMARY_FILE
)

missing <- required_inputs[
  !file.exists(
    required_inputs
  )
]

if (length(missing)) {
  stop(
    paste0(
      "Missing STEP11G3B input(s):\n",
      paste(
        names(missing),
        missing,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

READY_G3A <- fread(
  G3A_READINESS
)

if (
  nrow(
    READY_G3A
  ) != 19L ||
  !all(
    READY_G3A$pass %in%
      c(
        TRUE,
        "TRUE",
        1
      )
  )
) {

  print(
    READY_G3A
  )

  stop(
    "STEP11G3A is not fully PASS."
  )
}

SE <- readRDS(
  RSE_FILE
)

COUNTS <- SummarizedExperiment::assay(
  SE,
  "counts"
)

if (!inherits(
  COUNTS,
  "sparseMatrix"
)) {
  stop(
    "Peak-by-cell matrix is not sparse."
  )
}

if (
  nrow(COUNTS) !=
    EXPECTED_N_PEAKS ||
  ncol(COUNTS) !=
    EXPECTED_N_CELLS
) {
  stop(
    "Peak-by-cell matrix dimensions are not 212084 x 11986."
  )
}

CELL <- fread(
  CELL_SCORE_FILE
)

required_cell_cols <- c(
  "cell_id",
  "cell_type",
  "donor",
  "rhythm",
  "AF_gchromVAR_Z",
  "AF_SCAVENGE_seed"
)

if (!all(
  required_cell_cols %in%
    names(CELL)
)) {
  stop(
    "STEP11G3A cell-score file lacks required fields."
  )
}

if (
  nrow(CELL) !=
    EXPECTED_N_CELLS
) {
  stop(
    "Expected 11,986 cell-score rows."
  )
}

# G3A output is ranked by gchromVAR Z, whereas the ATAC matrix keeps the
# original donor/cell order. Re-align exactly to the matrix.
idx <- match(
  colnames(
    COUNTS
  ),
  CELL$cell_id
)

if (anyNA(idx)) {
  stop(
    "Could not align G3A cell scores to peak-by-cell matrix."
  )
}

CELL <- CELL[
  idx
]

if (!identical(
  CELL$cell_id,
  colnames(
    COUNTS
  )
)) {
  stop(
    "Cell-order alignment failed."
  )
}

if (
  uniqueN(
    CELL$cell_type
  ) !=
    EXPECTED_N_CELLTYPES
) {
  stop(
    "Expected 12 author cell types."
  )
}

if (
  uniqueN(
    CELL$donor
  ) !=
    7L
) {
  stop(
    "Expected 7 donors."
  )
}

if (
  sum(
    CELL$AF_SCAVENGE_seed
  ) !=
    599L
) {

  stop(
    "Expected 599 frozen G3A seed cells; observed ",
    sum(
      CELL$AF_SCAVENGE_seed
    )
  )
}

SCALE_SUM <- fread(
  SEED_SUMMARY_FILE
)

SCALE_FACTOR <- as.numeric(
  SCALE_SUM[
    metric ==
      "scale_factor",
    value
  ]
)

if (
  length(
    SCALE_FACTOR
  ) != 1L ||
  !is.finite(
    SCALE_FACTOR
  ) ||
  SCALE_FACTOR <= 0
) {
  stop(
    "Could not recover valid STEP11G3A scale factor."
  )
}

cat(
  "\nFrozen STEP11G3A input:\n",
  "  cells: ",
  nrow(CELL),
  "\n",
  "  seeds: ",
  sum(
    CELL$AF_SCAVENGE_seed
  ),
  "\n",
  "  scale factor: ",
  signif(
    SCALE_FACTOR,
    6
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 6. TF-IDF -> 30-D LSI (CHECKPOINT)
# ==============================================================================

if (file.exists(
  LSI_FILE
)) {

  LSI <- readRDS(
    LSI_FILE
  )

  lsi_valid <- (
    is.matrix(
      LSI
    ) &&
      nrow(
        LSI
      ) ==
        EXPECTED_N_CELLS &&
      ncol(
        LSI
      ) ==
        LSI_DIMS &&
      identical(
        rownames(
          LSI
        ),
        CELL$cell_id
      ) &&
      all(
        is.finite(
          LSI
        )
      )
  )

} else {

  lsi_valid <- FALSE
}

if (!lsi_valid) {

  cat(
    "\n============================================================\n",
    "TF-IDF -> LSI30\n",
    "============================================================\n",
    sep = ""
  )

  TFIDF <- sc_tfidf_source_equiv(
    bmat = COUNTS,
    mat_binary = TRUE,
    TF = TRUE,
    log_TF = TRUE
  )

  set.seed(
    9527
  )

  LSI <- sc_do_lsi_source_equiv(
    mat = TFIDF,
    dims = LSI_DIMS
  )

  LSI_METHOD <-
    "SCAVENGE official-source-equivalent local tfidf/do_lsi"

  if (
    nrow(LSI) !=
      EXPECTED_N_CELLS ||
    ncol(LSI) !=
      LSI_DIMS
  ) {
    stop(
      "LSI dimensions are incorrect."
    )
  }

  if (is.null(
    rownames(
      LSI
    )
  )) {

    rownames(
      LSI
    ) <- CELL$cell_id
  }

  if (!identical(
    rownames(
      LSI
    ),
    CELL$cell_id
  )) {
    stop(
      "LSI cell order differs from frozen cell order."
    )
  }

  if (!all(
    is.finite(
      LSI
    )
  )) {
    stop(
      "Raw LSI contains non-finite values."
    )
  }

  saveRDS(
    LSI,
    LSI_FILE,
    compress = TRUE
  )

  rm(
    TFIDF
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )

} else {

  LSI_METHOD <-
    "checkpoint_reused"
}

# ==============================================================================
# 7. HARMONY DONOR CORRECTION (CHECKPOINT)
# ==============================================================================

if (file.exists(
  HARMONY_FILE
)) {

  HARMONY <- readRDS(
    HARMONY_FILE
  )

  harmony_valid <- (
    is.matrix(
      HARMONY
    ) &&
      nrow(
        HARMONY
      ) ==
        EXPECTED_N_CELLS &&
      ncol(
        HARMONY
      ) ==
        LSI_DIMS &&
      identical(
        rownames(
          HARMONY
        ),
        CELL$cell_id
      ) &&
      all(
        is.finite(
          HARMONY
        )
      )
  )

} else {

  harmony_valid <- FALSE
}

if (!harmony_valid) {

  cat(
    "\n============================================================\n",
    "HARMONY — DONOR ONLY\n",
    "============================================================\n",
    sep = ""
  )

  HARMONY_META <- data.frame(
    donor = factor(
      CELL$donor
    ),
    row.names = CELL$cell_id
  )

  set.seed(
    9527
  )

  HARMONY <- harmony::RunHarmony(
    data_mat = LSI,
    meta_data = HARMONY_META,
    vars_use = HARMONY_VARIABLE,
    return_object = FALSE,
    ncores = HARMONY_NCORES,
    verbose = TRUE
  )

  HARMONY <- as.matrix(
    HARMONY
  )

  rownames(
    HARMONY
  ) <- CELL$cell_id

  colnames(
    HARMONY
  ) <- paste0(
    "Harmony_LSI_",
    seq_len(
      ncol(
        HARMONY
      )
    )
  )

  if (
    nrow(HARMONY) !=
      EXPECTED_N_CELLS ||
    ncol(HARMONY) !=
      LSI_DIMS ||
    !all(
      is.finite(
        HARMONY
      )
    )
  ) {
    stop(
      "Harmony output failed dimension/finite-value QC."
    )
  }

  saveRDS(
    HARMONY,
    HARMONY_FILE,
    compress = TRUE
  )
}

# Donor explanatory R2 in each dimension before and after Harmony.
donor_r2 <- function(
  embedding,
  donor
) {

  vapply(
    seq_len(
      ncol(
        embedding
      )
    ),
    function(j) {

      fit <- lm(
        embedding[
          ,
          j
        ] ~
          factor(
            donor
          )
      )

      max(
        0,
        summary(
          fit
        )$adj.r.squared
      )
    },
    numeric(1)
  )
}

R2_BEFORE <- donor_r2(
  LSI,
  CELL$donor
)

R2_AFTER <- donor_r2(
  HARMONY,
  CELL$donor
)

HARMONY_QC <- data.table(
  dimension =
    seq_len(
      LSI_DIMS
    ),
  donor_adjR2_before =
    R2_BEFORE,
  donor_adjR2_after =
    R2_AFTER
)

fwrite(
  HARMONY_QC,
  HARMONY_QC_FILE
)

# ==============================================================================
# 8. OPTIONAL UMAP FOR FINAL VISUALIZATION
# ==============================================================================

UMAP <- NULL

if (RUN_UMAP) {

  if (file.exists(
    UMAP_FILE
  )) {

    UMAP <- readRDS(
      UMAP_FILE
    )

    if (
      !is.matrix(
        UMAP
      ) ||
      nrow(
        UMAP
      ) !=
        EXPECTED_N_CELLS ||
      ncol(
        UMAP
      ) !=
        2L
    ) {
      UMAP <- NULL
    }
  }

  if (is.null(
    UMAP
  )) {

    cat(
      "\nComputing Harmony-UMAP for visualization only...\n"
    )

    set.seed(
      9527
    )

    UMAP <- uwot::umap(
      X = HARMONY,
      n_neighbors =
        UMAP_N_NEIGHBORS,
      min_dist =
        UMAP_MIN_DIST,
      metric = "cosine",
      n_components = 2L,
      n_threads = 1L,
      verbose = TRUE
    )

    rownames(
      UMAP
    ) <- CELL$cell_id

    colnames(
      UMAP
    ) <- c(
      "UMAP1",
      "UMAP2"
    )

    saveRDS(
      UMAP,
      UMAP_FILE,
      compress = TRUE
    )
  }
}

# ==============================================================================
# 9. MUTUAL kNN GRAPH k=30 (CHECKPOINT)
# ==============================================================================

if (file.exists(
  GRAPH_FILE
)) {

  MKNN <- readRDS(
    GRAPH_FILE
  )

  graph_valid <- (
    inherits(
      MKNN,
      "sparseMatrix"
    ) &&
      nrow(
        MKNN
      ) ==
        EXPECTED_N_CELLS &&
      ncol(
        MKNN
      ) ==
        EXPECTED_N_CELLS &&
      identical(
        rownames(
          MKNN
        ),
        CELL$cell_id
      )
  )

} else {

  graph_valid <- FALSE
}

if (!graph_valid) {

  cat(
    "\n============================================================\n",
    "MUTUAL kNN GRAPH — k=30\n",
    "============================================================\n",
    sep = ""
  )

  MKNN <- sc_getmutualknn_source_equiv(
    lsimat = HARMONY,
    num_k = K_MKNN
  )

  GRAPH_METHOD <-
    "SCAVENGE official-source-equivalent local getmutualknn"

  MKNN <- methods::as(
    MKNN,
    "dgCMatrix"
  )

  rownames(
    MKNN
  ) <- CELL$cell_id

  colnames(
    MKNN
  ) <- CELL$cell_id

  Matrix::diag(
    MKNN
  ) <- 0

  MKNN@x[] <- 1

  saveRDS(
    MKNN,
    GRAPH_FILE,
    compress = TRUE
  )

} else {

  GRAPH_METHOD <-
    "checkpoint_reused"
}

DEGREE <- as.numeric(
  Matrix::colSums(
    MKNN
  )
)

N_ISOLATED <- sum(
  DEGREE ==
    0
)

G_IGRAPH <- igraph::graph_from_adjacency_matrix(
  adjmatrix = MKNN,
  mode = "undirected",
  diag = FALSE
)

COMP <- igraph::components(
  G_IGRAPH
)

N_COMPONENTS <- COMP$no

LARGEST_COMPONENT <- max(
  COMP$csize
)

LARGEST_COMPONENT_FRACTION <-
  LARGEST_COMPONENT /
  EXPECTED_N_CELLS

N_EDGES <- igraph::gsize(
  G_IGRAPH
)

SEED_FULL <- as.logical(
  CELL$AF_SCAVENGE_seed
)

COMP_HAS_SEED <- tapply(
  SEED_FULL,
  COMP$membership,
  any
)

CELL_COMP_HAS_SEED <- COMP_HAS_SEED[
  as.character(
    COMP$membership
  )
]

GRAPH_QC <- data.table(
  metric = c(
    "n_cells",
    "k_requested",
    "n_edges_undirected",
    "min_degree",
    "median_degree",
    "mean_degree",
    "max_degree",
    "n_isolated_cells",
    "n_connected_components",
    "largest_component_cells",
    "largest_component_fraction",
    "n_components_with_seed",
    "n_cells_in_components_with_seed",
    "graph_method"
  ),
  value = c(
    EXPECTED_N_CELLS,
    K_MKNN,
    N_EDGES,
    min(
      DEGREE
    ),
    median(
      DEGREE
    ),
    mean(
      DEGREE
    ),
    max(
      DEGREE
    ),
    N_ISOLATED,
    N_COMPONENTS,
    LARGEST_COMPONENT,
    LARGEST_COMPONENT_FRACTION,
    sum(
      COMP_HAS_SEED
    ),
    sum(
      CELL_COMP_HAS_SEED
    ),
    GRAPH_METHOD
  )
)

fwrite(
  GRAPH_QC,
  GRAPH_QC_FILE
)

# ==============================================================================
# 10. REAL NETWORK PROPAGATION + TRS (CHECKPOINT)
# ==============================================================================

if (file.exists(
  REAL_RWR_FILE
)) {

  REAL <- readRDS(
    REAL_RWR_FILE
  )

  real_valid <- (
    is.list(
      REAL
    ) &&
      length(
        REAL$np_score_full
      ) ==
        EXPECTED_N_CELLS &&
      length(
        REAL$TRS_full
      ) ==
        EXPECTED_N_CELLS
  )

} else {

  real_valid <- FALSE
}

if (!real_valid) {

  cat(
    "\n============================================================\n",
    "REAL RWR — gamma=0.05\n",
    "============================================================\n",
    sep = ""
  )

  # Isolated cells cannot participate in a transition matrix.
  RWR_CANDIDATE <- DEGREE >
    0

  MKNN_RWR <- MKNN[
    RWR_CANDIDATE,
    RWR_CANDIDATE,
    drop = FALSE
  ]

  SEED_RWR <- SEED_FULL[
    RWR_CANDIDATE
  ]

  if (
    sum(
      SEED_RWR
    ) <=
      0
  ) {
    stop(
      "No seed cells remain in the non-isolated graph."
    )
  }

  W_RWR <- make_transition(
    MKNN_RWR
  )

  P0_REAL <- matrix(
    0,
    nrow = nrow(
      MKNN_RWR
    ),
    ncol = 1L
  )

  P0_REAL[
    SEED_RWR,
    1
  ] <- 1 /
    sum(
      SEED_RWR
    )

  REAL_RWR <- rwr_matrix(
    W = W_RWR,
    P0 = P0_REAL,
    gamma = RW_GAMMA,
    stationary_cutoff =
      RW_STATIONARY_CUTOFF,
    max_iter =
      RW_MAX_ITER,
    verbose = TRUE
  )

  NP_RWR <- as.numeric(
    REAL_RWR$score[
      ,
      1
    ]
  )

  NP_FULL <- numeric(
    EXPECTED_N_CELLS
  )

  NP_FULL[
    RWR_CANDIDATE
  ] <- NP_RWR

  # Official vignette removes cells with NP==0 (typically cells in
  # disconnected components not containing a seed), then assigns them
  # a unified score of 0 when recovering the full set.
  KEEP_NETWORK <- is.finite(
    NP_FULL
  ) &
    NP_FULL >
      0

  OMIT_NETWORK <- !KEEP_NETWORK

  if (
    sum(
      KEEP_NETWORK
    ) <
      0.90 *
        EXPECTED_N_CELLS
  ) {
    stop(
      "More than 10% of cells have zero/non-finite network propagation score."
    )
  }

  NP_KEEP <- NP_FULL[
    KEEP_NETWORK
  ]

  TRS_KEEP <-
    sc_cap_outlier_source_equiv(
      x = NP_KEEP,
      q_ceiling =
        TRS_CAP_QUANTILE
    )

  TRS_KEEP <-
    sc_max_min_source_equiv(
      TRS_KEEP
    )

  TRS_METHOD <-
    "SCAVENGE official-source-equivalent local capOutlierQuantile/max_min_scale"

  TRS_KEEP <- TRS_KEEP *
    SCALE_FACTOR

  TRS_FULL <- numeric(
    EXPECTED_N_CELLS
  )

  TRS_FULL[
    KEEP_NETWORK
  ] <- TRS_KEEP

  REAL <- list(
    np_score_full =
      NP_FULL,
    TRS_full =
      TRS_FULL,
    keep_network =
      KEEP_NETWORK,
    omit_network =
      OMIT_NETWORK,
    rwr_candidate =
      RWR_CANDIDATE,
    rwr_iterations =
      REAL_RWR$iterations,
    rwr_delta =
      REAL_RWR$delta,
    TRS_method =
      TRS_METHOD
  )

  saveRDS(
    REAL,
    REAL_RWR_FILE,
    compress = TRUE
  )
}

NP_FULL <- REAL$np_score_full
TRS_FULL <- REAL$TRS_full
KEEP_NETWORK <- REAL$keep_network
OMIT_NETWORK <- REAL$omit_network

N_OMITTED <- sum(
  OMIT_NETWORK
)

# ==============================================================================
# 11. PREPARE PRUNED GRAPH FOR DEGREE-MATCHED PERMUTATION
# ==============================================================================

MKNN_KEEP <- MKNN[
  KEEP_NETWORK,
  KEEP_NETWORK,
  drop = FALSE
]

CELL_KEEP <- CELL[
  KEEP_NETWORK
]

SEED_KEEP <- SEED_FULL[
  KEEP_NETWORK
]

NP_KEEP <- NP_FULL[
  KEEP_NETWORK
]

DEGREE_KEEP <- as.integer(
  Matrix::colSums(
    MKNN_KEEP
  )
)

if (any(
  DEGREE_KEEP <= 0
)) {
  stop(
    "Pruned network still contains degree-zero cells."
  )
}

if (
  sum(
    SEED_KEEP
  ) <=
    0
) {
  stop(
    "Pruned network has no seed cells."
  )
}

W_KEEP <- make_transition(
  MKNN_KEEP
)

# Degree-matched seed frequencies, exactly following SCAVENGE's
# permutation logic.
SEED_DEGREE_FREQ <- table(
  DEGREE_KEEP[
    SEED_KEEP
  ]
)

ALL_BY_DEGREE <- split(
  seq_along(
    DEGREE_KEEP
  ),
  DEGREE_KEEP
)

for (deg_name in names(
  SEED_DEGREE_FREQ
)) {

  n_need <- as.integer(
    SEED_DEGREE_FREQ[
      deg_name
    ]
  )

  n_have <- length(
    ALL_BY_DEGREE[[deg_name]]
  )

  if (
    is.null(
      n_have
    ) ||
    n_have <
      n_need
  ) {
    stop(
      "Degree-matched permutation impossible for degree ",
      deg_name
    )
  }
}

# ==============================================================================
# 12. 1000 DEGREE-MATCHED PERMUTATIONS — CHECKPOINTED BATCH RWR
# ==============================================================================

cat(
  "\n============================================================\n",
  "DEGREE-MATCHED PERMUTATION TEST — ",
  PERMUTATION_TIMES,
  " PERMUTATIONS\n",
  "============================================================\n",
  sep = ""
)

N_KEEP <- sum(
  KEEP_NETWORK
)

N_BATCHES <- ceiling(
  PERMUTATION_TIMES /
    PERM_BATCH_SIZE
)

set.seed(
  9527
)

PERM_SEEDS <- sample.int(
  .Machine$integer.max,
  PERMUTATION_TIMES,
  replace = FALSE
)

for (batch_id in seq_len(
  N_BATCHES
)) {

  first_perm <- (
    batch_id -
      1L
  ) *
    PERM_BATCH_SIZE +
    1L

  last_perm <- min(
    batch_id *
      PERM_BATCH_SIZE,
    PERMUTATION_TIMES
  )

  perm_ids <- first_perm:last_perm

  BATCH_FILE <- file.path(
    PERM_DIR,
    sprintf(
      "STEP11G3B_perm_batch_%02d_%04d-%04d.rds",
      batch_id,
      first_perm,
      last_perm
    )
  )

  if (file.exists(
    BATCH_FILE
  )) {

    bobj <- tryCatch(
      readRDS(
        BATCH_FILE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(
        bobj
      ) &&
      identical(
        bobj$perm_ids,
        perm_ids
      ) &&
      length(
        bobj$exceed_gt_count
      ) ==
        N_KEEP &&
      length(
        bobj$exceed_ge_count
      ) ==
        N_KEEP
    ) {

      cat(
        "Permutation batch ",
        batch_id,
        "/",
        N_BATCHES,
        " checkpoint PASS — skip\n",
        sep = ""
      )

      next
    }
  }

  cat(
    "\nPermutation batch ",
    batch_id,
    "/",
    N_BATCHES,
    " [",
    first_perm,
    "-",
    last_perm,
    "]\n",
    sep = ""
  )

  B <- length(
    perm_ids
  )

  P0 <- matrix(
    0,
    nrow = N_KEEP,
    ncol = B
  )

  for (jj in seq_along(
    perm_ids
  )) {

    pp <- perm_ids[
      jj
    ]

    set.seed(
      PERM_SEEDS[
        pp
      ]
    )

    sampled <- integer(
      0
    )

    for (deg_name in names(
      SEED_DEGREE_FREQ
    )) {

      pool <- ALL_BY_DEGREE[[deg_name]]

      n_need <- as.integer(
        SEED_DEGREE_FREQ[
          deg_name
        ]
      )

      sampled <- c(
        sampled,
        sample(
          pool,
          size = n_need,
          replace = FALSE
        )
      )
    }

    sampled <- unique(
      sampled
    )

    if (
      length(
        sampled
      ) !=
        sum(
          SEED_KEEP
        )
    ) {
      stop(
        "Degree-matched sampled seed count differs from real seed count."
      )
    }

    P0[
      sampled,
      jj
    ] <- 1 /
      length(
        sampled
      )
  }

  PR <- rwr_matrix(
    W = W_KEEP,
    P0 = P0,
    gamma = RW_GAMMA,
    stationary_cutoff =
      RW_STATIONARY_CUTOFF,
    max_iter =
      RW_MAX_ITER,
    verbose = FALSE
  )

  PERM_SCORE <- as.matrix(
    PR$score
  )

  EXCEED_GT <- rowSums(
    sweep(
      PERM_SCORE,
      1,
      NP_KEEP,
      FUN = ">"
    )
  )

  EXCEED_GE <- rowSums(
    sweep(
      PERM_SCORE,
      1,
      NP_KEEP,
      FUN = ">="
    )
  )

  bobj <- list(
    batch_id =
      batch_id,
    perm_ids =
      perm_ids,
    exceed_gt_count =
      as.integer(
        EXCEED_GT
      ),
    exceed_ge_count =
      as.integer(
        EXCEED_GE
      ),
    rwr_iterations =
      PR$iterations,
    rwr_delta =
      PR$delta
  )

  saveRDS(
    bobj,
    BATCH_FILE,
    compress = TRUE
  )

  rm(
    P0,
    PR,
    PERM_SCORE,
    EXCEED_GT,
    EXCEED_GE,
    bobj
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )
}

# ==============================================================================
# 13. COMBINE PERMUTATION CHECKPOINTS
# ==============================================================================

EXCEED_GT_TOTAL <- integer(
  N_KEEP
)

EXCEED_GE_TOTAL <- integer(
  N_KEEP
)

PERM_SEEN <- integer(
  0
)

PERM_BATCH_QC <- list()

for (batch_id in seq_len(
  N_BATCHES
)) {

  first_perm <- (
    batch_id -
      1L
  ) *
    PERM_BATCH_SIZE +
    1L

  last_perm <- min(
    batch_id *
      PERM_BATCH_SIZE,
    PERMUTATION_TIMES
  )

  BATCH_FILE <- file.path(
    PERM_DIR,
    sprintf(
      "STEP11G3B_perm_batch_%02d_%04d-%04d.rds",
      batch_id,
      first_perm,
      last_perm
    )
  )

  if (!file.exists(
    BATCH_FILE
  )) {
    stop(
      "Missing permutation batch checkpoint: ",
      BATCH_FILE
    )
  }

  bobj <- readRDS(
    BATCH_FILE
  )

  EXCEED_GT_TOTAL <- EXCEED_GT_TOTAL +
    as.integer(
      bobj$exceed_gt_count
    )

  EXCEED_GE_TOTAL <- EXCEED_GE_TOTAL +
    as.integer(
      bobj$exceed_ge_count
    )

  PERM_SEEN <- c(
    PERM_SEEN,
    bobj$perm_ids
  )

  PERM_BATCH_QC[[
    length(
      PERM_BATCH_QC
    ) +
      1L
  ]] <- data.table(
    batch_id =
      batch_id,
    first_perm =
      min(
        bobj$perm_ids
      ),
    last_perm =
      max(
        bobj$perm_ids
      ),
    n_permutations =
      length(
        bobj$perm_ids
      ),
    rwr_iterations =
      bobj$rwr_iterations,
    rwr_delta =
      bobj$rwr_delta
  )
}

PERM_BATCH_QC <- rbindlist(
  PERM_BATCH_QC
)

fwrite(
  PERM_BATCH_QC,
  file.path(
    QC_DIR,
    "STEP11G3B_permutation_batch_QC.csv"
  )
)

if (
  length(
    unique(
      PERM_SEEN
    )
  ) !=
    PERMUTATION_TIMES ||
  !setequal(
    PERM_SEEN,
    seq_len(
      PERMUTATION_TIMES
    )
  )
) {
  stop(
    "Permutation checkpoint coverage is not exactly 1:1000."
  )
}

# ------------------------------------------------------------------
# Primary classification = exact public SCAVENGE package implementation:
#   rowSums(null > real) <= alpha * B
#
# For a reportable non-zero P value, floor at 1/B.
# ------------------------------------------------------------------

EMPIRICAL_P_KEEP <- pmax(
  EXCEED_GT_TOTAL /
    PERMUTATION_TIMES,
  1 /
    PERMUTATION_TIMES
)

SIG_KEEP <- EXCEED_GT_TOTAL <=
  (
    TRUE_CELL_ALPHA *
      PERMUTATION_TIMES
  )

# ------------------------------------------------------------------
# Sensitivity = formula printed in the Nat Biotechnol paper:
#   (1 + sum(null >= real)) / (1 + B)
# and enriched cells use P < 0.05.
# This does NOT replace the package-code primary result.
# ------------------------------------------------------------------

PAPER_P_KEEP <- (
  1 +
    EXCEED_GE_TOTAL
) /
  (
    1 +
      PERMUTATION_TIMES
  )

PAPER_SIG_KEEP <- PAPER_P_KEEP <
  TRUE_CELL_ALPHA

# Full-cell recovery: omitted cells are unified to 0 score and not
# considered significant, matching the vignette's recommendation.
EMPIRICAL_P_FULL <- rep(
  NA_real_,
  EXPECTED_N_CELLS
)

EMPIRICAL_P_FULL[
  KEEP_NETWORK
] <- EMPIRICAL_P_KEEP

SIG_FULL <- rep(
  FALSE,
  EXPECTED_N_CELLS
)

SIG_FULL[
  KEEP_NETWORK
] <- SIG_KEEP

PAPER_P_FULL <- rep(
  NA_real_,
  EXPECTED_N_CELLS
)

PAPER_P_FULL[
  KEEP_NETWORK
] <- PAPER_P_KEEP

PAPER_SIG_FULL <- rep(
  FALSE,
  EXPECTED_N_CELLS
)

PAPER_SIG_FULL[
  KEEP_NETWORK
] <- PAPER_SIG_KEEP

# ==============================================================================
# 14. FINAL CELL-LEVEL RESULT
# ==============================================================================

FINAL <- copy(
  CELL
)

FINAL[
  ,
  `:=`(
    graph_degree =
      DEGREE,
    graph_component =
      as.integer(
        COMP$membership
      ),
    network_propagation_score =
      NP_FULL,
    AF_TRS =
      TRS_FULL,
    network_eligible =
      KEEP_NETWORK,
    network_omitted_zero_score =
      OMIT_NETWORK,
    permutation_empirical_P =
      EMPIRICAL_P_FULL,
    AF_SCAVENGE_significant_cell =
      SIG_FULL,
    paper_formula_empirical_P =
      PAPER_P_FULL,
    AF_SCAVENGE_significant_cell_paper_formula =
      PAPER_SIG_FULL
  )
]

if (!is.null(
  UMAP
)) {

  FINAL[
    ,
    `:=`(
      UMAP1 =
        UMAP[
          ,
          1
        ],
      UMAP2 =
        UMAP[
          ,
          2
        ]
    )
  ]
}

FINAL[
  ,
  AF_TRS_rank :=
    frank(
      -AF_TRS,
      ties.method = "average"
    )
]

FINAL[
  ,
  AF_TRS_percentile :=
    1 -
    (
      AF_TRS_rank -
        1
    ) /
    (
      .N -
        1
    )
]

fwrite(
  FINAL,
  CELL_RESULT_FILE,
  compress = "gzip"
)

# ==============================================================================
# 15. CELL-TYPE SUMMARY + FORMAL ENRICHMENT OF SIGNIFICANT CELLS
# ==============================================================================

CELLTYPE_SUMMARY <- FINAL[
  ,
  .(
    n_cells =
      .N,
    n_seed =
      sum(
        AF_SCAVENGE_seed
      ),
    seed_fraction =
      mean(
        AF_SCAVENGE_seed
      ),
    n_significant =
      sum(
        AF_SCAVENGE_significant_cell
      ),
    significant_fraction =
      mean(
        AF_SCAVENGE_significant_cell
      ),
    mean_gchromVAR_Z =
      mean(
        AF_gchromVAR_Z,
        na.rm = TRUE
      ),
    median_gchromVAR_Z =
      median(
        AF_gchromVAR_Z,
        na.rm = TRUE
      ),
    mean_TRS =
      mean(
        AF_TRS,
        na.rm = TRUE
      ),
    median_TRS =
      median(
        AF_TRS,
        na.rm = TRUE
      ),
    q90_TRS =
      as.numeric(
        quantile(
          AF_TRS,
          0.90,
          na.rm = TRUE
        )
      ),
    q95_TRS =
      as.numeric(
        quantile(
          AF_TRS,
          0.95,
          na.rm = TRUE
        )
      )
  ),
  by = cell_type
]

setorder(
  CELLTYPE_SUMMARY,
  -significant_fraction,
  -median_TRS
)

fwrite(
  CELLTYPE_SUMMARY,
  CELLTYPE_RESULT_FILE
)

# One-vs-rest Fisher enrichment of permutation-significant cells.
ENRICH_LIST <- lapply(
  unique(
    FINAL$cell_type
  ),
  function(ct) {

    in_ct <- FINAL$cell_type ==
      ct

    sig <- FINAL$AF_SCAVENGE_significant_cell

    a <- sum(
      in_ct &
        sig
    )

    b <- sum(
      in_ct &
        !sig
    )

    c <- sum(
      !in_ct &
        sig
    )

    d <- sum(
      !in_ct &
        !sig
    )

    ft <- fisher.test(
      matrix(
        c(
          a,
          b,
          c,
          d
        ),
        nrow = 2,
        byrow = TRUE
      ),
      alternative = "greater"
    )

    data.table(
      cell_type =
        ct,
      n_cells =
        a +
        b,
      n_significant =
        a,
      significant_fraction =
        a /
        (
          a +
            b
        ),
      odds_ratio =
        unname(
          ft$estimate
        ),
      fisher_P =
        ft$p.value
    )
  }
)

CELLTYPE_ENRICH <- rbindlist(
  ENRICH_LIST
)

CELLTYPE_ENRICH[
  ,
  FDR_BH :=
    p.adjust(
      fisher_P,
      method = "BH"
    )
]

setorder(
  CELLTYPE_ENRICH,
  fisher_P
)

fwrite(
  CELLTYPE_ENRICH,
  CELLTYPE_ENRICH_FILE
)

# ==============================================================================
# 16. DONOR CONSISTENCY + RHYTHM DESCRIPTIVE SUMMARY
# ==============================================================================

DONOR_CELLTYPE <- FINAL[
  ,
  .(
    n_cells =
      .N,
    n_seed =
      sum(
        AF_SCAVENGE_seed
      ),
    seed_fraction =
      mean(
        AF_SCAVENGE_seed
      ),
    n_significant =
      sum(
        AF_SCAVENGE_significant_cell
      ),
    significant_fraction =
      mean(
        AF_SCAVENGE_significant_cell
      ),
    mean_TRS =
      mean(
        AF_TRS
      ),
    median_TRS =
      median(
        AF_TRS
      ),
    mean_Z =
      mean(
        AF_gchromVAR_Z
      ),
    median_Z =
      median(
        AF_gchromVAR_Z
      )
  ),
  by = .(
    donor,
    rhythm,
    cell_type
  )
]

fwrite(
  DONOR_CELLTYPE,
  DONOR_CELLTYPE_FILE
)

# Descriptive only: only 7 donors, and rhythm is donor-level.
# Do not use per-cell pseudo-replication for rhythm inference.
RHYTHM_CELLTYPE <- FINAL[
  ,
  .(
    n_cells =
      .N,
    n_donors =
      uniqueN(
        donor
      ),
    seed_fraction =
      mean(
        AF_SCAVENGE_seed
      ),
    significant_fraction =
      mean(
        AF_SCAVENGE_significant_cell
      ),
    mean_TRS =
      mean(
        AF_TRS
      ),
    median_TRS =
      median(
        AF_TRS
      )
  ),
  by = .(
    rhythm,
    cell_type
  )
]

fwrite(
  RHYTHM_CELLTYPE,
  RHYTHM_CELLTYPE_FILE
)

# ==============================================================================
# 17. PRIMARY QC
# ==============================================================================

N_SIG <- sum(
  FINAL$AF_SCAVENGE_significant_cell
)

N_SIG_PAPER <- sum(
  FINAL$AF_SCAVENGE_significant_cell_paper_formula
)

PACKAGE_PAPER_AGREEMENT <- mean(
  FINAL$AF_SCAVENGE_significant_cell ==
    FINAL$AF_SCAVENGE_significant_cell_paper_formula
)

PACKAGE_PAPER_JACCARD <- {
  a <- FINAL$AF_SCAVENGE_significant_cell
  b <- FINAL$AF_SCAVENGE_significant_cell_paper_formula
  den <- sum(a | b)
  if (den == 0) NA_real_ else sum(a & b) / den
}

SIG_FRACTION <- mean(
  FINAL$AF_SCAVENGE_significant_cell
)

SEED_RECOVERY <- sum(
  FINAL$AF_SCAVENGE_significant_cell &
    FINAL$AF_SCAVENGE_seed
) /
  sum(
    FINAL$AF_SCAVENGE_seed
  )

FOLD_SIG_OVER_SEED <-
  N_SIG /
  sum(
    FINAL$AF_SCAVENGE_seed
  )

CM_DONOR <- DONOR_CELLTYPE[
  cell_type ==
    "CM"
]

CM_ALL_DONORS_PRESENT <- (
  nrow(
    CM_DONOR
  ) ==
    7L
)

CM_MEDIAN_TRS_POSITIVE_ALL_DONORS <- (
  CM_ALL_DONORS_PRESENT &&
    all(
      CM_DONOR$median_TRS >
        0
    )
)

PRIMARY_SUMMARY <- data.table(
  metric = c(
    "n_cells_total",
    "n_seed_cells",
    "n_network_eligible",
    "n_network_omitted_zero_score",
    "n_permutation_significant_cells",
    "n_significant_cells_paper_formula",
    "package_vs_paper_significance_agreement",
    "package_vs_paper_significance_Jaccard",
    "significant_cell_fraction",
    "seed_cells_recovered_as_significant_fraction",
    "fold_significant_cells_over_seed_cells",
    "scale_factor",
    "rwr_gamma",
    "rwr_iterations_real",
    "rwr_final_delta_real",
    "permutation_times",
    "cell_significance_threshold",
    "CM_all_7_donors_present",
    "CM_median_TRS_positive_all_donors"
  ),
  value = c(
    EXPECTED_N_CELLS,
    sum(
      FINAL$AF_SCAVENGE_seed
    ),
    sum(
      KEEP_NETWORK
    ),
    N_OMITTED,
    N_SIG,
    N_SIG_PAPER,
    PACKAGE_PAPER_AGREEMENT,
    PACKAGE_PAPER_JACCARD,
    SIG_FRACTION,
    SEED_RECOVERY,
    FOLD_SIG_OVER_SEED,
    SCALE_FACTOR,
    RW_GAMMA,
    REAL$rwr_iterations,
    REAL$rwr_delta,
    PERMUTATION_TIMES,
    TRUE_CELL_ALPHA,
    CM_ALL_DONORS_PRESENT,
    CM_MEDIAN_TRS_POSITIVE_ALL_DONORS
  )
)

fwrite(
  PRIMARY_SUMMARY,
  file.path(
    FINAL_DIR,
    "STEP11G3B_AF_SCAVENGE_primary_summary.csv"
  )
)

# ==============================================================================
# 17B. COMPARISON AGAINST PREVIOUS FALLBACK RUN (IF PRESENT)
# ==============================================================================

OLD_V2_RESULT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "05_SCAVENGE_NETWORK",
  "05_FINAL",
  "STEP11G3B_AF_SCAVENGE_cell_results.csv.gz"
)

if (file.exists(
  OLD_V2_RESULT
)) {

  OLD <- fread(
    OLD_V2_RESULT,
    select = c(
      "cell_id",
      "AF_TRS",
      "AF_SCAVENGE_significant_cell"
    )
  )

  setnames(
    OLD,
    c(
      "AF_TRS",
      "AF_SCAVENGE_significant_cell"
    ),
    c(
      "AF_TRS_V2_old",
      "SIG_V2_old"
    )
  )

  CMP <- merge(
    FINAL[
      ,
      .(
        cell_id,
        AF_TRS_V3 =
          AF_TRS,
        SIG_V3 =
          AF_SCAVENGE_significant_cell
      )
    ],
    OLD,
    by = "cell_id",
    all = FALSE
  )

  TRS_SPEARMAN_V2_V3 <- suppressWarnings(
    cor(
      CMP$AF_TRS_V2_old,
      CMP$AF_TRS_V3,
      method = "spearman"
    )
  )

  sig_union <- sum(
    CMP$SIG_V2_old |
      CMP$SIG_V3
  )

  SIG_JACCARD_V2_V3 <- if (
    sig_union >
      0
  ) {
    sum(
      CMP$SIG_V2_old &
        CMP$SIG_V3
    ) /
      sig_union
  } else {
    NA_real_
  }

  V2_V3_QC <- data.table(
    metric = c(
      "n_cells_compared",
      "TRS_Spearman_oldV2_vs_official_equivV3",
      "significant_cell_Jaccard_oldV2_vs_official_equivV3",
      "oldV2_n_significant",
      "official_equivV3_n_significant"
    ),
    value = c(
      nrow(
        CMP
      ),
      TRS_SPEARMAN_V2_V3,
      SIG_JACCARD_V2_V3,
      sum(
        CMP$SIG_V2_old
      ),
      sum(
        CMP$SIG_V3
      )
    )
  )

  fwrite(
    V2_V3_QC,
    file.path(
      QC_DIR,
      "STEP11G3B_V3_oldV2_comparison_QC.csv"
    )
  )

  rm(
    OLD,
    CMP
  )
}

# ==============================================================================
# 18. FIGURES — SINGLE FIGURE FILES, NO LEGEND INSIDE
# ==============================================================================

CT_ORDER <- CELLTYPE_SUMMARY$cell_type

FINAL[
  ,
  cell_type_plot :=
    factor(
      cell_type,
      levels =
        rev(
          CT_ORDER
        )
    )
]

CT_LEVELS <- unique(
  as.character(
    FINAL$cell_type_plot
  )
)

CT_LEVELS <- CT_LEVELS[
  !is.na(
    CT_LEVELS
  )
]

# Stable muted journal-style palette.
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

# Figure 1 — TRS boxplot.
p_trs <- ggplot(
  FINAL,
  aes(
    x = AF_TRS,
    y = cell_type_plot,
    fill = cell_type_plot
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.68,
    linewidth = 0.45
  ) +
  scale_fill_manual(
    values = CT_COLORS
  ) +
  labs(
    x = "AF trait relevance score (TRS)",
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
    "Figure_STEP11G3B_AF_TRS_by_celltype.tiff"
  ),
  p_trs,
  width = 6.5,
  height = 5.3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3B_AF_TRS_by_celltype.pdf"
  ),
  p_trs,
  width = 6.5,
  height = 5.3,
  units = "in"
)

# Figure 2 — permutation-significant fraction.
PLOT_SIG <- copy(
  CELLTYPE_SUMMARY
)

PLOT_SIG[
  ,
  cell_type :=
    factor(
      cell_type,
      levels =
        cell_type[
          order(
            significant_fraction
          )
        ]
    )
]

p_sig <- ggplot(
  PLOT_SIG,
  aes(
    x = significant_fraction,
    y = cell_type,
    fill = cell_type
  )
) +
  geom_col(
    width = 0.70,
    linewidth = 0.35,
    colour = "black"
  ) +
  scale_fill_manual(
    values = CT_COLORS
  ) +
  scale_x_continuous(
    labels = function(x) {
      paste0(
        round(
          100 *
            x,
          1
        ),
        "%"
      )
    },
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    x = "Permutation-significant AF-relevant cells",
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
    "Figure_STEP11G3B_AF_significant_cell_fraction.tiff"
  ),
  p_sig,
  width = 6.5,
  height = 5.3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G3B_AF_significant_cell_fraction.pdf"
  ),
  p_sig,
  width = 6.5,
  height = 5.3,
  units = "in"
)

# UMAP figures — visualization only.
if (
  !is.null(
    UMAP
  )
) {

  fwrite(
    FINAL[
      ,
      .(
        cell_id,
        cell_type,
        donor,
        rhythm,
        UMAP1,
        UMAP2,
        AF_gchromVAR_Z,
        AF_TRS,
        AF_SCAVENGE_significant_cell
      )
    ],
    file.path(
      FINAL_DIR,
      "STEP11G3B_UMAP_figure_source.csv.gz"
    ),
    compress = "gzip"
  )

  p_umap_trs <- ggplot(
    FINAL,
    aes(
      x = UMAP1,
      y = UMAP2,
      colour = AF_TRS
    )
  ) +
    geom_point(
      size = 0.35,
      alpha = 0.80
    ) +
    scale_colour_viridis_c() +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
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
        linewidth = 0.55,
        colour = "black"
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G3B_AF_TRS_UMAP.tiff"
    ),
    p_umap_trs,
    width = 5.5,
    height = 4.8,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G3B_AF_TRS_UMAP.pdf"
    ),
    p_umap_trs,
    width = 5.5,
    height = 4.8,
    units = "in"
  )

  FINAL[
    ,
    sig_plot :=
      factor(
        ifelse(
          AF_SCAVENGE_significant_cell,
          "Significant",
          "Other"
        ),
        levels = c(
          "Other",
          "Significant"
        )
      )
  ]

  p_umap_sig <- ggplot(
    FINAL,
    aes(
      x = UMAP1,
      y = UMAP2
    )
  ) +
    geom_point(
      data = FINAL[
        sig_plot ==
          "Other"
      ],
      size = 0.25,
      alpha = 0.35,
      colour = "grey75"
    ) +
    geom_point(
      data = FINAL[
        sig_plot ==
          "Significant"
      ],
      size = 0.45,
      alpha = 0.90,
      colour = "#B2182B"
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
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
        linewidth = 0.55,
        colour = "black"
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G3B_AF_significant_cells_UMAP.tiff"
    ),
    p_umap_sig,
    width = 5.5,
    height = 4.8,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G3B_AF_significant_cells_UMAP.pdf"
    ),
    p_umap_sig,
    width = 5.5,
    height = 4.8,
    units = "in"
  )
}

# ==============================================================================
# 19. READINESS
# ==============================================================================

CM_ENRICH <- CELLTYPE_ENRICH[
  cell_type ==
    "CM"
]

READINESS <- data.table(
  check = c(
    "STEP11G3A_all_19_pass",
    "11986_cells_aligned",
    "30_LSI_dimensions_created",
    "Harmony_donor_correction_created",
    "Harmony_all_values_finite",
    "mutualKNN_k30_created",
    "graph_has_edges",
    "network_propagation_converged",
    "network_eligible_ge_90pct_cells",
    "TRS_all_finite",
    "1000_permutations_complete",
    "permutation_empirical_P_created",
    "package_vs_paper_significance_agreement_ge_99pct",
    "significant_cells_created",
    "12_celltype_summary_created",
    "12_celltype_Fisher_tests_created",
    "CM_present_in_all_7_donors",
    "final_cell_result_created",
    "STEP11G3B_complete"
  ),
  pass = c(
    nrow(
      READY_G3A
    ) ==
      19L &&
      all(
        READY_G3A$pass %in%
          c(
            TRUE,
            "TRUE",
            1
          )
      ),
    nrow(
      FINAL
    ) ==
      EXPECTED_N_CELLS &&
      identical(
        FINAL$cell_id,
        colnames(
          COUNTS
        )
      ),
    nrow(
      LSI
    ) ==
      EXPECTED_N_CELLS &&
      ncol(
        LSI
      ) ==
        LSI_DIMS,
    nrow(
      HARMONY
    ) ==
      EXPECTED_N_CELLS &&
      ncol(
        HARMONY
      ) ==
        LSI_DIMS,
    all(
      is.finite(
        HARMONY
      )
    ),
    nrow(
      MKNN
    ) ==
      EXPECTED_N_CELLS &&
      ncol(
        MKNN
      ) ==
        EXPECTED_N_CELLS,
    N_EDGES >
      0,
    is.finite(
      REAL$rwr_delta
    ) &&
      REAL$rwr_delta <=
        RW_STATIONARY_CUTOFF,
    mean(
      KEEP_NETWORK
    ) >=
      0.90,
    all(
      is.finite(
        FINAL$AF_TRS
      )
    ),
    length(
      unique(
        PERM_SEEN
      )
    ) ==
      PERMUTATION_TIMES,
    all(
      is.finite(
        FINAL$permutation_empirical_P[
          FINAL$network_eligible
        ]
      )
    ),
    PACKAGE_PAPER_AGREEMENT >=
      0.99,
    N_SIG >
      0,
    nrow(
      CELLTYPE_SUMMARY
    ) ==
      EXPECTED_N_CELLTYPES,
    nrow(
      CELLTYPE_ENRICH
    ) ==
      EXPECTED_N_CELLTYPES,
    CM_ALL_DONORS_PRESENT,
    file.exists(
      CELL_RESULT_FILE
    ),
    TRUE
  )
)

fwrite(
  READINESS,
  READINESS_FILE
)

METHOD <- data.table(
  field = c(
    "dataset",
    "cells",
    "cell_types",
    "TFIDF",
    "LSI_dimensions",
    "batch_correction",
    "Harmony_covariate",
    "rhythm_corrected",
    "mutual_kNN",
    "RWR_gamma",
    "RWR_stationary_cutoff",
    "TRS_cap_quantile",
    "TRS_scale_factor_source",
    "permutation_null",
    "permutation_times",
    "cell_significance_rule",
    "SCAVENGE_package_available",
    "source_equivalence_audit",
    "TFIDF_log_transform",
    "mutualKNN_zero_reciprocal_rule",
    "permutation_primary_rule",
    "permutation_paper_sensitivity_rule",
    "LSI_method",
    "graph_method",
    "TRS_method"
  ),
  value = c(
    "GSE238242 left atrial appendage snATAC",
    as.character(
      EXPECTED_N_CELLS
    ),
    as.character(
      EXPECTED_N_CELLTYPES
    ),
    "binary=TRUE; TF=TRUE; log_TF=TRUE",
    as.character(
      LSI_DIMS
    ),
    "Harmony",
    "donor only",
    "No",
    paste0(
      "mutual kNN k=",
      K_MKNN
    ),
    as.character(
      RW_GAMMA
    ),
    as.character(
      RW_STATIONARY_CUTOFF
    ),
    as.character(
      TRS_CAP_QUANTILE
    ),
    paste0(
      "STEP11G3A top-1% gchromVAR Z mean; scale_factor=",
      signif(
        SCALE_FACTOR,
        6
      )
    ),
    "seed cells randomly resampled within exact graph-degree strata",
    as.character(
      PERMUTATION_TIMES
    ),
    paste0(
      "SCAVENGE convention: count(null NP > real NP)/",
      PERMUTATION_TIMES,
      " <= ",
      TRUE_CELL_ALPHA
    ),
    as.character(
      HAS_SCAVENGE
    ),
    "Local functions transcribed to match current public SCAVENGE R source",
    "log1p(TF * 100000), matching SCAVENGE::tfidf",
    "If no reciprocal neighbour, connect to nearest neighbour, matching SCAVENGE::getmutualknn",
    "rowSums(null > real) <= 0.05*B, matching public get_sigcell_simple.R",
    "(1 + sum(null >= real))/(1+B) < 0.05, Nat Biotechnol formula sensitivity",
    LSI_METHOD,
    GRAPH_METHOD,
    REAL$TRS_method
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G3B_method_provenance.csv"
  )
)

PKG <- data.table(
  package = c(
    "R",
    "Matrix",
    "irlba",
    "RANN",
    "igraph",
    "harmony",
    "SCAVENGE"
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
        "irlba"
      )
    ),
    as.character(
      packageVersion(
        "RANN"
      )
    ),
    as.character(
      packageVersion(
        "igraph"
      )
    ),
    as.character(
      packageVersion(
        "harmony"
      )
    ),
    if (
      HAS_SCAVENGE
    ) {
      as.character(
        packageVersion(
          "SCAVENGE"
        )
      )
    } else {
      "not installed; exact-compatible core fallback used"
    }
  )
)

fwrite(
  PKG,
  file.path(
    QC_DIR,
    "STEP11G3B_package_versions.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G3B_sessionInfo.txt"
  )
)

# ==============================================================================
# 20. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G3B V3 OFFICIAL-SOURCE-EQUIVALENT COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Cells: ",
  EXPECTED_N_CELLS,
  "\n",
  "Seeds: ",
  sum(
    FINAL$AF_SCAVENGE_seed
  ),
  "\n",
  "Graph edges: ",
  format(
    N_EDGES,
    big.mark = ","
  ),
  "\n",
  "Connected components: ",
  N_COMPONENTS,
  "\n",
  "Network-eligible cells: ",
  sum(
    KEEP_NETWORK
  ),
  " (",
  sprintf(
    "%.2f%%",
    100 *
      mean(
        KEEP_NETWORK
      )
  ),
  ")\n",
  "Permutation-significant cells: ",
  N_SIG,
  " (",
  sprintf(
    "%.2f%%",
    100 *
      SIG_FRACTION
  ),
  ")\n",
  "Seed cells recovered as significant: ",
  sprintf(
    "%.2f%%",
    100 *
      SEED_RECOVERY
  ),
  "\n",
  sep = ""
)

cat(
  "\nCell-type summary:\n"
)

print(
  CELLTYPE_SUMMARY
)

cat(
  "\nCell-type Fisher enrichment:\n"
)

print(
  CELLTYPE_ENRICH
)

cat(
  "\nCM donor consistency:\n"
)

print(
  CM_DONOR
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

if (
  all(
    READINESS$pass
  )
) {

  cat(
    "\nSTEP11G3B V3 OFFICIAL-SOURCE-EQUIVALENT PASSED.\n",
    "Upload the entire 00_QC folder plus 05_FINAL results.\n",
    "The next stage will be STEP11G4: biological interpretation, ",
    "CM-state dissection, robustness/sensitivity, and final publication panels.\n",
    sep = ""
  )

} else {

  cat(
    "\nOne or more STEP11G3B readiness checks failed.\n",
    "Do NOT proceed to biological interpretation until reviewed.\n",
    sep = ""
  )
}
