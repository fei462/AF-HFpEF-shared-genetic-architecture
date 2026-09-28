# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G4A_V2_FIX_HITS_NAMESPACE_SENSITIVITY(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G4A V2 — LOW-REFERENCE-MATCH BLOCK SENSITIVITY
# AF SCAVENGE robustness analysis
#
# PURPOSE
#   Re-run only the trait-weighted part of the frozen analysis after excluding
#   the 3 fine-mapping blocks previously flagged with
#   reference_match_fraction < 0.90.
#
# IMPORTANT
#   This DOES NOT rebuild:
#     - peak-by-cell matrix
#     - TF-IDF / LSI
#     - Harmony
#     - mutual-kNN graph
#
#   It reuses the finalized V3 official-source-equivalent network and repeats:
#     PIP weights -> gchromVAR Z -> seed cells -> RWR/TRS ->
#     1000 degree-matched permutations.
#
# PRIMARY QUESTION
#   Does exclusion of the 3 low-reference-match blocks materially change:
#     1) gchromVAR Z ranking?
#     2) seed-cell set?
#     3) TRS ranking?
#     4) significant-cell set?
#     5) CM prioritization?
#
# OUTPUT
#   06_STEP11G4A_LOW_MATCH_SENSITIVITY/
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))
set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

EXPECTED_N_CELLS <- 11986L
EXPECTED_N_PEAKS <- 212084L
EXPECTED_LOW_BLOCKS <- 3L

PIP_MIN <- 0.001
BACKGROUND_ITERATIONS <- 200L

SEED_PERCENT <- 0.05
SCALE_PERCENT <- 0.01

RW_GAMMA <- 0.05
RW_STATIONARY_CUTOFF <- 1e-5
RW_MAX_ITER <- 10000L

TRS_CAP_QUANTILE <- 0.95

PERMUTATION_TIMES <- 1000L
PERM_BATCH_SIZE <- 100L
CELL_ALPHA <- 0.05

# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT <- file.path(DATA_ROOT, "STEP11_GSE238242")

G2B_ROOT <- file.path(
  ROOT, "04_SCAVENGE", "03_AF_SUSIE_FINEMAP"
)

PIP_FILE <- file.path(
  G2B_ROOT, "03_FINAL_PIP",
  "AF_SCAVENGE_PIP_PRIMARY_ALL.tsv.gz"
)

G3A_ROOT <- file.path(
  ROOT, "04_SCAVENGE", "04_GCHROMVAR"
)

LOW_BLOCK_FILE <- file.path(
  G3A_ROOT, "00_QC",
  "STEP11G3A_low_reference_match_blocks_for_sensitivity.csv"
)

SE_GC_FILE <- file.path(
  G3A_ROOT, "03_GCHROMVAR",
  "GSE238242_SCAVENGE_SE_with_GCbias.rds"
)

BACKGROUND_FILE <- file.path(
  G3A_ROOT, "03_GCHROMVAR",
  "GSE238242_gchromVAR_background200.rds"
)

G3B_V3_ROOT <- file.path(
  ROOT, "04_SCAVENGE",
  "05_SCAVENGE_NETWORK_V3_OFFICIAL_EQUIV"
)

GRAPH_FILE <- file.path(
  G3B_V3_ROOT, "02_NETWORK",
  "STEP11G3B_mutualKNN_k30.rds"
)

PRIMARY_CELL_FILE <- file.path(
  G3B_V3_ROOT, "05_FINAL",
  "STEP11G3B_AF_SCAVENGE_cell_results.csv.gz"
)

PRIMARY_READINESS_FILE <- file.path(
  G3B_V3_ROOT, "00_QC",
  "STEP11G3B_readiness.csv"
)

OUT_ROOT <- file.path(
  ROOT, "04_SCAVENGE",
  "06_STEP11G4A_LOW_MATCH_SENSITIVITY"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
GCV_DIR <- file.path(OUT_ROOT, "01_GCHROMVAR")
RWR_DIR <- file.path(OUT_ROOT, "02_RWR_TRS")
PERM_DIR <- file.path(OUT_ROOT, "03_PERMUTATION")
FINAL_DIR <- file.path(OUT_ROOT, "04_FINAL")
FIG_DIR <- file.path(OUT_ROOT, "05_FIGURES")

for (d in c(QC_DIR, GCV_DIR, RWR_DIR, PERM_DIR, FINAL_DIR, FIG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

SENS_DEV_FILE <- file.path(
  GCV_DIR,
  "STEP11G4A_lowmatch_excluded_gchromVAR_deviations.rds"
)

SENS_SCORE_FILE <- file.path(
  GCV_DIR,
  "STEP11G4A_lowmatch_excluded_gchromVAR_scores.csv.gz"
)

REAL_RWR_FILE <- file.path(
  RWR_DIR,
  "STEP11G4A_lowmatch_excluded_real_RWR_TRS.rds"
)

FINAL_CELL_FILE <- file.path(
  FINAL_DIR,
  "STEP11G4A_lowmatch_excluded_cell_results.csv.gz"
)

COMPARISON_FILE <- file.path(
  FINAL_DIR,
  "STEP11G4A_primary_vs_lowmatch_excluded_comparison.csv"
)

CELLTYPE_FILE <- file.path(
  FINAL_DIR,
  "STEP11G4A_lowmatch_excluded_celltype_summary.csv"
)

DONOR_CELLTYPE_FILE <- file.path(
  FINAL_DIR,
  "STEP11G4A_lowmatch_excluded_donor_celltype_summary.csv"
)

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

required_pkgs <- c(
  "data.table",
  "Matrix",
  "SummarizedExperiment",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "gchromVAR"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_pkgs, collapse = ", "),
      "\nThese packages were already required by STEP11G3A. ",
      "Please use the same R library/environment that successfully ran G3A."
    )
  )
}

library(data.table)
library(Matrix)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

jaccard_bool <- function(a, b) {
  a <- as.logical(a)
  b <- as.logical(b)
  den <- sum(a | b, na.rm = TRUE)
  if (den == 0) return(NA_real_)
  sum(a & b, na.rm = TRUE) / den
}

make_seed_index <- function(z_score, percent_cut = 0.05) {

  finite <- is.finite(z_score)

  p_one <- rep(NA_real_, length(z_score))
  p_one[finite] <- pnorm(
    z_score[finite],
    lower.tail = FALSE
  )

  initial <- finite & p_one <= 0.05

  if (
    100 * mean(initial) >
      100 * percent_cut
  ) {

    r <- rank(
      -z_score,
      na.last = "keep"
    )

    seed <- finite &
      r <= (
        percent_cut *
          length(z_score)
      )

  } else {

    seed <- initial
  }

  seed[is.na(seed)] <- FALSE

  list(
    seed = seed,
    initial = initial,
    p_one = p_one
  )
}

make_scale_factor <- function(z_score, percent_cut = 0.01) {

  finite <- is.finite(z_score)
  z <- z_score[finite]

  r <- rank(-z)

  idx <- r <= (
    percent_cut *
      length(z)
  )

  mean(z[idx])
}

make_transition <- function(adj) {

  adj <- methods::as(
    adj,
    "dgCMatrix"
  )

  cs <- Matrix::colSums(adj)

  if (any(!is.finite(cs) | cs <= 0)) {
    stop("Graph contains zero/invalid degree cells.")
  }

  Matrix::t(
    Matrix::t(adj) /
      cs
  )
}

rwr_matrix <- function(
  W,
  P0,
  gamma = 0.05,
  stationary_cutoff = 1e-5,
  max_iter = 10000L
) {

  P <- as.matrix(P0)

  cs <- colSums(P)

  if (any(!is.finite(cs) | cs <= 0)) {
    stop("Every P0 column must have positive seed mass.")
  }

  P <- sweep(P, 2, cs, "/")

  for (iter in seq_len(max_iter)) {

    P_new <- (
      1 - gamma
    ) * (
      W %*% P
    ) + gamma * P0

    delta_each <- colSums(
      abs(
        as.matrix(P_new) -
          P
      )
    )

    delta <- max(delta_each)

    P <- as.matrix(P_new)

    if (
      is.finite(delta) &&
      delta <= stationary_cutoff
    ) {

      return(
        list(
          score = P,
          iterations = iter,
          delta = delta
        )
      )
    }
  }

  stop("RWR did not converge.")
}

cap95_scale <- function(x, scale_factor) {

  q <- as.numeric(
    quantile(
      x,
      probs = TRS_CAP_QUANTILE,
      na.rm = TRUE,
      names = FALSE
    )
  )

  y <- pmin(x, q)

  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)

  if (
    !is.finite(ymin) ||
    !is.finite(ymax) ||
    ymax <= ymin
  ) {
    stop("TRS min-max scaling failed.")
  }

  y <- (
    y - ymin
  ) / (
    ymax - ymin
  )

  y * scale_factor
}

# ==============================================================================
# 4. HARD INPUT CHECK
# ==============================================================================

required_inputs <- c(
  PIP_FILE,
  LOW_BLOCK_FILE,
  SE_GC_FILE,
  BACKGROUND_FILE,
  GRAPH_FILE,
  PRIMARY_CELL_FILE,
  PRIMARY_READINESS_FILE
)

if (!all(file.exists(required_inputs))) {
  stop(
    paste(
      "Missing input(s):",
      paste(
        required_inputs[!file.exists(required_inputs)],
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}

PRIMARY_READY <- fread(PRIMARY_READINESS_FILE)

if (!all(PRIMARY_READY$pass %in% c(TRUE, "TRUE", 1))) {
  stop("Frozen STEP11G3B V3 readiness is not fully PASS.")
}

LOW_BLOCKS <- fread(LOW_BLOCK_FILE)

if (
  !"selected_block_label" %in%
    names(LOW_BLOCKS)
) {
  stop("Low-reference-match block table lacks selected_block_label.")
}

LOW_LABELS <- unique(
  as.character(
    LOW_BLOCKS$selected_block_label
  )
)

if (length(LOW_LABELS) != EXPECTED_LOW_BLOCKS) {
  stop(
    "Expected exactly 3 low-reference-match blocks; observed ",
    length(LOW_LABELS)
  )
}

fwrite(
  LOW_BLOCKS,
  file.path(
    QC_DIR,
    "STEP11G4A_excluded_low_reference_match_blocks.csv"
  )
)

cat(
  "\nExcluded low-reference-match blocks:\n"
)
print(LOW_BLOCKS)

# ==============================================================================
# 5. LOAD FROZEN OBJECTS
# ==============================================================================

SE_GC <- readRDS(SE_GC_FILE)
BG <- readRDS(BACKGROUND_FILE)
MKNN <- readRDS(GRAPH_FILE)
PRIMARY <- fread(PRIMARY_CELL_FILE)

if (
  nrow(SE_GC) != EXPECTED_N_PEAKS ||
  ncol(SE_GC) != EXPECTED_N_CELLS
) {
  stop("SE_GC dimensions mismatch.")
}

if (
  nrow(MKNN) != EXPECTED_N_CELLS ||
  ncol(MKNN) != EXPECTED_N_CELLS
) {
  stop("Frozen V3 graph dimensions mismatch.")
}

if (
  nrow(BG) != EXPECTED_N_PEAKS ||
  ncol(BG) != BACKGROUND_ITERATIONS
) {
  stop("Background peak matrix dimensions mismatch.")
}

if (nrow(PRIMARY) != EXPECTED_N_CELLS) {
  stop("Primary V3 cell result count mismatch.")
}

# Align primary cells to graph.
idx_primary <- match(
  rownames(MKNN),
  PRIMARY$cell_id
)

if (anyNA(idx_primary)) {
  stop("Could not align primary V3 result to frozen graph.")
}

PRIMARY <- PRIMARY[idx_primary]

if (!identical(PRIMARY$cell_id, rownames(MKNN))) {
  stop("Primary V3 cell order alignment failed.")
}

# ==============================================================================
# 6. EXCLUDE LOW-MATCH BLOCKS AND BUILD GRCh38 PIP WEIGHTS
# ==============================================================================

PIP <- fread(
  PIP_FILE,
  select = c(
    "selected_block_label",
    "CHR_HG19",
    "BP_GRCh38",
    "SNP",
    "PIP_primary"
  )
)

PIP_ALL_GE <- PIP[
  is.finite(PIP_primary) &
    PIP_primary >= PIP_MIN &
    is.finite(BP_GRCh38) &
    BP_GRCh38 > 0
]

PIP_SENS <- PIP_ALL_GE[
  !selected_block_label %in%
    LOW_LABELS
]

PIP_QC <- data.table(
  metric = c(
    "PIP_threshold",
    "primary_variants_GE_threshold",
    "sensitivity_variants_GE_threshold",
    "variants_removed",
    "primary_blocks_GE_threshold",
    "sensitivity_blocks_GE_threshold",
    "excluded_blocks",
    "primary_retained_PIP_mass_GE_threshold",
    "sensitivity_retained_PIP_mass_GE_threshold",
    "removed_PIP_mass_GE_threshold",
    "removed_PIP_mass_fraction_of_GE_threshold"
  ),
  value = c(
    PIP_MIN,
    nrow(PIP_ALL_GE),
    nrow(PIP_SENS),
    nrow(PIP_ALL_GE) - nrow(PIP_SENS),
    uniqueN(PIP_ALL_GE$selected_block_label),
    uniqueN(PIP_SENS$selected_block_label),
    length(LOW_LABELS),
    sum(PIP_ALL_GE$PIP_primary),
    sum(PIP_SENS$PIP_primary),
    sum(PIP_ALL_GE$PIP_primary) -
      sum(PIP_SENS$PIP_primary),
    (
      sum(PIP_ALL_GE$PIP_primary) -
        sum(PIP_SENS$PIP_primary)
    ) /
      sum(PIP_ALL_GE$PIP_primary)
  )
)

fwrite(
  PIP_QC,
  file.path(
    QC_DIR,
    "STEP11G4A_PIP_exclusion_QC.csv"
  )
)

PEAK_GR <- SummarizedExperiment::rowRanges(SE_GC)

PIP_GR <- GenomicRanges::GRanges(
  seqnames = paste0(
    "chr",
    as.integer(
      PIP_SENS$CHR_HG19
    )
  ),
  ranges = IRanges::IRanges(
    start = as.integer(
      PIP_SENS$BP_GRCh38
    ),
    width = 1L
  )
)

hits <- GenomicRanges::findOverlaps(
  PEAK_GR,
  PIP_GR,
  ignore.strand = TRUE
)

WEIGHT <- numeric(
  EXPECTED_N_PEAKS
)

if (length(hits)) {

  tmp <- data.table(
    peak_index =
      S4Vectors::queryHits(hits),
    weight =
      PIP_SENS$PIP_primary[
        S4Vectors::subjectHits(hits)
      ]
  )[
    ,
    .(
      weight = sum(weight)
    ),
    by = peak_index
  ]

  WEIGHT[
    tmp$peak_index
  ] <- tmp$weight
}

if (sum(WEIGHT) <= 0) {
  stop("No sensitivity PIP overlaps the frozen ATAC peak universe.")
}

WEIGHT_MAT <- Matrix::Matrix(
  WEIGHT,
  ncol = 1L,
  sparse = TRUE
)

colnames(WEIGHT_MAT) <- "AF_lowmatch_excluded"

OVERLAP_QC <- data.table(
  metric = c(
    "n_sensitivity_PIP_variants",
    "n_variants_overlapping_ATAC_peaks",
    "variant_overlap_fraction",
    "n_weighted_ATAC_peaks",
    "sum_peak_PIP_weight"
  ),
  value = c(
    nrow(PIP_SENS),
    uniqueN(
      S4Vectors::subjectHits(hits)
    ),
    uniqueN(
      S4Vectors::subjectHits(hits)
    ) /
      nrow(PIP_SENS),
    sum(WEIGHT > 0),
    sum(WEIGHT)
  )
)

fwrite(
  OVERLAP_QC,
  file.path(
    QC_DIR,
    "STEP11G4A_PIP_peak_overlap_QC.csv"
  )
)

# ==============================================================================
# 7. gchromVAR SENSITIVITY — REUSE SAME GC OBJECT + SAME 200 BACKGROUNDS
# ==============================================================================

if (file.exists(SENS_DEV_FILE)) {

  DEV <- readRDS(SENS_DEV_FILE)

} else {

  cat(
    "\nRunning gchromVAR after excluding low-match blocks...\n"
  )

  DEV <- gchromVAR::computeWeightedDeviations(
    object = SE_GC,
    weights = methods::as(
      WEIGHT_MAT,
      "dgCMatrix"
    ),
    background_peaks = BG
  )

  saveRDS(
    DEV,
    SENS_DEV_FILE,
    compress = TRUE
  )
}

if (!"z" %in% SummarizedExperiment::assayNames(DEV)) {
  stop("Sensitivity gchromVAR result lacks 'z' assay.")
}

ZMAT <- SummarizedExperiment::assay(
  DEV,
  "z"
)

if (
  nrow(ZMAT) == 1L &&
  ncol(ZMAT) == EXPECTED_N_CELLS
) {

  Z <- as.numeric(
    ZMAT[1, ]
  )

} else if (
  ncol(ZMAT) == 1L &&
  nrow(ZMAT) == EXPECTED_N_CELLS
) {

  Z <- as.numeric(
    ZMAT[, 1]
  )

} else {

  stop(
    "Unexpected sensitivity gchromVAR Z dimensions: ",
    paste(dim(ZMAT), collapse = " x ")
  )
}

if (mean(is.finite(Z)) < 0.99) {
  stop("<99% cells have finite sensitivity gchromVAR Z.")
}

SEED_RES <- make_seed_index(
  Z,
  SEED_PERCENT
)

SEED <- SEED_RES$seed

SCALE_FACTOR <- make_scale_factor(
  Z,
  SCALE_PERCENT
)

SENS_GCV <- data.table(
  cell_id = rownames(MKNN),
  AF_gchromVAR_Z_lowmatch_excluded =
    Z,
  AF_gchromVAR_P_lowmatch_excluded =
    SEED_RES$p_one,
  AF_seed_lowmatch_excluded =
    SEED
)

fwrite(
  SENS_GCV,
  SENS_SCORE_FILE,
  compress = "gzip"
)

# ==============================================================================
# 8. REAL RWR/TRS ON EXACT SAME FROZEN V3 GRAPH
# ==============================================================================

if (file.exists(REAL_RWR_FILE)) {

  REAL <- readRDS(REAL_RWR_FILE)

} else {

  W <- make_transition(MKNN)

  P0 <- matrix(
    0,
    nrow = EXPECTED_N_CELLS,
    ncol = 1L
  )

  P0[
    SEED,
    1
  ] <- 1 / sum(SEED)

  RR <- rwr_matrix(
    W = W,
    P0 = P0,
    gamma = RW_GAMMA,
    stationary_cutoff =
      RW_STATIONARY_CUTOFF,
    max_iter =
      RW_MAX_ITER
  )

  NP <- as.numeric(
    RR$score[, 1]
  )

  KEEP <- is.finite(NP) &
    NP > 0

  if (mean(KEEP) < 0.90) {
    stop("<90% cells are network-eligible in sensitivity analysis.")
  }

  TRS <- numeric(
    EXPECTED_N_CELLS
  )

  TRS[
    KEEP
  ] <- cap95_scale(
    NP[KEEP],
    SCALE_FACTOR
  )

  REAL <- list(
    NP = NP,
    TRS = TRS,
    KEEP = KEEP,
    iterations = RR$iterations,
    delta = RR$delta,
    scale_factor = SCALE_FACTOR,
    n_seed = sum(SEED)
  )

  saveRDS(
    REAL,
    REAL_RWR_FILE,
    compress = TRUE
  )
}

NP <- REAL$NP
TRS <- REAL$TRS
KEEP <- REAL$KEEP

# ==============================================================================
# 9. DEGREE-MATCHED 1000 PERMUTATIONS
# ==============================================================================

MKNN_KEEP <- MKNN[
  KEEP,
  KEEP,
  drop = FALSE
]

SEED_KEEP <- SEED[
  KEEP
]

NP_KEEP <- NP[
  KEEP
]

DEGREE_KEEP <- as.integer(
  Matrix::colSums(
    MKNN_KEEP
  )
)

if (any(DEGREE_KEEP <= 0)) {
  stop("Pruned sensitivity graph contains degree-zero cells.")
}

if (sum(SEED_KEEP) <= 0) {
  stop("No sensitivity seed cell remains in pruned graph.")
}

W_KEEP <- make_transition(
  MKNN_KEEP
)

SEED_DEGREE_FREQ <- table(
  DEGREE_KEEP[
    SEED_KEEP
  ]
)

ALL_BY_DEGREE <- split(
  seq_along(DEGREE_KEEP),
  DEGREE_KEEP
)

for (deg_name in names(SEED_DEGREE_FREQ)) {

  n_need <- as.integer(
    SEED_DEGREE_FREQ[[deg_name]]
  )

  pool <- ALL_BY_DEGREE[[deg_name]]

  if (
    is.null(pool) ||
    length(pool) < n_need
  ) {
    stop(
      "Degree-matched permutation impossible for degree ",
      deg_name
    )
  }
}

N_KEEP <- sum(KEEP)
N_BATCHES <- ceiling(
  PERMUTATION_TIMES /
    PERM_BATCH_SIZE
)

set.seed(9527)

PERM_SEEDS <- sample.int(
  .Machine$integer.max,
  PERMUTATION_TIMES,
  replace = FALSE
)

for (batch_id in seq_len(N_BATCHES)) {

  first_perm <- (
    batch_id - 1L
  ) * PERM_BATCH_SIZE + 1L

  last_perm <- min(
    batch_id *
      PERM_BATCH_SIZE,
    PERMUTATION_TIMES
  )

  perm_ids <- first_perm:last_perm

  BATCH_FILE <- file.path(
    PERM_DIR,
    sprintf(
      "STEP11G4A_perm_batch_%02d_%04d-%04d.rds",
      batch_id,
      first_perm,
      last_perm
    )
  )

  checkpoint_ok <- FALSE

  if (file.exists(BATCH_FILE)) {

    bobj <- tryCatch(
      readRDS(BATCH_FILE),
      error = function(e) NULL
    )

    if (
      !is.null(bobj) &&
      identical(
        bobj$perm_ids,
        perm_ids
      ) &&
      length(
        bobj$exceed_gt_count
      ) == N_KEEP &&
      length(
        bobj$exceed_ge_count
      ) == N_KEEP
    ) {
      checkpoint_ok <- TRUE
    }
  }

  if (checkpoint_ok) {

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

  B <- length(perm_ids)

  P0 <- matrix(
    0,
    nrow = N_KEEP,
    ncol = B
  )

  for (jj in seq_along(perm_ids)) {

    pp <- perm_ids[jj]

    set.seed(
      PERM_SEEDS[pp]
    )

    sampled <- integer(0)

    for (deg_name in names(SEED_DEGREE_FREQ)) {

      pool <- ALL_BY_DEGREE[[deg_name]]

      n_need <- as.integer(
        SEED_DEGREE_FREQ[[deg_name]]
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

    if (
      length(unique(sampled)) !=
        sum(SEED_KEEP)
    ) {
      stop("Sensitivity permuted seed count mismatch.")
    }

    P0[
      sampled,
      jj
    ] <- 1 / length(sampled)
  }

  PR <- rwr_matrix(
    W = W_KEEP,
    P0 = P0,
    gamma = RW_GAMMA,
    stationary_cutoff =
      RW_STATIONARY_CUTOFF,
    max_iter =
      RW_MAX_ITER
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

  saveRDS(
    list(
      perm_ids = perm_ids,
      exceed_gt_count =
        as.integer(EXCEED_GT),
      exceed_ge_count =
        as.integer(EXCEED_GE),
      rwr_iterations =
        PR$iterations,
      rwr_delta =
        PR$delta
    ),
    BATCH_FILE,
    compress = TRUE
  )

  rm(P0, PR, PERM_SCORE)
  invisible(gc(verbose = FALSE))
}

# ==============================================================================
# 10. COMBINE PERMUTATIONS
# ==============================================================================

GT_TOTAL <- integer(N_KEEP)
GE_TOTAL <- integer(N_KEEP)
PERM_SEEN <- integer(0)
PERM_QC_LIST <- list()

for (batch_id in seq_len(N_BATCHES)) {

  first_perm <- (
    batch_id - 1L
  ) * PERM_BATCH_SIZE + 1L

  last_perm <- min(
    batch_id *
      PERM_BATCH_SIZE,
    PERMUTATION_TIMES
  )

  BATCH_FILE <- file.path(
    PERM_DIR,
    sprintf(
      "STEP11G4A_perm_batch_%02d_%04d-%04d.rds",
      batch_id,
      first_perm,
      last_perm
    )
  )

  bobj <- readRDS(BATCH_FILE)

  GT_TOTAL <- GT_TOTAL +
    bobj$exceed_gt_count

  GE_TOTAL <- GE_TOTAL +
    bobj$exceed_ge_count

  PERM_SEEN <- c(
    PERM_SEEN,
    bobj$perm_ids
  )

  PERM_QC_LIST[[length(PERM_QC_LIST) + 1L]] <-
    data.table(
      batch_id = batch_id,
      first_perm = min(bobj$perm_ids),
      last_perm = max(bobj$perm_ids),
      n_permutations =
        length(bobj$perm_ids),
      rwr_iterations =
        bobj$rwr_iterations,
      rwr_delta =
        bobj$rwr_delta
    )
}

if (
  length(unique(PERM_SEEN)) !=
    PERMUTATION_TIMES
) {
  stop("Sensitivity permutation coverage != 1000.")
}

PERM_QC <- rbindlist(PERM_QC_LIST)

fwrite(
  PERM_QC,
  file.path(
    QC_DIR,
    "STEP11G4A_permutation_batch_QC.csv"
  )
)

P_KEEP <- pmax(
  GT_TOTAL /
    PERMUTATION_TIMES,
  1 /
    PERMUTATION_TIMES
)

SIG_KEEP <- GT_TOTAL <=
  (
    CELL_ALPHA *
      PERMUTATION_TIMES
  )

PAPER_P_KEEP <- (
  1 + GE_TOTAL
) /
  (
    1 + PERMUTATION_TIMES
  )

PAPER_SIG_KEEP <- PAPER_P_KEEP <
  CELL_ALPHA

P_FULL <- rep(
  NA_real_,
  EXPECTED_N_CELLS
)

SIG_FULL <- rep(
  FALSE,
  EXPECTED_N_CELLS
)

PAPER_P_FULL <- rep(
  NA_real_,
  EXPECTED_N_CELLS
)

PAPER_SIG_FULL <- rep(
  FALSE,
  EXPECTED_N_CELLS
)

P_FULL[KEEP] <- P_KEEP
SIG_FULL[KEEP] <- SIG_KEEP
PAPER_P_FULL[KEEP] <- PAPER_P_KEEP
PAPER_SIG_FULL[KEEP] <- PAPER_SIG_KEEP

# ==============================================================================
# 11. FINAL SENSITIVITY CELL TABLE
# ==============================================================================

SENS <- data.table(
  cell_id =
    rownames(MKNN),
  cell_type =
    PRIMARY$cell_type,
  donor =
    PRIMARY$donor,
  rhythm =
    PRIMARY$rhythm,
  AF_gchromVAR_Z_lowmatch_excluded =
    Z,
  AF_seed_lowmatch_excluded =
    SEED,
  network_score_lowmatch_excluded =
    NP,
  AF_TRS_lowmatch_excluded =
    TRS,
  network_eligible_lowmatch_excluded =
    KEEP,
  empirical_P_lowmatch_excluded =
    P_FULL,
  significant_lowmatch_excluded =
    SIG_FULL,
  paper_formula_P_lowmatch_excluded =
    PAPER_P_FULL,
  significant_paper_formula_lowmatch_excluded =
    PAPER_SIG_FULL
)

fwrite(
  SENS,
  FINAL_CELL_FILE,
  compress = "gzip"
)

# ==============================================================================
# 12. PRIMARY vs SENSITIVITY COMPARISON
# ==============================================================================

PRIMARY_SEED <- as.logical(
  PRIMARY$AF_SCAVENGE_seed
)

PRIMARY_SIG <- as.logical(
  PRIMARY$AF_SCAVENGE_significant_cell
)

Z_SPEARMAN <- suppressWarnings(
  cor(
    PRIMARY$AF_gchromVAR_Z,
    Z,
    method = "spearman",
    use = "complete.obs"
  )
)

TRS_SPEARMAN <- suppressWarnings(
  cor(
    PRIMARY$AF_TRS,
    TRS,
    method = "spearman",
    use = "complete.obs"
  )
)

SEED_JACCARD <- jaccard_bool(
  PRIMARY_SEED,
  SEED
)

SIG_JACCARD <- jaccard_bool(
  PRIMARY_SIG,
  SIG_FULL
)

PRIMARY_PAPER_AGREEMENT <- mean(
  SIG_FULL ==
    PAPER_SIG_FULL
)

CM_PRIMARY <- PRIMARY[
  cell_type == "CM"
]

CM_SENS <- SENS[
  cell_type == "CM"
]

PRIMARY_CM_SIG_FRAC <- mean(
  CM_PRIMARY$AF_SCAVENGE_significant_cell
)

SENS_CM_SIG_FRAC <- mean(
  CM_SENS$significant_lowmatch_excluded
)

COMPARISON <- data.table(
  metric = c(
    "excluded_low_reference_match_blocks",
    "gchromVAR_Z_Spearman_primary_vs_sensitivity",
    "seed_Jaccard_primary_vs_sensitivity",
    "TRS_Spearman_primary_vs_sensitivity",
    "significant_cell_Jaccard_primary_vs_sensitivity",
    "primary_n_seed",
    "sensitivity_n_seed",
    "primary_n_significant",
    "sensitivity_n_significant",
    "primary_CM_significant_fraction",
    "sensitivity_CM_significant_fraction",
    "absolute_CM_significant_fraction_change",
    "sensitivity_package_vs_paper_significance_agreement",
    "primary_network_eligible",
    "sensitivity_network_eligible",
    "sensitivity_RWR_iterations",
    "sensitivity_RWR_delta"
  ),
  value = c(
    length(LOW_LABELS),
    Z_SPEARMAN,
    SEED_JACCARD,
    TRS_SPEARMAN,
    SIG_JACCARD,
    sum(PRIMARY_SEED),
    sum(SEED),
    sum(PRIMARY_SIG),
    sum(SIG_FULL),
    PRIMARY_CM_SIG_FRAC,
    SENS_CM_SIG_FRAC,
    abs(
      PRIMARY_CM_SIG_FRAC -
        SENS_CM_SIG_FRAC
    ),
    PRIMARY_PAPER_AGREEMENT,
    sum(PRIMARY$network_eligible),
    sum(KEEP),
    REAL$iterations,
    REAL$delta
  )
)

fwrite(
  COMPARISON,
  COMPARISON_FILE
)

# ==============================================================================
# 13. CELL-TYPE + DONOR CONSISTENCY
# ==============================================================================

CELLTYPE <- SENS[
  ,
  .(
    n_cells = .N,
    n_seed =
      sum(
        AF_seed_lowmatch_excluded
      ),
    seed_fraction =
      mean(
        AF_seed_lowmatch_excluded
      ),
    n_significant =
      sum(
        significant_lowmatch_excluded
      ),
    significant_fraction =
      mean(
        significant_lowmatch_excluded
      ),
    mean_Z =
      mean(
        AF_gchromVAR_Z_lowmatch_excluded,
        na.rm = TRUE
      ),
    median_Z =
      median(
        AF_gchromVAR_Z_lowmatch_excluded,
        na.rm = TRUE
      ),
    mean_TRS =
      mean(
        AF_TRS_lowmatch_excluded,
        na.rm = TRUE
      ),
    median_TRS =
      median(
        AF_TRS_lowmatch_excluded,
        na.rm = TRUE
      )
  ),
  by = cell_type
]

setorder(
  CELLTYPE,
  -significant_fraction,
  -median_TRS
)

fwrite(
  CELLTYPE,
  CELLTYPE_FILE
)

DONOR_CT <- SENS[
  ,
  .(
    n_cells = .N,
    seed_fraction =
      mean(
        AF_seed_lowmatch_excluded
      ),
    significant_fraction =
      mean(
        significant_lowmatch_excluded
      ),
    median_Z =
      median(
        AF_gchromVAR_Z_lowmatch_excluded,
        na.rm = TRUE
      ),
    median_TRS =
      median(
        AF_TRS_lowmatch_excluded,
        na.rm = TRUE
      )
  ),
  by = .(
    donor,
    rhythm,
    cell_type
  )
]

DONOR_CT[
  ,
  rank_significant_fraction :=
    frank(
      -significant_fraction,
      ties.method = "min"
    ),
  by = donor
]

DONOR_CT[
  ,
  rank_median_TRS :=
    frank(
      -median_TRS,
      ties.method = "min"
    ),
  by = donor
]

fwrite(
  DONOR_CT,
  DONOR_CELLTYPE_FILE
)

CM_DONOR <- DONOR_CT[
  cell_type == "CM"
]

SAME_TOP_CELLTYPE_SIG <- (
  CELLTYPE$cell_type[1] ==
    "CM"
)

SAME_TOP_CELLTYPE_TRS <- (
  CELLTYPE[
    which.max(
      median_TRS
    ),
    cell_type
  ] ==
    "CM"
)

CM_RANK1_SIG_ALL_DONORS <- (
  nrow(CM_DONOR) == 7L &&
    all(
      CM_DONOR$rank_significant_fraction ==
        1L
    )
)

CM_RANK1_TRS_ALL_DONORS <- (
  nrow(CM_DONOR) == 7L &&
    all(
      CM_DONOR$rank_median_TRS ==
        1L
    )
)

# ==============================================================================
# 14. TECHNICAL READINESS
# ==============================================================================

READINESS <- data.table(
  check = c(
    "frozen_G3B_V3_primary_all_pass",
    "exactly_3_low_reference_match_blocks_excluded",
    "same_212084_peak_GC_object_reused",
    "same_200_background_iterations_reused",
    "same_frozen_V3_graph_reused",
    "sensitivity_gchromVAR_Z_finite_ge_99pct",
    "sensitivity_seed_cells_created",
    "sensitivity_RWR_converged",
    "sensitivity_network_eligible_ge_90pct",
    "1000_degree_matched_permutations_complete",
    "sensitivity_package_vs_paper_agreement_ge_99pct",
    "CM_remains_top_by_significant_fraction",
    "CM_remains_top_by_median_TRS",
    "CM_rank1_significant_fraction_all_7_donors",
    "CM_rank1_median_TRS_all_7_donors",
    "STEP11G4A_complete"
  ),
  pass = c(
    all(
      PRIMARY_READY$pass %in%
        c(TRUE, "TRUE", 1)
    ),
    length(LOW_LABELS) ==
      EXPECTED_LOW_BLOCKS,
    nrow(SE_GC) ==
      EXPECTED_N_PEAKS &&
      ncol(SE_GC) ==
        EXPECTED_N_CELLS,
    nrow(BG) ==
      EXPECTED_N_PEAKS &&
      ncol(BG) ==
        BACKGROUND_ITERATIONS,
    nrow(MKNN) ==
      EXPECTED_N_CELLS &&
      ncol(MKNN) ==
        EXPECTED_N_CELLS,
    mean(
      is.finite(Z)
    ) >= 0.99,
    sum(SEED) > 0,
    is.finite(
      REAL$delta
    ) &&
      REAL$delta <=
        RW_STATIONARY_CUTOFF,
    mean(KEEP) >=
      0.90,
    length(
      unique(PERM_SEEN)
    ) ==
      PERMUTATION_TIMES,
    PRIMARY_PAPER_AGREEMENT >=
      0.99,
    SAME_TOP_CELLTYPE_SIG,
    SAME_TOP_CELLTYPE_TRS,
    CM_RANK1_SIG_ALL_DONORS,
    CM_RANK1_TRS_ALL_DONORS,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G4A_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "sensitivity_question",
    "blocks_excluded_rule",
    "n_blocks_excluded",
    "PIP_threshold",
    "gchromVAR_GC_object",
    "gchromVAR_background",
    "network",
    "network_rebuilt",
    "RWR_gamma",
    "permutations",
    "primary_significance_rule",
    "paper_formula_sensitivity_rule"
  ),
  value = c(
    "Exclude fine-mapping blocks with reference_match_fraction <0.90",
    "Frozen STEP11G3A low-reference-match manifest",
    as.character(
      length(LOW_LABELS)
    ),
    as.character(PIP_MIN),
    "Same STEP11G3A GC-corrected peak-by-cell object",
    "Same STEP11G3A 200 matched-background iterations",
    "Same frozen STEP11G3B V3 official-source-equivalent mutual-kNN graph",
    "No",
    as.character(RW_GAMMA),
    as.character(PERMUTATION_TIMES),
    "count(null > real)/1000 <= 0.05",
    "(1 + count(null >= real))/(1+1000) < 0.05"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G4A_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G4A_sessionInfo.txt"
  )
)

# ==============================================================================
# 15. SIMPLE ROBUSTNESS FIGURES — ONE FIGURE PER FILE, NO LEGEND
# ==============================================================================

PLOT <- data.table(
  primary_TRS =
    PRIMARY$AF_TRS,
  sensitivity_TRS =
    TRS,
  cell_type =
    PRIMARY$cell_type
)

p1 <- ggplot2::ggplot(
  PLOT,
  ggplot2::aes(
    x = primary_TRS,
    y = sensitivity_TRS
  )
) +
  ggplot2::geom_point(
    size = 0.35,
    alpha = 0.25
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::labs(
    x = "Primary AF TRS",
    y = "TRS after excluding 3 low-match blocks"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  ) +
  ggplot2::theme(
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G4A_primary_vs_lowmatch_excluded_TRS.tiff"
  ),
  p1,
  width = 5.2,
  height = 4.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================================================================
# 16. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G4A COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Excluded blocks: ",
  length(LOW_LABELS),
  "\n",
  "gchromVAR Z Spearman: ",
  signif(Z_SPEARMAN, 5),
  "\n",
  "Seed Jaccard: ",
  signif(SEED_JACCARD, 5),
  "\n",
  "TRS Spearman: ",
  signif(TRS_SPEARMAN, 5),
  "\n",
  "Significant-cell Jaccard: ",
  signif(SIG_JACCARD, 5),
  "\n",
  "Primary CM significant fraction: ",
  sprintf("%.2f%%", 100 * PRIMARY_CM_SIG_FRAC),
  "\n",
  "Sensitivity CM significant fraction: ",
  sprintf("%.2f%%", 100 * SENS_CM_SIG_FRAC),
  "\n\n",
  sep = ""
)

cat(
  "Sensitivity cell-type summary:\n"
)
print(CELLTYPE)

cat(
  "\nCM donor consistency:\n"
)
print(CM_DONOR)

cat(
  "\nReadiness:\n"
)
print(READINESS)

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  FINAL_DIR,
  "\n",
  sep = ""
)
