# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G4C_V2_FIX_DATA_TABLE_ORDERING(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G4C V2 — CM INTERNAL MOTIF / TF REGULATORY ANALYSIS
# Project: AF–HFpEF–BMI–OSA shared genetics
# Dataset: GSE238242 left atrial appendage snATAC
#
# INPUT:
#   Frozen STEP11G4B CM High-vs-Low TRS pseudobulk differential accessibility
#   (39,016 tested peaks; 12 primary DARs; 11 robust DARs)
#
# WHY THIS ANALYSIS:
#   Direct overlap between robust CM-state DARs and AF fine-mapping PIP peaks
#   is zero. Therefore this step DOES NOT force a peak-level colocalization
#   narrative. Instead it asks whether regulatory motifs / TF programs converge:
#
#   A) CM state axis:
#      JASPAR2024 motif matches across ALL 39,016 tested peaks
#      -> rank-based motif enrichment using signed edgeR statistic
#      -> donor-direction consistency across 7 donors
#
#   B) AF genetic axis:
#      among the same 39,016 peaks, test whether motif-bearing peaks are
#      enriched among frozen AF PIP-linked peaks
#      -> logistic regression adjusted for GC, peak width, and mean accessibility
#
#   C) Integration:
#      identify motifs/TFs supported by both:
#        CM TRS-associated chromatin direction
#        + AF fine-mapping peak enrichment
#
#   D) Robust DAR network:
#      map prioritized motifs -> robust DARs -> nearest genes
#
# INTERPRETATION:
#   This is downstream regulatory characterization, not independent validation.
#   TRS and chromatin accessibility originate from the same snATAC dataset.
#
# OUTPUT ROOT:
#   D:/A/data/STEP11_GSE238242/04_SCAVENGE/
#   08_STEP11G4C_CM_MOTIF_TF/
# ==============================================================================

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = max(3600, getOption("timeout"))
)

set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

MOTIF_P_CUTOFF <- 5e-5

FGSEA_MIN_SIZE <- 20L
FGSEA_MAX_SIZE <- 10000L
STATE_FDR_CUTOFF <- 0.05
DONOR_CONSISTENCY_MIN <- 5L

GENETIC_FDR_CUTOFF <- 0.05
MIN_MOTIF_PEAKS_FOR_GENETIC_TEST <- 20L

TOP_N_STATE_PLOT <- 20L
TOP_N_GENETIC_PLOT <- 20L
TOP_N_NETWORK_MOTIFS <- 30L

# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

G4B_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "07_STEP11G4B_CM_TRS_CHROMATIN"
)

G4B_READINESS <- file.path(
  G4B_ROOT,
  "00_QC",
  "STEP11G4B_readiness.csv"
)

DA_ALL_FILE <- file.path(
  G4B_ROOT,
  "03_DIFFERENTIAL_ACCESSIBILITY",
  "STEP11G4B_CM_TRS_DA_all_tested_peaks.csv.gz"
)

DA_ROBUST_FILE <- file.path(
  G4B_ROOT,
  "03_DIFFERENTIAL_ACCESSIBILITY",
  "STEP11G4B_CM_TRS_robust_DAR_5of7donors.csv.gz"
)

DA_SIG_FILE <- file.path(
  G4B_ROOT,
  "03_DIFFERENTIAL_ACCESSIBILITY",
  "STEP11G4B_CM_TRS_DAR_FDR05_logFC05.csv.gz"
)

G4B_DA_QC <- file.path(
  G4B_ROOT,
  "00_QC",
  "STEP11G4B_differential_accessibility_QC.csv"
)

OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "08_STEP11G4C_CM_MOTIF_TF"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

MOTIF_DIR <- file.path(
  OUT_ROOT,
  "01_MOTIF_MATCH"
)

STATE_DIR <- file.path(
  OUT_ROOT,
  "02_STATE_MOTIF_ENRICHMENT"
)

GENETIC_DIR <- file.path(
  OUT_ROOT,
  "03_AF_GENETIC_MOTIF_ENRICHMENT"
)

INTEGRATION_DIR <- file.path(
  OUT_ROOT,
  "04_INTEGRATION"
)

NETWORK_DIR <- file.path(
  OUT_ROOT,
  "05_TF_PEAK_GENE_NETWORK"
)

FIG_DIR <- file.path(
  OUT_ROOT,
  "06_FIGURES"
)

for (d in c(
  OUT_ROOT,
  QC_DIR,
  MOTIF_DIR,
  STATE_DIR,
  GENETIC_DIR,
  INTEGRATION_DIR,
  NETWORK_DIR,
  FIG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

MOTIF_MATCH_FILE <- file.path(
  MOTIF_DIR,
  "STEP11G4C_JASPAR2024_motif_matches_39016peaks.rds"
)

MOTIF_META_FILE <- file.path(
  MOTIF_DIR,
  "STEP11G4C_JASPAR2024_motif_metadata.csv"
)

STATE_RESULT_FILE <- file.path(
  STATE_DIR,
  "STEP11G4C_CM_state_ranked_motif_enrichment.csv"
)

STATE_ROBUST_FILE <- file.path(
  STATE_DIR,
  "STEP11G4C_CM_state_robust_motifs.csv"
)

STATE_DONOR_FILE <- file.path(
  STATE_DIR,
  "STEP11G4C_CM_state_motif_donor_effects.csv"
)

GENETIC_RESULT_FILE <- file.path(
  GENETIC_DIR,
  "STEP11G4C_AF_PIP_motif_logistic_enrichment.csv"
)

INTEGRATED_FILE <- file.path(
  INTEGRATION_DIR,
  "STEP11G4C_integrated_TF_motif_prioritization.csv"
)

TF_SUMMARY_FILE <- file.path(
  INTEGRATION_DIR,
  "STEP11G4C_integrated_TF_summary.csv"
)

NETWORK_FILE <- file.path(
  NETWORK_DIR,
  "STEP11G4C_TF_robustDAR_nearestGene_edges.csv"
)

READINESS_FILE <- file.path(
  QC_DIR,
  "STEP11G4C_readiness.csv"
)

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

cran_pkgs <- c(
  "data.table",
  "Matrix",
  "RSQLite",
  "ggplot2"
)

for (p in cran_pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

bioc_pkgs <- c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg38",
  "TFBSTools",
  "JASPAR2024",
  "motifmatchr",
  "SummarizedExperiment",
  "fgsea"
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

still_missing <- c(
  cran_pkgs,
  bioc_pkgs
)[
  !vapply(
    c(
      cran_pkgs,
      bioc_pkgs
    ),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(still_missing)) {
  stop(
    paste0(
      "Required package(s) could not be loaded:\n",
      paste(
        still_missing,
        collapse = "\n"
      )
    )
  )
}

library(data.table)
library(Matrix)
library(ggplot2)

HG38 <-
  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# ==============================================================================
# 3. HARD INPUT CHECK
# ==============================================================================

required_inputs <- c(
  G4B_READINESS,
  DA_ALL_FILE,
  DA_ROBUST_FILE,
  DA_SIG_FILE,
  G4B_DA_QC
)

if (!all(file.exists(required_inputs))) {

  stop(
    paste0(
      "Missing STEP11G4C input(s):\n",
      paste(
        required_inputs[
          !file.exists(
            required_inputs
          )
        ],
        collapse = "\n"
      )
    )
  )
}

READY_B <- fread(
  G4B_READINESS
)

if (
  nrow(READY_B) != 15L ||
  !all(
    READY_B$pass %in%
      c(
        TRUE,
        "TRUE",
        1
      )
  )
) {

  print(
    READY_B
  )

  stop(
    "STEP11G4B readiness is not fully PASS."
  )
}

DA <- fread(
  DA_ALL_FILE,
  showProgress = TRUE
)

ROBUST <- fread(
  DA_ROBUST_FILE
)

PRIMARY_DAR <- fread(
  DA_SIG_FILE
)

required_da_cols <- c(
  "peak_id",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR",
  "chr",
  "start",
  "end",
  "AF_PIP_peak_weight",
  "AF_PIP_linked"
)

missing_da_cols <- setdiff(
  required_da_cols,
  names(DA)
)

if (length(missing_da_cols)) {

  stop(
    paste0(
      "DA table missing required columns: ",
      paste(
        missing_da_cols,
        collapse = ", "
      )
    )
  )
}

DELTA_COLS <- grep(
  "^delta_logCPM_",
  names(DA),
  value = TRUE
)

if (length(DELTA_COLS) != 7L) {

  stop(
    "Expected 7 donor delta_logCPM columns; observed ",
    length(DELTA_COLS)
  )
}

DONORS <- sub(
  "^delta_logCPM_",
  "",
  DELTA_COLS
)

if (nrow(DA) < 1000L) {
  stop(
    "Too few tested peaks for motif analysis."
  )
}

if (anyDuplicated(
  DA$peak_id
)) {
  stop(
    "DA table peak_id is not unique."
  )
}

# Preserve the original frozen G4B peak order.
# No sorting is required here because downstream motif matrix rows must remain
# exactly aligned to the all-tested DA table.

cat(
  "\nFrozen G4B input:\n",
  "  tested peaks: ",
  format(
    nrow(DA),
    big.mark = ","
  ),
  "\n",
  "  primary DARs: ",
  nrow(PRIMARY_DAR),
  "\n",
  "  robust DARs: ",
  nrow(ROBUST),
  "\n",
  "  AF PIP-linked tested peaks: ",
  sum(
    DA$AF_PIP_peak_weight >
      0
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 4. TESTED PEAK GRANGES + SEQUENCE COVARIATES
# ==============================================================================

PEAK_GR <- GenomicRanges::GRanges(
  seqnames =
    as.character(
      DA$chr
    ),
  ranges = IRanges::IRanges(
    start =
      as.integer(
        DA$start
      ),
    end =
      as.integer(
        DA$end
      )
  )
)

names(
  PEAK_GR
) <- DA$peak_id

# Keep only standard autosomal/tested chromosomes with valid hg38 sequence.
VALID_CHR <- as.character(
  GenomeInfoDb::seqnames(
    HG38
  )
)

VALID <- as.character(
  GenomicRanges::seqnames(
    PEAK_GR
  )
) %in%
  VALID_CHR

if (!all(VALID)) {

  stop(
    "At least one tested peak chromosome is absent from hg38 BSgenome."
  )
}

SEQ <- BSgenome::getSeq(
  HG38,
  PEAK_GR
)

GC <- as.numeric(
  Biostrings::letterFrequency(
    SEQ,
    letters = c(
      "G",
      "C"
    ),
    as.prob = TRUE
  )[
    ,
    1
  ] +
    Biostrings::letterFrequency(
      SEQ,
      letters = c(
        "G",
        "C"
      ),
      as.prob = TRUE
    )[
      ,
      2
    ]
)

PEAK_WIDTH <- width(
  PEAK_GR
)

COVAR <- data.table(
  peak_id =
    DA$peak_id,
  GC =
    GC,
  peak_width =
    PEAK_WIDTH,
  logCPM =
    DA$logCPM
)

fwrite(
  COVAR,
  file.path(
    QC_DIR,
    "STEP11G4C_peak_covariates.csv.gz"
  ),
  compress = "gzip"
)

rm(
  SEQ
)

invisible(
  gc(
    verbose = FALSE
  )
)

# ==============================================================================
# 5. LOAD JASPAR2024 HUMAN CORE MOTIFS
# ==============================================================================

cat(
  "\nLoading JASPAR2024 Homo sapiens CORE motifs...\n"
)

J24 <- JASPAR2024::JASPAR2024()

J24_CONN <- RSQLite::dbConnect(
  RSQLite::SQLite(),
  JASPAR2024::db(
    J24
  )
)

on.exit(
  try(
    RSQLite::dbDisconnect(
      J24_CONN
    ),
    silent = TRUE
  ),
  add = TRUE
)

PFM <- TFBSTools::getMatrixSet(
  x = J24_CONN,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    species = "Homo sapiens",
    all_versions = FALSE,
    matrixtype = "PFM"
  )
)

if (length(PFM) < 100L) {
  stop(
    "Unexpectedly few JASPAR2024 human CORE motifs: ",
    length(PFM)
  )
}

MOTIF_ID <- names(
  PFM
)

MOTIF_NAME <- vapply(
  PFM,
  function(x) {
    as.character(
      TFBSTools::name(
        x
      )
    )
  },
  character(1)
)

MOTIF_CLASS <- vapply(
  PFM,
  function(x) {

    z <- tryCatch(
      as.character(
        TFBSTools::tags(
          x
        )$class
      ),
      error = function(e) NA_character_
    )

    if (!length(z)) {
      NA_character_
    } else {
      paste(
        z,
        collapse = ";"
      )
    }
  },
  character(1)
)

MOTIF_META <- data.table(
  motif_id =
    MOTIF_ID,
  motif_name =
    MOTIF_NAME,
  motif_class =
    MOTIF_CLASS
)

fwrite(
  MOTIF_META,
  MOTIF_META_FILE
)

cat(
  "JASPAR motifs loaded: ",
  nrow(
    MOTIF_META
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 6. MOTIF MATCHING — CHECKPOINT
# ==============================================================================

MOTIF_MATCH <- NULL

if (file.exists(
  MOTIF_MATCH_FILE
)) {

  obj <- tryCatch(
    readRDS(
      MOTIF_MATCH_FILE
    ),
    error = function(e) NULL
  )

  if (
    !is.null(obj) &&
    inherits(
      obj,
      "sparseMatrix"
    ) &&
    nrow(obj) ==
      nrow(DA) &&
    ncol(obj) ==
      nrow(MOTIF_META) &&
    identical(
      rownames(obj),
      DA$peak_id
    ) &&
    identical(
      colnames(obj),
      MOTIF_META$motif_id
    )
  ) {
    MOTIF_MATCH <- obj
  }
}

if (is.null(
  MOTIF_MATCH
)) {

  cat(
    "\nMatching JASPAR2024 motifs to ",
    format(
      nrow(DA),
      big.mark = ","
    ),
    " hg38 tested peaks...\n",
    sep = ""
  )

  MATCH_SE <- motifmatchr::matchMotifs(
    pwms = PFM,
    subject = PEAK_GR,
    genome = HG38,
    out = "matches",
    p.cutoff = MOTIF_P_CUTOFF
  )

  MOTIF_MATCH <- motifmatchr::motifMatches(
    MATCH_SE
  )

  MOTIF_MATCH <- methods::as(
    MOTIF_MATCH,
    "lgCMatrix"
  )

  rownames(
    MOTIF_MATCH
  ) <- DA$peak_id

  colnames(
    MOTIF_MATCH
  ) <- MOTIF_META$motif_id

  saveRDS(
    MOTIF_MATCH,
    MOTIF_MATCH_FILE,
    compress = TRUE
  )
}

N_MOTIF_PEAKS <- as.integer(
  Matrix::colSums(
    MOTIF_MATCH
  )
)

MOTIF_META[
  ,
  n_tested_peaks_with_motif :=
    N_MOTIF_PEAKS
]

fwrite(
  MOTIF_META,
  MOTIF_META_FILE
)

MOTIF_MATCH_QC <- data.table(
  metric = c(
    "tested_peaks",
    "JASPAR2024_human_CORE_motifs",
    "motifs_with_at_least_20_peaks",
    "median_peaks_per_motif",
    "min_peaks_per_motif",
    "max_peaks_per_motif",
    "motif_match_p_cutoff"
  ),
  value = c(
    nrow(DA),
    ncol(MOTIF_MATCH),
    sum(
      N_MOTIF_PEAKS >=
        FGSEA_MIN_SIZE
    ),
    median(
      N_MOTIF_PEAKS
    ),
    min(
      N_MOTIF_PEAKS
    ),
    max(
      N_MOTIF_PEAKS
    ),
    MOTIF_P_CUTOFF
  )
)

fwrite(
  MOTIF_MATCH_QC,
  file.path(
    QC_DIR,
    "STEP11G4C_motif_match_QC.csv"
  )
)

# ==============================================================================
# 7. CM STATE MOTIF ENRICHMENT — FULL RANKED PEAK SET
# ==============================================================================

cat(
  "\nRunning ranked motif enrichment across all tested peaks...\n"
)

# For a one-degree-of-freedom edgeR QL contrast, sign(logFC)*sqrt(F)
# behaves as a signed test statistic and retains effect direction.
SIGNED_STAT <- sign(
  DA$logFC
) *
  sqrt(
    pmax(
      DA$F,
      0
    )
  )

names(
  SIGNED_STAT
) <- DA$peak_id

SIGNED_STAT <- sort(
  SIGNED_STAT,
  decreasing = TRUE
)

VALID_MOTIFS <- which(
  N_MOTIF_PEAKS >=
    FGSEA_MIN_SIZE &
    N_MOTIF_PEAKS <=
      FGSEA_MAX_SIZE
)

PATHWAYS <- lapply(
  VALID_MOTIFS,
  function(j) {

    DA$peak_id[
      as.logical(
        MOTIF_MATCH[
          ,
          j
        ]
      )
    ]
  }
)

names(
  PATHWAYS
) <- MOTIF_META$motif_id[
  VALID_MOTIFS
]

FG <- fgsea::fgseaMultilevel(
  pathways = PATHWAYS,
  stats = SIGNED_STAT,
  minSize = FGSEA_MIN_SIZE,
  maxSize = FGSEA_MAX_SIZE,
  eps = 0
)

FG <- as.data.table(
  FG
)

setnames(
  FG,
  "pathway",
  "motif_id"
)

# leadingEdge is a list-column and should be exported separately/flattened.
FG[
  ,
  leading_edge_peaks :=
    vapply(
      leadingEdge,
      function(x) {
        paste(
          x,
          collapse = ";"
        )
      },
      character(1)
    )
]

FG[
  ,
  leadingEdge :=
    NULL
]

FG <- merge(
  FG,
  MOTIF_META,
  by = "motif_id",
  all.x = TRUE,
  sort = FALSE
)

# ==============================================================================
# 8. DONOR DIRECTION CONSISTENCY FOR EACH MOTIF
# ==============================================================================

DELTA <- as.matrix(
  DA[
    ,
    ..DELTA_COLS
  ]
)

storage.mode(
  DELTA
) <- "double"

# Convert logical sparse motif incidence to numeric 0/1 sparse matrix.
M_NUM <- methods::as(
  MOTIF_MATCH,
  "dgCMatrix"
)

M_NUM@x[] <- 1

MOTIF_COUNTS <- Matrix::colSums(
  M_NUM
)

# motif x donor mean peak delta-logCPM
MOTIF_DONOR_MEAN <- as.matrix(
  Matrix::crossprod(
    M_NUM,
    DELTA
  )
)

MOTIF_DONOR_MEAN <- sweep(
  MOTIF_DONOR_MEAN,
  1,
  MOTIF_COUNTS,
  "/"
)

colnames(
  MOTIF_DONOR_MEAN
) <- DONORS

rownames(
  MOTIF_DONOR_MEAN
) <- colnames(
  M_NUM
)

DONOR_LONG <- as.data.table(
  as.table(
    MOTIF_DONOR_MEAN
  )
)

setnames(
  DONOR_LONG,
  c(
    "motif_id",
    "donor",
    "mean_delta_logCPM"
  )
)

DONOR_LONG <- merge(
  DONOR_LONG,
  MOTIF_META[
    ,
    .(
      motif_id,
      motif_name
    )
  ],
  by = "motif_id",
  all.x = TRUE
)

fwrite(
  DONOR_LONG,
  STATE_DONOR_FILE
)

FG[
  ,
  n_donors_positive :=
    vapply(
      motif_id,
      function(m) {
        sum(
          MOTIF_DONOR_MEAN[
            m,
          ] >
            0,
          na.rm = TRUE
        )
      },
      integer(1)
    )
]

FG[
  ,
  n_donors_negative :=
    vapply(
      motif_id,
      function(m) {
        sum(
          MOTIF_DONOR_MEAN[
            m,
          ] <
            0,
          na.rm = TRUE
        )
      },
      integer(1)
    )
]

FG[
  ,
  state_direction :=
    fifelse(
      NES >
        0,
      "HighTRS_open",
      "LowTRS_open"
    )
]

FG[
  ,
  n_donors_concordant :=
    fifelse(
      NES >
        0,
      n_donors_positive,
      n_donors_negative
    )
]

FG[
  ,
  donor_consistency_fraction :=
    n_donors_concordant /
    7
]

FG[
  ,
  state_robust :=
    padj <
      STATE_FDR_CUTOFF &
    n_donors_concordant >=
      DONOR_CONSISTENCY_MIN
]

FG[, abs_NES_for_order := abs(NES)]

setorderv(
  FG,
  cols = c("padj", "abs_NES_for_order"),
  order = c(1L, -1L),
  na.last = TRUE
)

FG[, abs_NES_for_order := NULL]

fwrite(
  FG,
  STATE_RESULT_FILE
)

STATE_ROBUST <- FG[
  state_robust ==
    TRUE
]

fwrite(
  STATE_ROBUST,
  STATE_ROBUST_FILE
)

# ==============================================================================
# 9. AF PIP-LINKED MOTIF ENRICHMENT
#    Logistic model adjusts for GC, peak width, and mean accessibility.
# ==============================================================================

cat(
  "\nTesting motif enrichment among AF PIP-linked peaks...\n"
)

Y <- as.integer(
  DA$AF_PIP_peak_weight >
    0
)

if (sum(Y) < 20L) {
  stop(
    "Too few AF PIP-linked tested peaks for motif enrichment."
  )
}

Z_GC <- as.numeric(
  scale(
    GC
  )
)

Z_WIDTH <- as.numeric(
  scale(
    log1p(
      PEAK_WIDTH
    )
  )
)

Z_LOGCPM <- as.numeric(
  scale(
    DA$logCPM
  )
)

GEN_LIST <- vector(
  "list",
  ncol(
    MOTIF_MATCH
  )
)

for (j in seq_len(
  ncol(
    MOTIF_MATCH
  )
)) {

  xmotif <- as.numeric(
    MOTIF_MATCH[
      ,
      j
    ]
  )

  n_pos <- sum(
    xmotif >
      0
  )

  n_pip_with_motif <- sum(
    Y ==
      1 &
      xmotif >
        0
  )

  n_nonpip_with_motif <- sum(
    Y ==
      0 &
      xmotif >
        0
  )

  if (
    n_pos <
      MIN_MOTIF_PEAKS_FOR_GENETIC_TEST ||
    n_pos >
      (
        length(
          Y
        ) -
          MIN_MOTIF_PEAKS_FOR_GENETIC_TEST
      )
  ) {

    GEN_LIST[[
      j
    ]] <- data.table(
      motif_id =
        colnames(
          MOTIF_MATCH
        )[
          j
        ],
      n_motif_peaks =
        n_pos,
      n_PIP_linked_with_motif =
        n_pip_with_motif,
      beta_motif =
        NA_real_,
      SE_motif =
        NA_real_,
      OR_motif =
        NA_real_,
      P_motif =
        NA_real_,
      model_converged =
        FALSE
    )

    next
  }

  X <- cbind(
    intercept = 1,
    motif = xmotif,
    GC = Z_GC,
    width = Z_WIDTH,
    logCPM = Z_LOGCPM
  )

  FIT <- tryCatch(
    glm.fit(
      x = X,
      y = Y,
      family = binomial()
    ),
    error = function(e) NULL
  )

  if (
    is.null(FIT) ||
    !isTRUE(
      FIT$converged
    ) ||
    FIT$rank <
      ncol(
        X
      )
  ) {

    GEN_LIST[[
      j
    ]] <- data.table(
      motif_id =
        colnames(
          MOTIF_MATCH
        )[
          j
        ],
      n_motif_peaks =
        n_pos,
      n_PIP_linked_with_motif =
        n_pip_with_motif,
      beta_motif =
        NA_real_,
      SE_motif =
        NA_real_,
      OR_motif =
        NA_real_,
      P_motif =
        NA_real_,
      model_converged =
        FALSE
    )

    next
  }

  # glm.fit can be converted to class glm for standard coefficient SEs.
  class(
    FIT
  ) <- c(
    "glm",
    "lm"
  )

  SM <- tryCatch(
    summary(
      FIT
    ),
    error = function(e) NULL
  )

  if (
    is.null(SM) ||
    !"motif" %in%
      rownames(
        SM$coefficients
      )
  ) {

    beta <- se <- pval <- NA_real_

  } else {

    beta <- SM$coefficients[
      "motif",
      "Estimate"
    ]

    se <- SM$coefficients[
      "motif",
      "Std. Error"
    ]

    pval <- SM$coefficients[
      "motif",
      "Pr(>|z|)"
    ]
  }

  GEN_LIST[[
    j
  ]] <- data.table(
    motif_id =
      colnames(
        MOTIF_MATCH
      )[
        j
      ],
    n_motif_peaks =
      n_pos,
    n_PIP_linked_with_motif =
      n_pip_with_motif,
    beta_motif =
      beta,
    SE_motif =
      se,
    OR_motif =
      exp(
        beta
      ),
    P_motif =
      pval,
    model_converged =
      isTRUE(
        FIT$converged
      )
  )
}

GEN <- rbindlist(
  GEN_LIST,
  fill = TRUE
)

GEN[
  ,
  FDR_BH :=
    p.adjust(
      P_motif,
      method = "BH"
    )
]

GEN[
  ,
  genetic_enriched :=
    is.finite(
      FDR_BH
    ) &
    FDR_BH <
      GENETIC_FDR_CUTOFF &
    is.finite(
      beta_motif
    ) &
    beta_motif >
      0
]

GEN <- merge(
  GEN,
  MOTIF_META,
  by = "motif_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  GEN,
  FDR_BH,
  -OR_motif
)

fwrite(
  GEN,
  GENETIC_RESULT_FILE
)

# ==============================================================================
# 10. INTEGRATE STATE + GENETIC MOTIF EVIDENCE
# ==============================================================================

INTEGRATED <- merge(
  FG[
    ,
    .(
      motif_id,
      motif_name,
      motif_class,
      size,
      NES,
      pval,
      padj,
      state_direction,
      n_donors_positive,
      n_donors_negative,
      n_donors_concordant,
      donor_consistency_fraction,
      state_robust,
      leading_edge_peaks
    )
  ],
  GEN[
    ,
    .(
      motif_id,
      n_motif_peaks,
      n_PIP_linked_with_motif,
      beta_motif,
      SE_motif,
      OR_motif,
      P_motif,
      FDR_BH,
      genetic_enriched
    )
  ],
  by = "motif_id",
  all = TRUE,
  sort = FALSE
)

# Recover metadata where motif fell outside the fgsea size range.
INTEGRATED <- merge(
  MOTIF_META[
    ,
    .(
      motif_id,
      motif_name_meta =
        motif_name,
      motif_class_meta =
        motif_class
    )
  ],
  INTEGRATED,
  by = "motif_id",
  all.y = TRUE,
  sort = FALSE
)

INTEGRATED[
  ,
  motif_name :=
    fifelse(
      is.na(
        motif_name
      ) |
        motif_name ==
          "",
      motif_name_meta,
      motif_name
    )
]

INTEGRATED[
  ,
  motif_class :=
    fifelse(
      is.na(
        motif_class
      ) |
        motif_class ==
          "",
      motif_class_meta,
      motif_class
    )
]

INTEGRATED[
  ,
  c(
    "motif_name_meta",
    "motif_class_meta"
  ) :=
    NULL
]

INTEGRATED[
  ,
  convergent_state_and_genetic :=
    state_robust %
      in% TRUE &
    genetic_enriched %
      in% TRUE
]

# Ranking score is descriptive and is NOT a formal probability.
INTEGRATED[
  ,
  integration_score :=
    fifelse(
      is.finite(
        padj
      ),
      -log10(
        pmax(
          padj,
          1e-300
        )
      ),
      0
    ) +
    fifelse(
      is.finite(
        FDR_BH
      ),
      -log10(
        pmax(
          FDR_BH,
          1e-300
        )
      ),
      0
    ) +
    fifelse(
      is.finite(
        donor_consistency_fraction
      ),
      donor_consistency_fraction,
      0
    )
]

setorder(
  INTEGRATED,
  -convergent_state_and_genetic,
  -integration_score,
  padj
)

fwrite(
  INTEGRATED,
  INTEGRATED_FILE
)

# TF-level collapse: JASPAR can contain more than one motif per TF.
TF_SUMMARY <- INTEGRATED[
  !is.na(
    motif_name
  ) &
    motif_name !=
      "",
  {

    state_idx <- which.min(
      fifelse(
        is.finite(
          padj
        ),
        padj,
        Inf
      )
    )

    genetic_idx <- which.min(
      fifelse(
        is.finite(
          FDR_BH
        ),
        FDR_BH,
        Inf
      )
    )

    .(
      n_motifs =
        .N,
      best_state_motif_id =
        motif_id[
          state_idx
        ],
      best_state_NES =
        NES[
          state_idx
        ],
      best_state_FDR =
        padj[
          state_idx
        ],
      best_state_direction =
        state_direction[
          state_idx
        ],
      best_state_donor_concordance =
        n_donors_concordant[
          state_idx
        ],
      any_state_robust =
        any(
          state_robust %
            in% TRUE,
          na.rm = TRUE
        ),
      best_genetic_motif_id =
        motif_id[
          genetic_idx
        ],
      best_genetic_OR =
        OR_motif[
          genetic_idx
        ],
      best_genetic_FDR =
        FDR_BH[
          genetic_idx
        ],
      any_genetic_enriched =
        any(
          genetic_enriched %
            in% TRUE,
          na.rm = TRUE
        ),
      any_convergent =
        any(
          convergent_state_and_genetic %
            in% TRUE,
          na.rm = TRUE
        ),
      max_integration_score =
        max(
          integration_score,
          na.rm = TRUE
        )
    )
  },
  by = motif_name
]

setorder(
  TF_SUMMARY,
  -any_convergent,
  -max_integration_score,
  best_state_FDR
)

fwrite(
  TF_SUMMARY,
  TF_SUMMARY_FILE
)

# ==============================================================================
# 11. PRIORITIZED TF -> ROBUST DAR -> NEAREST GENE NETWORK
# ==============================================================================

TOP_MOTIFS <- unique(
  c(
    INTEGRATED[
      convergent_state_and_genetic %
        in% TRUE,
      motif_id
    ],
    INTEGRATED[
      state_robust %
        in% TRUE
    ][
      order(
        padj,
        -abs(
          NES
        )
      ),
      head(
        motif_id,
        TOP_N_NETWORK_MOTIFS
      )
    ]
  )
)

ROBUST_IDX <- match(
  ROBUST$peak_id,
  DA$peak_id
)

if (anyNA(
  ROBUST_IDX
)) {
  stop(
    "At least one robust DAR is absent from all-tested peak table."
  )
}

NETWORK <- data.table()

if (
  length(
    TOP_MOTIFS
  ) >
    0 &&
  nrow(
    ROBUST
  ) >
    0
) {

  motif_col_idx <- match(
    TOP_MOTIFS,
    colnames(
      MOTIF_MATCH
    )
  )

  motif_col_idx <- motif_col_idx[
    !is.na(
      motif_col_idx
    )
  ]

  NET_LIST <- list()

  for (jj in motif_col_idx) {

    present <- as.logical(
      MOTIF_MATCH[
        ROBUST_IDX,
        jj
      ]
    )

    if (!any(
      present
    )) {
      next
    }

    tmp <- ROBUST[
      present
    ]

    tmp[
      ,
      `:=`(
        motif_id =
          colnames(
            MOTIF_MATCH
          )[
            jj
          ],
        motif_name =
          MOTIF_META[
            motif_id ==
              colnames(
                MOTIF_MATCH
              )[
                jj
              ],
            motif_name
          ]
      )
    ]

    NET_LIST[[
      length(
        NET_LIST
      ) +
        1L
    ]] <- tmp[
      ,
      .(
        motif_id,
        motif_name,
        peak_id,
        peak_direction =
          direction,
        peak_logFC =
          logFC,
        peak_FDR =
          FDR,
        donor_consistency_fraction,
        AF_PIP_peak_weight,
        AF_PIP_linked,
        nearest_gene_symbol,
        nearest_gene_name,
        distance_to_nearest_TSS
      )
    ]
  }

  if (length(
    NET_LIST
  )) {

    NETWORK <- unique(
      rbindlist(
        NET_LIST,
        fill = TRUE
      )
    )
  }
}

fwrite(
  NETWORK,
  NETWORK_FILE
)

# ==============================================================================
# 12. QC / READINESS
# ==============================================================================

N_STATE_ROBUST <- nrow(
  STATE_ROBUST
)

N_GENETIC_ENRICHED <- sum(
  GEN$genetic_enriched,
  na.rm = TRUE
)

N_CONVERGENT <- sum(
  INTEGRATED$convergent_state_and_genetic,
  na.rm = TRUE
)

QC_SUMMARY <- data.table(
  metric = c(
    "tested_peaks",
    "primary_DARs",
    "robust_DARs",
    "robust_PIP_linked_DARs",
    "AF_PIP_linked_tested_peaks",
    "JASPAR2024_human_CORE_motifs",
    "motifs_tested_state_fgsea",
    "state_robust_motifs_FDR05_5of7",
    "motifs_tested_AF_genetic_logistic",
    "AF_genetic_enriched_motifs_FDR05",
    "convergent_state_and_genetic_motifs",
    "TF_peak_gene_network_edges"
  ),
  value = c(
    nrow(DA),
    nrow(PRIMARY_DAR),
    nrow(ROBUST),
    sum(
      ROBUST$AF_PIP_linked,
      na.rm = TRUE
    ),
    sum(
      DA$AF_PIP_peak_weight >
        0
    ),
    nrow(MOTIF_META),
    nrow(FG),
    N_STATE_ROBUST,
    sum(
      GEN$model_converged,
      na.rm = TRUE
    ),
    N_GENETIC_ENRICHED,
    N_CONVERGENT,
    nrow(NETWORK)
  )
)

fwrite(
  QC_SUMMARY,
  file.path(
    QC_DIR,
    "STEP11G4C_analysis_QC.csv"
  )
)

READINESS <- data.table(
  check = c(
    "STEP11G4B_all_15_checks_pass",
    "39016_tested_peaks_or_more_than_1000",
    "JASPAR2024_human_CORE_motifs_loaded",
    "motif_match_matrix_created",
    "motif_match_rows_align_with_DA_peaks",
    "state_ranked_motif_enrichment_completed",
    "7_donor_motif_direction_summary_completed",
    "AF_PIP_genetic_motif_models_completed",
    "integrated_motif_table_created",
    "TF_level_summary_created",
    "robust_DAR_network_exported",
    "STEP11G4C_complete"
  ),
  pass = c(
    nrow(READY_B) ==
      15L &&
      all(
        READY_B$pass %in%
          c(
            TRUE,
            "TRUE",
            1
          )
      ),
    nrow(DA) >
      1000L,
    nrow(MOTIF_META) >
      100L,
    inherits(
      MOTIF_MATCH,
      "sparseMatrix"
    ),
    identical(
      rownames(
        MOTIF_MATCH
      ),
      DA$peak_id
    ),
    nrow(FG) >
      0L,
    nrow(DONOR_LONG) ==
      ncol(
        MOTIF_MATCH
      ) *
        7L,
    sum(
      GEN$model_converged,
      na.rm = TRUE
    ) >
      100L,
    file.exists(
      INTEGRATED_FILE
    ),
    file.exists(
      TF_SUMMARY_FILE
    ),
    file.exists(
      NETWORK_FILE
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
    "analysis_scope",
    "motif_database",
    "motif_subset",
    "motif_matcher",
    "motif_match_p_cutoff",
    "state_motif_test",
    "state_rank_statistic",
    "state_robust_rule",
    "genetic_motif_test",
    "genetic_covariates",
    "genetic_positive_definition",
    "integration_rule",
    "robust_DAR_PIP_overlap",
    "interpretation_caution"
  ),
  value = c(
    "CM High-vs-Low AF-TRS downstream regulatory characterization",
    "JASPAR2024",
    "Homo sapiens CORE; vertebrates; latest versions only",
    "motifmatchr::matchMotifs",
    as.character(
      MOTIF_P_CUTOFF
    ),
    "fgseaMultilevel across all G4B-tested peaks",
    "sign(logFC)*sqrt(edgeR F)",
    paste0(
      "BH FDR<",
      STATE_FDR_CUTOFF,
      " plus motif-average direction concordant in >=",
      DONOR_CONSISTENCY_MIN,
      "/7 donors"
    ),
    "Peak-level binomial logistic regression",
    "motif + GC + log(peak width) + mean logCPM",
    "Frozen AF_PIP_peak_weight > 0 among the same G4B-tested peaks",
    "State-robust motif AND AF-genetic-enriched motif",
    as.character(
      sum(
        ROBUST$AF_PIP_linked,
        na.rm = TRUE
      )
    ),
    paste0(
      "Motif-level convergence is not peak-level colocalization; ",
      "TRS/state accessibility derive from the same snATAC dataset."
    )
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G4C_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G4C_sessionInfo.txt"
  )
)

# ==============================================================================
# 13. FIGURES — SINGLE FIGURE PER FILE, NO LEGEND
# ==============================================================================

# ------------------------------------------------------------------
# Figure 1: top state-associated motifs
# ------------------------------------------------------------------

P_STATE <- FG[
  is.finite(
    NES
  ) &
    is.finite(
      padj
    )
]

P_STATE <- rbind(
  head(
    P_STATE[
      NES >
        0
    ][
      order(
        padj,
        -NES
      )
    ],
    ceiling(
      TOP_N_STATE_PLOT /
        2
    )
  ),
  head(
    P_STATE[
      NES <
        0
    ][
      order(
        padj,
        NES
      )
    ],
    floor(
      TOP_N_STATE_PLOT /
        2
    )
  ),
  fill = TRUE
)

if (nrow(
  P_STATE
)) {

  P_STATE[
    ,
    label :=
      paste0(
        motif_name,
        " [",
        motif_id,
        "]"
      )
  ]

  P_STATE[
    ,
    label :=
      factor(
        label,
        levels =
          label[
            order(
              NES
            )
          ]
      )
  ]

  p1 <- ggplot(
    P_STATE,
    aes(
      x = NES,
      y = label,
      fill = state_direction
    )
  ) +
    geom_col(
      width = 0.72,
      linewidth = 0.35,
      colour = "black"
    ) +
    scale_fill_manual(
      values = c(
        HighTRS_open = "#B34D4D",
        LowTRS_open = "#557DA1"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.45
    ) +
    labs(
      x = "Rank-based motif enrichment NES",
      y = NULL
    ) +
    theme_classic(
      base_size = 10.5
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(
        colour = "black"
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_CM_state_top_motif_NES.tiff"
    ),
    p1,
    width = 6.3,
    height = 6.0,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_CM_state_top_motif_NES.pdf"
    ),
    p1,
    width = 6.3,
    height = 6.0,
    units = "in"
  )
}

# ------------------------------------------------------------------
# Figure 2: AF PIP-linked motif enrichment
# ------------------------------------------------------------------

P_GEN <- GEN[
  model_converged ==
    TRUE &
    is.finite(
      OR_motif
    ) &
    is.finite(
      FDR_BH
    )
]

P_GEN <- head(
  P_GEN[
    order(
      FDR_BH,
      -OR_motif
    )
  ],
  TOP_N_GENETIC_PLOT
)

if (nrow(
  P_GEN
)) {

  P_GEN[
    ,
    label :=
      paste0(
        motif_name,
        " [",
        motif_id,
        "]"
      )
  ]

  P_GEN[
    ,
    label :=
      factor(
        label,
        levels =
          rev(
            label
          )
      )
  ]

  p2 <- ggplot(
    P_GEN,
    aes(
      x = OR_motif,
      y = label
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = 2,
      linewidth = 0.45
    ) +
    geom_point(
      size = 2.2
    ) +
    labs(
      x = "Adjusted OR for AF PIP-linked peak",
      y = NULL
    ) +
    theme_classic(
      base_size = 10.5
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(
        colour = "black"
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_AF_PIP_motif_adjusted_OR.tiff"
    ),
    p2,
    width = 6.2,
    height = 6.0,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_AF_PIP_motif_adjusted_OR.pdf"
    ),
    p2,
    width = 6.2,
    height = 6.0,
    units = "in"
  )
}

# ------------------------------------------------------------------
# Figure 3: integrated state NES vs genetic log2 OR
# ------------------------------------------------------------------

P_INT <- INTEGRATED[
  is.finite(
    NES
  ) &
    is.finite(
      OR_motif
    ) &
    OR_motif >
      0
]

if (nrow(
  P_INT
)) {

  P_INT[
    ,
    category :=
      fifelse(
        convergent_state_and_genetic %
          in% TRUE,
        "Convergent",
        fifelse(
          state_robust %
            in% TRUE,
          "State only",
          fifelse(
            genetic_enriched %
              in% TRUE,
            "Genetic only",
            "Other"
          )
        )
      )
  ]

  p3 <- ggplot(
    P_INT,
    aes(
      x = NES,
      y = log2(
        OR_motif
      )
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = 2
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      linetype = 2
    ) +
    geom_point(
      aes(
        shape = category
      ),
      size = 1.7,
      alpha = 0.70
    ) +
    labs(
      x = "CM High-vs-Low TRS motif NES",
      y = "log2 adjusted OR for AF PIP-linked peak"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(
        colour = "black"
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_state_vs_genetic_motif_convergence.tiff"
    ),
    p3,
    width = 5.5,
    height = 4.8,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4C_state_vs_genetic_motif_convergence.pdf"
    ),
    p3,
    width = 5.5,
    height = 4.8,
    units = "in"
  )
}

# ==============================================================================
# 14. FINAL SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G4C COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Tested peaks: ",
  format(
    nrow(DA),
    big.mark = ","
  ),
  "\n",
  "JASPAR2024 human CORE motifs: ",
  nrow(MOTIF_META),
  "\n",
  "State-robust motifs: ",
  N_STATE_ROBUST,
  "\n",
  "AF-genetic-enriched motifs: ",
  N_GENETIC_ENRICHED,
  "\n",
  "Convergent state+genetic motifs: ",
  N_CONVERGENT,
  "\n",
  "TF->robustDAR->gene edges: ",
  nrow(NETWORK),
  "\n\n",
  sep = ""
)

cat(
  "Readiness:\n"
)

print(
  READINESS
)

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  STATE_DIR,
  "\n",
  GENETIC_DIR,
  "\n",
  INTEGRATION_DIR,
  "\n",
  NETWORK_DIR,
  "\n",
  sep = ""
)
