# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G4B_CM_INTERNAL_TRS_CHROMATIN_PSEUDOBULK(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G4B — CM-INTERNAL AF-TRS CHROMATIN STATE DISSECTION
# Project: AF–HFpEF–BMI–OSA shared genetics
# Dataset: GSE238242 left atrial appendage snATAC
#
# PURPOSE
#   After freezing the cell-type result (CM is the dominant AF-relevant cell
#   population), characterize chromatin accessibility associated with high
#   vs low AF trait relevance WITHIN cardiomyocytes.
#
# IMPORTANT INTERPRETATION
#   TRS is derived from the same scATAC data. Therefore, this is a downstream
#   regulatory characterization / mechanism-generating analysis, NOT an
#   independent validation of SCAVENGE.
#
# PRIMARY DESIGN
#   1) Keep CM only.
#   2) Within EACH donor, select exactly:
#        top 25% CM by frozen V3 AF_TRS  = High
#        bottom 25% CM by frozen V3 AF_TRS = Low
#      Middle 50% are not used in the primary differential-accessibility test.
#   3) Aggregate raw ATAC counts to donor × state pseudobulk:
#        7 donors × 2 states = 14 pseudobulk samples.
#   4) edgeR paired quasi-likelihood model:
#        ~ donor + state
#      coefficient: stateHigh
#   5) Annotate peaks to nearest hg38 TSS.
#   6) Overlay frozen AF fine-mapping peak weights from STEP11G3A.
#   7) Quantify donor-direction consistency for every tested peak.
#
# PRIMARY DAR DEFINITION
#   FDR < 0.05 AND |log2FC| >= 0.5
#
# ROBUST DAR DEFINITION
#   primary DAR + same direction in >= 5/7 donors
#
# OUTPUT
#   07_STEP11G4B_CM_TRS_CHROMATIN/
#
# NEXT STEP AFTER REVIEW
#   STEP11G4C: motif / TF enrichment using robust High-TRS CM DARs,
#   plus regulatory interpretation and final publication panels.
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
EXPECTED_N_CM <- 3049L
EXPECTED_N_DONORS <- 7L
EXPECTED_N_PEAKS <- 212084L

CM_LABEL <- "CM"

EXTREME_FRACTION <- 0.25

DAR_FDR <- 0.05
DAR_ABS_LOGFC <- 0.5

ROBUST_MIN_DONORS <- 5L

TOP_PEAKS_FOR_TABLE <- 500L
TOP_GENES_FOR_FIGURE <- 20L

# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

# Frozen G3A peak-by-cell object.
G3A_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "04_GCHROMVAR"
)

RSE_FILE <- file.path(
  G3A_ROOT,
  "01_SCATAC_MATRIX",
  "GSE238242_SCAVENGE_peak_by_cell_GRCh38.rds"
)

TRAIT_WEIGHT_FILE <- file.path(
  G3A_ROOT,
  "03_GCHROMVAR",
  "AF_gchromVAR_trait_weights.rds"
)

# Frozen official-source-equivalent V3 primary SCAVENGE result.
G3B_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "05_SCAVENGE_NETWORK_V3_OFFICIAL_EQUIV"
)

G3B_READINESS <- file.path(
  G3B_ROOT,
  "00_QC",
  "STEP11G3B_readiness.csv"
)

PRIMARY_CELL_FILE <- file.path(
  G3B_ROOT,
  "05_FINAL",
  "STEP11G3B_AF_SCAVENGE_cell_results.csv.gz"
)

# G4A robustness must already be closed before G4B.
G4A_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "06_STEP11G4A_LOW_MATCH_SENSITIVITY"
)

G4A_READINESS <- file.path(
  G4A_ROOT,
  "00_QC",
  "STEP11G4A_readiness.csv"
)

# New output.
OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "07_STEP11G4B_CM_TRS_CHROMATIN"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

STATE_DIR <- file.path(
  OUT_ROOT,
  "01_CM_STATE"
)

PB_DIR <- file.path(
  OUT_ROOT,
  "02_PSEUDOBULK"
)

DA_DIR <- file.path(
  OUT_ROOT,
  "03_DIFFERENTIAL_ACCESSIBILITY"
)

GENE_DIR <- file.path(
  OUT_ROOT,
  "04_GENE_ANNOTATION"
)

FIG_DIR <- file.path(
  OUT_ROOT,
  "05_FIGURES"
)

for (d in c(
  OUT_ROOT,
  QC_DIR,
  STATE_DIR,
  PB_DIR,
  DA_DIR,
  GENE_DIR,
  FIG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

CM_STATE_FILE <- file.path(
  STATE_DIR,
  "STEP11G4B_CM_high_low_cell_manifest.csv.gz"
)

PB_COUNTS_FILE <- file.path(
  PB_DIR,
  "STEP11G4B_CM_high_low_pseudobulk_counts.rds"
)

PB_META_FILE <- file.path(
  PB_DIR,
  "STEP11G4B_CM_high_low_pseudobulk_sample_metadata.csv"
)

DA_ALL_FILE <- file.path(
  DA_DIR,
  "STEP11G4B_CM_TRS_DA_all_tested_peaks.csv.gz"
)

DA_SIG_FILE <- file.path(
  DA_DIR,
  "STEP11G4B_CM_TRS_DAR_FDR05_logFC05.csv.gz"
)

DA_ROBUST_FILE <- file.path(
  DA_DIR,
  "STEP11G4B_CM_TRS_robust_DAR_5of7donors.csv.gz"
)

DA_PIP_FILE <- file.path(
  DA_DIR,
  "STEP11G4B_CM_TRS_robust_PIP_linked_DAR.csv.gz"
)

GENE_ALL_FILE <- file.path(
  GENE_DIR,
  "STEP11G4B_CM_TRS_candidate_gene_summary_all.csv"
)

GENE_HIGH_FILE <- file.path(
  GENE_DIR,
  "STEP11G4B_CM_TRS_candidate_genes_HighTRS.csv"
)

GENE_LOW_FILE <- file.path(
  GENE_DIR,
  "STEP11G4B_CM_TRS_candidate_genes_LowTRS.csv"
)

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

cran_pkgs <- c(
  "data.table",
  "Matrix",
  "ggplot2"
)

for (p in cran_pkgs) {

  if (!requireNamespace(
    p,
    quietly = TRUE
  )) {

    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }

  if (!requireNamespace(
    p,
    quietly = TRUE
  )) {
    stop(
      "Failed to install/load CRAN package: ",
      p
    )
  }
}

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
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
  "edgeR",
  "GenomicFeatures",
  "AnnotationDbi",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

missing_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(
  missing_bioc
)) {

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

still_missing <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(
  still_missing
)) {
  stop(
    paste0(
      "Missing Bioconductor package(s): ",
      paste(
        still_missing,
        collapse = ", "
      )
    )
  )
}

library(data.table)
library(Matrix)
library(ggplot2)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

bool_pass <- function(x) {
  x %in% c(
    TRUE,
    "TRUE",
    1,
    "1"
  )
}

safe_spearman <- function(x, y) {
  suppressWarnings(
    cor(
      x,
      y,
      method = "spearman",
      use = "complete.obs"
    )
  )
}

# ==============================================================================
# 4. HARD INPUT CHECK / FROZEN ANALYSIS CHECK
# ==============================================================================

required_inputs <- c(
  RSE_FILE,
  TRAIT_WEIGHT_FILE,
  G3B_READINESS,
  PRIMARY_CELL_FILE,
  G4A_READINESS
)

if (!all(
  file.exists(
    required_inputs
  )
)) {
  stop(
    paste0(
      "Missing required input(s):\n",
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

READY_G3B <- fread(
  G3B_READINESS
)

READY_G4A <- fread(
  G4A_READINESS
)

if (!all(
  bool_pass(
    READY_G3B$pass
  )
)) {
  stop(
    "Frozen STEP11G3B V3 is not fully PASS."
  )
}

if (!all(
  bool_pass(
    READY_G4A$pass
  )
)) {
  stop(
    "STEP11G4A sensitivity is not fully PASS."
  )
}

# ==============================================================================
# 5. LOAD FROZEN PRIMARY CELL RESULTS + PEAK-BY-CELL MATRIX
# ==============================================================================

PRIMARY <- fread(
  PRIMARY_CELL_FILE
)

required_primary_cols <- c(
  "cell_id",
  "cell_type",
  "donor",
  "rhythm",
  "AF_gchromVAR_Z",
  "AF_TRS",
  "AF_SCAVENGE_significant_cell"
)

if (!all(
  required_primary_cols %in%
    names(
      PRIMARY
    )
)) {
  stop(
    paste0(
      "Primary V3 cell table lacks: ",
      paste(
        setdiff(
          required_primary_cols,
          names(
            PRIMARY
          )
        ),
        collapse = ", "
      )
    )
  )
}

if (
  nrow(
    PRIMARY
  ) !=
    EXPECTED_N_CELLS
) {
  stop(
    "Expected 11,986 primary cells; observed ",
    nrow(
      PRIMARY
    )
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
    "Frozen peak-by-cell counts are not sparse."
  )
}

if (
  nrow(
    COUNTS
  ) !=
    EXPECTED_N_PEAKS ||
  ncol(
    COUNTS
  ) !=
    EXPECTED_N_CELLS
) {
  stop(
    "Frozen peak-by-cell matrix dimensions differ from 212084 x 11986."
  )
}

idx <- match(
  colnames(
    COUNTS
  ),
  PRIMARY$cell_id
)

if (anyNA(
  idx
)) {
  stop(
    "Could not align V3 cell results to the frozen ATAC matrix."
  )
}

PRIMARY <- PRIMARY[
  idx
]

if (!identical(
  PRIMARY$cell_id,
  colnames(
    COUNTS
  )
)) {
  stop(
    "Cell order alignment failed."
  )
}

# ==============================================================================
# 6. DEFINE DONOR-BALANCED HIGH / LOW TRS CM STATES
# ==============================================================================

CM <- copy(
  PRIMARY[
    cell_type ==
      CM_LABEL
  ]
)

if (
  nrow(
    CM
  ) !=
    EXPECTED_N_CM
) {
  stop(
    "Expected 3,049 CM cells; observed ",
    nrow(
      CM
    )
  )
}

if (
  uniqueN(
    CM$donor
  ) !=
    EXPECTED_N_DONORS
) {
  stop(
    "Expected CM cells from 7 donors."
  )
}

# Number of cells available per donor.
CM[
  ,
  n_CM_donor :=
    .N,
  by = donor
]

CM[
  ,
  n_extreme :=
    pmax(
      1L,
      floor(
        n_CM_donor *
          EXTREME_FRACTION
      )
    )
]

# High rank: deterministic ordering by TRS, then gchromVAR Z, then cell id.
setorder(
  CM,
  donor,
  -AF_TRS,
  -AF_gchromVAR_Z,
  cell_id
)

CM[
  ,
  rank_high :=
    seq_len(
      .N
    ),
  by = donor
]

# Low rank.
setorder(
  CM,
  donor,
  AF_TRS,
  AF_gchromVAR_Z,
  cell_id
)

CM[
  ,
  rank_low :=
    seq_len(
      .N
    ),
  by = donor
]

CM[
  ,
  TRS_state :=
    fifelse(
      rank_high <=
        n_extreme,
      "High",
      fifelse(
        rank_low <=
          n_extreme,
        "Low",
        "Middle"
      )
    )
]

CM[
  ,
  TRS_state :=
    factor(
      TRS_state,
      levels = c(
        "Low",
        "Middle",
        "High"
      )
    )
]

STATE_QC <- CM[
  ,
  .(
    n_CM =
      .N,
    n_High =
      sum(
        TRS_state ==
          "High"
      ),
    n_Low =
      sum(
        TRS_state ==
          "Low"
      ),
    n_Middle =
      sum(
        TRS_state ==
          "Middle"
      ),
    High_median_TRS =
      median(
        AF_TRS[
          TRS_state ==
            "High"
        ]
      ),
    Low_median_TRS =
      median(
        AF_TRS[
          TRS_state ==
            "Low"
        ]
      ),
    High_median_Z =
      median(
        AF_gchromVAR_Z[
          TRS_state ==
            "High"
        ]
      ),
    Low_median_Z =
      median(
        AF_gchromVAR_Z[
          TRS_state ==
            "Low"
        ]
      )
  ),
  by = .(
    donor,
    rhythm
  )
]

if (any(
  STATE_QC$n_High !=
    STATE_QC$n_Low
)) {
  stop(
    "High and Low CM counts are not balanced within donor."
  )
}

if (any(
  STATE_QC$High_median_TRS <=
    STATE_QC$Low_median_TRS
)) {
  stop(
    "At least one donor does not have High TRS > Low TRS."
  )
}

fwrite(
  STATE_QC,
  file.path(
    QC_DIR,
    "STEP11G4B_CM_state_definition_QC.csv"
  )
)

fwrite(
  CM,
  CM_STATE_FILE,
  compress = "gzip"
)

# ==============================================================================
# 7. BUILD 14 DONOR x STATE PSEUDOBULK SAMPLES
# ==============================================================================

CM_EXTREME <- CM[
  TRS_state %in%
    c(
      "Low",
      "High"
    )
]

CM_EXTREME[
  ,
  sample_id :=
    paste(
      donor,
      as.character(
        TRS_state
      ),
      sep = "__"
    )
]

cm_col_idx <- match(
  CM_EXTREME$cell_id,
  colnames(
    COUNTS
  )
)

if (anyNA(
  cm_col_idx
)) {
  stop(
    "Could not align High/Low CM cells to ATAC matrix."
  )
}

COUNTS_CM <- COUNTS[
  ,
  cm_col_idx,
  drop = FALSE
]

colnames(
  COUNTS_CM
) <- CM_EXTREME$cell_id

SAMPLES <- unique(
  CM_EXTREME[
    ,
    .(
      sample_id,
      donor,
      rhythm,
      state =
        as.character(
          TRS_state
        )
    )
  ]
)

SAMPLES[
  ,
  state :=
    factor(
      state,
      levels = c(
        "Low",
        "High"
      )
    )
]

SAMPLES[
  ,
  donor :=
    factor(
      donor
    )
]

setorder(
  SAMPLES,
  donor,
  state
)

if (
  nrow(
    SAMPLES
  ) !=
    2L *
      EXPECTED_N_DONORS
) {
  stop(
    "Expected 14 donor x state pseudobulk samples."
  )
}

PB_LIST <- vector(
  "list",
  nrow(
    SAMPLES
  )
)

for (
  j in seq_len(
    nrow(
      SAMPLES
    )
  )
) {

  sid <- SAMPLES$sample_id[
    j
  ]

  cell_ids <- CM_EXTREME[
    sample_id ==
      sid,
    cell_id
  ]

  jj <- match(
    cell_ids,
    colnames(
      COUNTS_CM
    )
  )

  PB_LIST[[
    j
  ]] <- Matrix::rowSums(
    COUNTS_CM[
      ,
      jj,
      drop = FALSE
    ]
  )
}

PB <- do.call(
  cbind,
  PB_LIST
)

storage.mode(
  PB
) <- "integer"

colnames(
  PB
) <- SAMPLES$sample_id

rownames(
  PB
) <- rownames(
  COUNTS
)

SAMPLES[
  ,
  n_cells :=
    vapply(
      sample_id,
      function(x) {
        sum(
          CM_EXTREME$sample_id ==
            x
        )
      },
      integer(1)
    )
]

SAMPLES[
  ,
  library_size_raw :=
    colSums(
      PB
    )
]

if (any(
  SAMPLES$library_size_raw <=
    0
)) {
  stop(
    "At least one pseudobulk sample has zero library size."
  )
}

saveRDS(
  PB,
  PB_COUNTS_FILE,
  compress = TRUE
)

fwrite(
  SAMPLES,
  PB_META_FILE
)

# ==============================================================================
# 8. PAIRED PSEUDOBULK DIFFERENTIAL ACCESSIBILITY — edgeR QL
# ==============================================================================

Y <- edgeR::DGEList(
  counts = PB,
  samples = as.data.frame(
    SAMPLES
  )
)

DESIGN <- model.matrix(
  ~ donor + state,
  data = SAMPLES
)

if (
  qr(
    DESIGN
  )$rank !=
    ncol(
      DESIGN
    )
) {
  stop(
    "Paired donor + state design matrix is not full rank."
  )
}

if (!"stateHigh" %in%
    colnames(
      DESIGN
    )) {
  stop(
    "Expected stateHigh coefficient was not created."
  )
}

KEEP_PEAK <- edgeR::filterByExpr(
  Y,
  design = DESIGN
)

if (sum(
  KEEP_PEAK
) <
    1000L) {
  stop(
    "Too few peaks passed edgeR expression/accessibility filtering."
  )
}

YF <- Y[
  KEEP_PEAK,
  ,
  keep.lib.sizes = FALSE
]

YF <- edgeR::calcNormFactors(
  YF,
  method = "TMM"
)

YF <- edgeR::estimateDisp(
  YF,
  DESIGN,
  robust = TRUE
)

FIT <- edgeR::glmQLFit(
  YF,
  DESIGN,
  robust = TRUE
)

QLF <- edgeR::glmQLFTest(
  FIT,
  coef = "stateHigh"
)

TT <- edgeR::topTags(
  QLF,
  n = Inf,
  sort.by = "PValue"
)$table

DA <- as.data.table(
  TT,
  keep.rownames = "peak_id"
)

DA[
  ,
  direction :=
    fifelse(
      logFC >
        0,
      "HighTRS_open",
      "LowTRS_open"
    )
]

DA[
  ,
  primary_DAR :=
    FDR <
      DAR_FDR &
    abs(
      logFC
    ) >=
      DAR_ABS_LOGFC
]

# ==============================================================================
# 9. DONOR-DIRECTION CONSISTENCY
# ==============================================================================

LOGCPM <- edgeR::cpm(
  YF,
  log = TRUE,
  prior.count = 2
)

DONORS <- levels(
  SAMPLES$donor
)

DIFF_LIST <- lapply(
  DONORS,
  function(d) {

    hi <- which(
      SAMPLES$donor ==
        d &
      SAMPLES$state ==
        "High"
    )

    lo <- which(
      SAMPLES$donor ==
        d &
      SAMPLES$state ==
        "Low"
    )

    if (
      length(
        hi
      ) != 1L ||
      length(
        lo
      ) != 1L
    ) {
      stop(
        "Expected exactly one High and one Low pseudobulk for donor ",
        d
      )
    }

    LOGCPM[
      ,
      hi
    ] -
      LOGCPM[
        ,
        lo
      ]
  }
)

DIFF_MAT <- do.call(
  cbind,
  DIFF_LIST
)

colnames(
  DIFF_MAT
) <- paste0(
  "delta_logCPM_",
  DONORS
)

DA_ROW <- match(
  DA$peak_id,
  rownames(
    DIFF_MAT
  )
)

if (anyNA(
  DA_ROW
)) {
  stop(
    "Could not align differential-accessibility table to donor delta matrix."
  )
}

DIFF_MAT <- DIFF_MAT[
  DA_ROW,
  ,
  drop = FALSE
]

DA[
  ,
  n_donors_positive :=
    rowSums(
      DIFF_MAT >
        0
    )
]

DA[
  ,
  n_donors_negative :=
    rowSums(
      DIFF_MAT <
        0
    )
]

DA[
  ,
  n_donors_concordant :=
    fifelse(
      logFC >
        0,
      n_donors_positive,
      n_donors_negative
    )
]

DA[
  ,
  donor_consistency_fraction :=
    n_donors_concordant /
      EXPECTED_N_DONORS
]

DA[
  ,
  robust_DAR :=
    primary_DAR &
      n_donors_concordant >=
        ROBUST_MIN_DONORS
]

# Add individual donor deltas for source-data auditing.
for (
  j in seq_along(
    DONORS
  )
) {

  new_col <- paste0(
    "delta_logCPM_",
    DONORS[
      j
    ]
  )

  DA[
    ,
    (new_col) :=
      DIFF_MAT[
        ,
        j
      ]
  ]
}

# ==============================================================================
# 10. ADD GRCh38 PEAK COORDINATES
# ==============================================================================

PEAK_GR <- SummarizedExperiment::rowRanges(
  SE
)

PEAK_IDS <- rownames(
  SE
)

if (is.null(
  PEAK_IDS
)) {
  PEAK_IDS <- rownames(
    COUNTS
  )
}

coord_idx <- match(
  DA$peak_id,
  PEAK_IDS
)

if (anyNA(
  coord_idx
)) {
  stop(
    "Could not align edgeR peaks to rowRanges."
  )
}

DA[
  ,
  `:=`(
    chr =
      as.character(
        GenomicRanges::seqnames(
          PEAK_GR[
            coord_idx
          ]
        )
      ),
    start =
      GenomicRanges::start(
        PEAK_GR[
          coord_idx
        ]
      ),
    end =
      GenomicRanges::end(
        PEAK_GR[
          coord_idx
        ]
      )
  )
]

# ==============================================================================
# 11. ADD FROZEN AF PIP PEAK WEIGHT
# ==============================================================================

TRAIT_WEIGHT <- readRDS(
  TRAIT_WEIGHT_FILE
)

if (inherits(
  TRAIT_WEIGHT,
  "SummarizedExperiment"
)) {

  W <- SummarizedExperiment::assay(
    TRAIT_WEIGHT,
    1
  )

} else if (
  inherits(
    TRAIT_WEIGHT,
    "Matrix"
  ) ||
  is.matrix(
    TRAIT_WEIGHT
  )
) {

  W <- TRAIT_WEIGHT

} else {

  stop(
    "Unsupported frozen AF trait-weight object."
  )
}

W <- as.numeric(
  W[
    ,
    1
  ]
)

if (
  length(
    W
  ) !=
    EXPECTED_N_PEAKS
) {
  stop(
    "Frozen AF PIP peak-weight vector length mismatch."
  )
}

DA[
  ,
  AF_PIP_peak_weight :=
    W[
      coord_idx
    ]
]

DA[
  ,
  AF_PIP_linked :=
    is.finite(
      AF_PIP_peak_weight
    ) &
      AF_PIP_peak_weight >
        0
]

DA[
  ,
  AF_PIP_x_abs_logFC :=
    AF_PIP_peak_weight *
      abs(
        logFC
      )
]

# ==============================================================================
# 12. NEAREST hg38 TSS ANNOTATION
# ==============================================================================

TXDB <-
  TxDb.Hsapiens.UCSC.hg38.knownGene::
  TxDb.Hsapiens.UCSC.hg38.knownGene

GENES <- GenomicFeatures::genes(
  TXDB,
  single.strand.genes.only = TRUE
)

gene_entrez <- names(
  GENES
)

gene_strand <- as.character(
  GenomicRanges::strand(
    GENES
  )
)

gene_tss <- ifelse(
  gene_strand ==
    "-",
  GenomicRanges::end(
    GENES
  ),
  GenomicRanges::start(
    GENES
  )
)

TSS_GR <- GenomicRanges::GRanges(
  seqnames =
    GenomicRanges::seqnames(
      GENES
    ),
  ranges = IRanges::IRanges(
    start = gene_tss,
    width = 1L
  ),
  strand =
    GenomicRanges::strand(
      GENES
    )
)

S4Vectors::mcols(
  TSS_GR
)$ENTREZID <- gene_entrez

PEAK_TEST_GR <- PEAK_GR[
  coord_idx
]

PEAK_MID <- GenomicRanges::resize(
  PEAK_TEST_GR,
  width = 1L,
  fix = "center"
)

NEAR <- GenomicRanges::distanceToNearest(
  PEAK_MID,
  TSS_GR,
  ignore.strand = TRUE
)

qhit <- S4Vectors::queryHits(
  NEAR
)

shit <- S4Vectors::subjectHits(
  NEAR
)

NEAREST_ENTREZ <- rep(
  NA_character_,
  nrow(
    DA
  )
)

DIST_TSS <- rep(
  NA_integer_,
  nrow(
    DA
  )
)

NEAREST_ENTREZ[
  qhit
] <- S4Vectors::mcols(
  TSS_GR
)$ENTREZID[
  shit
]

DIST_TSS[
  qhit
] <- S4Vectors::mcols(
  NEAR
)$distance

unique_entrez <- unique(
  NEAREST_ENTREZ[
    !is.na(
      NEAREST_ENTREZ
    )
  ]
)

SYMBOL_MAP <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

GENENAME_MAP <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique_entrez,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)

DA[
  ,
  nearest_ENTREZID :=
    NEAREST_ENTREZ
]

DA[
  ,
  nearest_gene_symbol :=
    unname(
      SYMBOL_MAP[
        nearest_ENTREZID
      ]
    )
]

DA[
  ,
  nearest_gene_name :=
    unname(
      GENENAME_MAP[
        nearest_ENTREZID
      ]
    )
]

DA[
  ,
  distance_to_nearest_TSS :=
    DIST_TSS
]

# ==============================================================================
# 13. FINAL DAR TABLES
# ==============================================================================

DA[
  ,
  abs_logFC :=
    abs(
      logFC
    )
]

setorder(
  DA,
  FDR,
  -abs_logFC
)

fwrite(
  DA,
  DA_ALL_FILE,
  compress = "gzip"
)

SIG <- DA[
  primary_DAR ==
    TRUE
]

ROBUST <- DA[
  robust_DAR ==
    TRUE
]

ROBUST_PIP <- ROBUST[
  AF_PIP_linked ==
    TRUE
]

fwrite(
  SIG,
  DA_SIG_FILE,
  compress = "gzip"
)

fwrite(
  ROBUST,
  DA_ROBUST_FILE,
  compress = "gzip"
)

fwrite(
  ROBUST_PIP,
  DA_PIP_FILE,
  compress = "gzip"
)

# ==============================================================================
# 14. CANDIDATE GENE SUMMARIES
# ==============================================================================

GENE_BASE <- ROBUST[
  !is.na(
    nearest_gene_symbol
  ) &
    nzchar(
      nearest_gene_symbol
    )
]

if (nrow(
  GENE_BASE
)) {

  GENE_SUM <- GENE_BASE[
    ,
    .(
      n_robust_DARs =
        .N,
      n_HighTRS_open_DARs =
        sum(
          logFC >
            0
        ),
      n_LowTRS_open_DARs =
        sum(
          logFC <
            0
        ),
      n_AF_PIP_linked_DARs =
        sum(
          AF_PIP_linked
        ),
      min_FDR =
        min(
          FDR,
          na.rm = TRUE
        ),
      max_abs_logFC =
        max(
          abs(
            logFC
          ),
          na.rm = TRUE
        ),
      max_HighTRS_logFC =
        if (
          any(
            logFC >
              0
          )
        ) {
          max(
            logFC[
              logFC >
                0
            ]
          )
        } else {
          NA_real_
        },
      min_LowTRS_logFC =
        if (
          any(
            logFC <
              0
          )
        ) {
          min(
            logFC[
              logFC <
                0
            ]
          )
        } else {
          NA_real_
        },
      closest_TSS_distance =
        min(
          distance_to_nearest_TSS,
          na.rm = TRUE
        ),
      max_AF_PIP_peak_weight =
        max(
          AF_PIP_peak_weight,
          na.rm = TRUE
        ),
      max_AF_PIP_x_abs_logFC =
        max(
          AF_PIP_x_abs_logFC,
          na.rm = TRUE
        ),
      max_donor_consistency_fraction =
        max(
          donor_consistency_fraction,
          na.rm = TRUE
        )
    ),
    by = .(
      nearest_ENTREZID,
      nearest_gene_symbol,
      nearest_gene_name
    )
  ]

  setorder(
    GENE_SUM,
    -n_AF_PIP_linked_DARs,
    -max_AF_PIP_x_abs_logFC,
    min_FDR
  )

} else {

  GENE_SUM <- data.table()
}

fwrite(
  GENE_SUM,
  GENE_ALL_FILE
)

GENE_HIGH <- if (nrow(
  GENE_SUM
)) {

  GENE_SUM[
    n_HighTRS_open_DARs >
      0
  ][
    order(
      -n_AF_PIP_linked_DARs,
      -max_AF_PIP_x_abs_logFC,
      -max_HighTRS_logFC,
      min_FDR
    )
  ]

} else {
  data.table()
}

GENE_LOW <- if (nrow(
  GENE_SUM
)) {

  GENE_SUM[
    n_LowTRS_open_DARs >
      0
  ][
    order(
      -n_AF_PIP_linked_DARs,
      -max_AF_PIP_x_abs_logFC,
      min_LowTRS_logFC,
      min_FDR
    )
  ]

} else {
  data.table()
}

fwrite(
  GENE_HIGH,
  GENE_HIGH_FILE
)

fwrite(
  GENE_LOW,
  GENE_LOW_FILE
)

# ==============================================================================
# 15. QC SUMMARIES
# ==============================================================================

DA_QC <- data.table(
  metric = c(
    "CM_cells_total",
    "donors",
    "High_CM_cells",
    "Low_CM_cells",
    "Middle_CM_cells",
    "pseudobulk_samples",
    "peaks_total",
    "peaks_tested_edgeR",
    "primary_DARs",
    "primary_HighTRS_open_DARs",
    "primary_LowTRS_open_DARs",
    "robust_DARs_5of7",
    "robust_HighTRS_open_DARs",
    "robust_LowTRS_open_DARs",
    "robust_PIP_linked_DARs",
    "candidate_genes_robust",
    "candidate_genes_with_PIP_linked_DAR",
    "edgeR_common_dispersion"
  ),
  value = c(
    nrow(
      CM
    ),
    uniqueN(
      CM$donor
    ),
    sum(
      CM$TRS_state ==
        "High"
    ),
    sum(
      CM$TRS_state ==
        "Low"
    ),
    sum(
      CM$TRS_state ==
        "Middle"
    ),
    nrow(
      SAMPLES
    ),
    nrow(
      COUNTS
    ),
    nrow(
      DA
    ),
    nrow(
      SIG
    ),
    sum(
      SIG$logFC >
        0
    ),
    sum(
      SIG$logFC <
        0
    ),
    nrow(
      ROBUST
    ),
    sum(
      ROBUST$logFC >
        0
    ),
    sum(
      ROBUST$logFC <
        0
    ),
    nrow(
      ROBUST_PIP
    ),
    nrow(
      GENE_SUM
    ),
    if (
      nrow(
        GENE_SUM
      )
    ) {
      sum(
        GENE_SUM$n_AF_PIP_linked_DARs >
          0
      )
    } else {
      0
    },
    YF$common.dispersion
  )
)

fwrite(
  DA_QC,
  file.path(
    QC_DIR,
    "STEP11G4B_differential_accessibility_QC.csv"
  )
)

PB_QC <- copy(
  SAMPLES
)

PB_QC[
  ,
  effective_library_size :=
    YF$samples$lib.size *
      YF$samples$norm.factors
]

fwrite(
  PB_QC,
  file.path(
    QC_DIR,
    "STEP11G4B_pseudobulk_sample_QC.csv"
  )
)

# ==============================================================================
# 16. FIGURES — SINGLE FILE PER FIGURE, NO LEGEND INSIDE
# ==============================================================================

# 16A. TRS state distribution by donor.
PLOT_CM <- CM[
  TRS_state %in%
    c(
      "Low",
      "High"
    )
]

p_state <- ggplot(
  PLOT_CM,
  aes(
    x = donor,
    y = AF_TRS,
    fill = TRS_state
  )
) +
  geom_boxplot(
    width = 0.72,
    outlier.shape = NA,
    linewidth = 0.45,
    position = position_dodge(
      width = 0.75
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low" = "#7E9CB7",
      "High" = "#C8574C"
    )
  ) +
  labs(
    x = NULL,
    y = "AF trait relevance score (TRS)"
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
    "Figure_STEP11G4B_CM_HighLow_TRS_by_donor.tiff"
  ),
  p_state,
  width = 6.2,
  height = 4.4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G4B_CM_HighLow_TRS_by_donor.pdf"
  ),
  p_state,
  width = 6.2,
  height = 4.4,
  units = "in"
)

# 16B. Volcano.
DA[
  ,
  volcano_group :=
    fifelse(
      robust_DAR &
        logFC >
          0,
      "HighTRS-open robust",
      fifelse(
        robust_DAR &
          logFC <
            0,
        "LowTRS-open robust",
        "Other"
      )
    )
]

VOLC_COL <- c(
  "Other" = "grey78",
  "HighTRS-open robust" = "#B34D4D",
  "LowTRS-open robust" = "#557DA1"
)

p_volcano <- ggplot(
  DA,
  aes(
    x = logFC,
    y = -log10(
      pmax(
        FDR,
        1e-300
      )
    ),
    colour = volcano_group
  )
) +
  geom_point(
    size = 0.55,
    alpha = 0.58
  ) +
  geom_vline(
    xintercept = c(
      -DAR_ABS_LOGFC,
      DAR_ABS_LOGFC
    ),
    linetype = 2,
    linewidth = 0.45,
    colour = "grey45"
  ) +
  geom_hline(
    yintercept =
      -log10(
        DAR_FDR
      ),
    linetype = 2,
    linewidth = 0.45,
    colour = "grey45"
  ) +
  scale_colour_manual(
    values = VOLC_COL
  ) +
  labs(
    x = "log2 fold change (High TRS vs Low TRS)",
    y = expression(
      -log[10](
        FDR
      )
    )
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
    "Figure_STEP11G4B_CM_TRS_DAR_volcano.tiff"
  ),
  p_volcano,
  width = 5.5,
  height = 4.7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP11G4B_CM_TRS_DAR_volcano.pdf"
  ),
  p_volcano,
  width = 5.5,
  height = 4.7,
  units = "in"
)

# 16C. Top High-TRS candidate genes, if available.
if (
  nrow(
    GENE_HIGH
  ) >
    0
) {

  TOPG <- head(
    GENE_HIGH,
    TOP_GENES_FOR_FIGURE
  )

  TOPG[
    ,
    gene_plot :=
      factor(
        nearest_gene_symbol,
        levels = rev(
          nearest_gene_symbol
        )
      )
  ]

  p_gene <- ggplot(
    TOPG,
    aes(
      x = max_HighTRS_logFC,
      y = gene_plot
    )
  ) +
    geom_col(
      width = 0.70,
      linewidth = 0.35,
      fill = "#B34D4D",
      colour = "black"
    ) +
    labs(
      x = "Maximum robust peak log2FC",
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
      "Figure_STEP11G4B_top_HighTRS_candidate_genes.tiff"
    ),
    p_gene,
    width = 5.6,
    height = 5.2,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    file.path(
      FIG_DIR,
      "Figure_STEP11G4B_top_HighTRS_candidate_genes.pdf"
    ),
    p_gene,
    width = 5.6,
    height = 5.2,
    units = "in"
  )
}

# ==============================================================================
# 17. READINESS — TECHNICAL ONLY
# ==============================================================================

READINESS <- data.table(
  check = c(
    "STEP11G3B_V3_primary_fully_passed",
    "STEP11G4A_sensitivity_fully_passed",
    "3049_CM_cells_present",
    "7_CM_donors_present",
    "High_Low_balanced_within_each_donor",
    "High_median_TRS_gt_Low_in_each_donor",
    "14_pseudobulk_samples_created",
    "all_pseudobulk_library_sizes_positive",
    "paired_design_full_rank",
    "edgeR_tested_gt_1000_peaks",
    "all_tested_peaks_have_hg38_coordinates",
    "AF_PIP_peak_weights_aligned",
    "nearest_TSS_annotation_completed",
    "all_primary_output_tables_written",
    "STEP11G4B_complete"
  ),
  pass = c(
    all(
      bool_pass(
        READY_G3B$pass
      )
    ),
    all(
      bool_pass(
        READY_G4A$pass
      )
    ),
    nrow(
      CM
    ) ==
      EXPECTED_N_CM,
    uniqueN(
      CM$donor
    ) ==
      EXPECTED_N_DONORS,
    all(
      STATE_QC$n_High ==
        STATE_QC$n_Low
    ),
    all(
      STATE_QC$High_median_TRS >
        STATE_QC$Low_median_TRS
    ),
    nrow(
      SAMPLES
    ) ==
      14L,
    all(
      SAMPLES$library_size_raw >
        0
    ),
    qr(
      DESIGN
    )$rank ==
      ncol(
        DESIGN
      ),
    nrow(
      DA
    ) >
      1000L,
    all(
      !is.na(
        DA$chr
      )
    ) &&
      all(
        !is.na(
          DA$start
        )
      ) &&
      all(
        !is.na(
          DA$end
        )
      ),
    length(
      W
    ) ==
      EXPECTED_N_PEAKS,
    mean(
      !is.na(
        DA$nearest_gene_symbol
      )
    ) >
      0.80,
    all(
      file.exists(
        c(
          DA_ALL_FILE,
          DA_SIG_FILE,
          DA_ROBUST_FILE,
          DA_PIP_FILE,
          GENE_ALL_FILE,
          GENE_HIGH_FILE,
          GENE_LOW_FILE
        )
      )
    ),
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G4B_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "biological_question",
    "interpretation_scope",
    "cell_population",
    "state_definition",
    "extreme_fraction_per_donor",
    "middle_cells_in_DA_test",
    "statistical_unit",
    "pseudobulk_samples",
    "DA_method",
    "design",
    "contrast",
    "primary_DAR",
    "robust_DAR",
    "gene_annotation",
    "genome_build",
    "AF_genetic_overlay"
  ),
  value = c(
    "Chromatin accessibility associated with AF-TRS within cardiomyocytes",
    "Downstream characterization; not independent validation because TRS derives from the same scATAC data",
    "Author-defined CM",
    "Within-donor top versus bottom TRS quartile",
    as.character(
      EXTREME_FRACTION
    ),
    "Excluded",
    "Donor-level pseudobulk",
    "14 = 7 donors x High/Low",
    "edgeR quasi-likelihood negative-binomial GLM",
    "~ donor + state",
    "stateHigh",
    paste0(
      "FDR<",
      DAR_FDR,
      " and |log2FC|>=",
      DAR_ABS_LOGFC
    ),
    paste0(
      "Primary DAR + concordant direction in >=",
      ROBUST_MIN_DONORS,
      "/7 donors"
    ),
    "Nearest hg38 knownGene TSS + Entrez/SYMBOL mapping",
    "GRCh38/hg38",
    "Frozen STEP11G3A AF fine-mapping peak PIP weight"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G4B_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G4B_sessionInfo.txt"
  )
)

# ==============================================================================
# 18. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G4B COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "CM cells: ",
  nrow(
    CM
  ),
  "\n",
  "High TRS CM: ",
  sum(
    CM$TRS_state ==
      "High"
  ),
  "\n",
  "Low TRS CM: ",
  sum(
    CM$TRS_state ==
      "Low"
  ),
  "\n",
  "Pseudobulk samples: ",
  nrow(
    SAMPLES
  ),
  "\n",
  "Peaks tested: ",
  format(
    nrow(
      DA
    ),
    big.mark = ","
  ),
  "\n",
  "Primary DARs: ",
  format(
    nrow(
      SIG
    ),
    big.mark = ","
  ),
  "\n",
  "Robust DARs (>=5/7 donors): ",
  format(
    nrow(
      ROBUST
    ),
    big.mark = ","
  ),
  "\n",
  "Robust PIP-linked DARs: ",
  format(
    nrow(
      ROBUST_PIP
    ),
    big.mark = ","
  ),
  "\n",
  "Candidate genes: ",
  format(
    nrow(
      GENE_SUM
    ),
    big.mark = ","
  ),
  "\n\n",
  sep = ""
)

cat(
  "CM state QC:\n"
)
print(
  STATE_QC
)

cat(
  "\nDifferential-accessibility QC:\n"
)
print(
  DA_QC
)

cat(
  "\nReadiness:\n"
)
print(
  READINESS
)

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  DA_DIR,
  "\n",
  GENE_DIR,
  "\n",
  sep = ""
)
