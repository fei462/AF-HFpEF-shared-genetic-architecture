# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP13B_V2_FIX_CONDITION_SCP3342_HFPEF_CM_PSEUDOBULK.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP13B V2 — SCP3342 HFpEF vs CTRL CARDIOMYOCYTE PATIENT-LEVEL PSEUDOBULK DE
#
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# ROLE IN THE STUDY
#   Independent HFpEF disease-tissue validation after:
#     shared genetics -> regulatory loci/genes -> AF cellular localization (CM)
#     -> AF disease-tissue analysis (GSE255612)
#
# PRIMARY QUESTION
#   In human HFpEF ventricular-septal myocardium, what transcriptional changes
#   occur in author-defined cardiomyocytes at the PATIENT level?
#
# DATASET
#   SCP3342
#   48,866 nuclei
#   19 HFpEF + 24 non-failing CTRL donors
#   14 author-defined cell types
#
# FROZEN INPUTS FROM STEP13A
#   metadata cell type = original_cell_type
#   primary cell type = Cardiomyocyte
#   raw count source  = layers/raw_count_cellranger
#
# ANALYSIS DESIGN
#   - Raw CellRanger counts are streamed directly from H5AD CSR storage.
#   - Cardiomyocyte counts are summed within each donor.
#   - limma-voom pseudobulk DE.
#   - Primary design: ~ sex + pool + condition
#       (matches the published SCP3342 HFpEF-vs-control DE model).
#   - Sensitivity design: ~ sex + condition
#   - FDR is controlled within the cardiomyocyte transcriptome.
#
# BACKGROUND-CONTAMINATION HEURISTIC
#   A per-gene cardiomyocyte background heuristic is calculated from the
#   CellRanger raw-count matrix:
#       heuristic = bkg_prob * nontarget_prob
#   with:
#       bkg_prob = empirical cumulative rank of global gene UMI fraction
#       nontarget_prob = 1 - mean(PPV(expr>0), PPV(expr>1))
#   PPVs are standardized to equal target/non-target prevalence.
#   Genes with heuristic > 0.40 are flagged as high-background.
#
# IMPORTANT LIMITATION
#   The released H5AD layer named raw_count_cellbender contains continuous
#   non-integer values with the same scale as X and cannot safely be treated as
#   recoverable raw CellBender counts for patient-level summation.
#   Therefore STEP13B is a reproducible CellRanger-raw-count reanalysis.
#   It does NOT claim exact reproduction of the publication's dual-matrix
#   CellBender+CellRanger significance definition.
#
# NO CELL-LEVEL DE.
# NO RECLUSTERING.
# NO REANNOTATION.
# NO UPSTREAM CANDIDATE-GENE FILTERING.
#
# NEXT
#   STEP13C will integrate the FROZEN HFpEF-CM DE result with:
#     - frozen upstream genetics/regulatory candidate genes
#     - frozen AF-CM disease-tissue result from GSE255612
#
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))
set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

H5AD_FILE <- file.path(
  DATA_ROOT,
  "HFpEF_snRNAseq_single_cell_portal_10.14.2025.h5ad"
)

STEP13_ROOT <- file.path(
  DATA_ROOT,
  "STEP13_HFpEF_SCP3342_snRNA"
)

STEP13A_READINESS <- file.path(
  STEP13_ROOT,
  "00_QC",
  "STEP13A_readiness.csv"
)

COUNT_SOURCE_LOCK <- file.path(
  STEP13_ROOT,
  "00_QC",
  "STEP13A_count_source_lock.csv"
)

META_FILE <- file.path(
  STEP13_ROOT,
  "01_METADATA",
  "SCP3342_author_metadata_LOCKED.csv.gz"
)

VAR_FILE <- file.path(
  STEP13_ROOT,
  "01_METADATA",
  "SCP3342_var_gene_metadata_RAW.csv.gz"
)

PRIMARY_CELLTYPE <- "Cardiomyocyte"

EXPECTED_CELLS <- 48866L
EXPECTED_GENES <- 36601L
EXPECTED_DONORS <- 43L
EXPECTED_HFPEF <- 19L
EXPECTED_CTRL <- 24L

BLOCK_SIZE <- 500L

DE_FDR <- 0.05
BACKGROUND_HEURISTIC_MAX <- 0.40

# ==============================================================================
# 1. OUTPUT
# ==============================================================================

OUT_ROOT <- file.path(
  STEP13_ROOT,
  "STEP13B_HFpEF_vs_CTRL_CM_DE"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
PB_DIR <- file.path(OUT_ROOT, "01_PSEUDOBULK")
DE_DIR <- file.path(OUT_ROOT, "02_PRIMARY_DE")
SENS_DIR <- file.path(OUT_ROOT, "03_SENSITIVITY")
FIG_DIR <- file.path(OUT_ROOT, "04_FIGURES")
SOURCE_DIR <- file.path(OUT_ROOT, "05_SOURCE_DATA")

for (d in c(
  OUT_ROOT,
  QC_DIR,
  PB_DIR,
  DE_DIR,
  SENS_DIR,
  FIG_DIR,
  SOURCE_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PB_FILE <- file.path(
  PB_DIR,
  "STEP13B_Cardiomyocyte_CellRanger_patient_pseudobulk_counts.rds"
)

PB_META_FILE <- file.path(
  PB_DIR,
  "STEP13B_Cardiomyocyte_patient_metadata.csv"
)

BKG_FILE <- file.path(
  PB_DIR,
  "STEP13B_Cardiomyocyte_background_contamination_heuristic.csv.gz"
)

STREAM_QC_FILE <- file.path(
  QC_DIR,
  "STEP13B_H5AD_streaming_QC.csv"
)

ALL_DE_FILE <- file.path(
  DE_DIR,
  "STEP13B_Cardiomyocyte_HFpEF_vs_CTRL_all_genes.csv.gz"
)

FDR_DE_FILE <- file.path(
  DE_DIR,
  "STEP13B_Cardiomyocyte_HFpEF_vs_CTRL_FDR05.csv"
)

FDR_BKG_DE_FILE <- file.path(
  DE_DIR,
  "STEP13B_Cardiomyocyte_HFpEF_vs_CTRL_FDR05_background_le040.csv"
)

SENS_DE_FILE <- file.path(
  SENS_DIR,
  "STEP13B_Cardiomyocyte_HFpEF_vs_CTRL_sensitivity_sex_only.csv.gz"
)

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

cran_pkgs <- c(
  "data.table",
  "Matrix",
  "ggplot2",
  "ggrepel"
)

for (p in cran_pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }

  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Could not install/load CRAN package: ", p)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_pkgs <- c(
  "rhdf5",
  "edgeR",
  "limma",
  "AnnotationDbi",
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

if (length(missing_bioc)) {
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

if (length(still_missing)) {
  stop(
    "Missing Bioconductor package(s): ",
    paste(still_missing, collapse = ", ")
  )
}

library(data.table)
library(Matrix)
library(ggplot2)
library(ggrepel)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

bool_pass <- function(x) {
  x %in% c(TRUE, "TRUE", 1, "1")
}

norm_text <- function(x) {
  trimws(as.character(x))
}

norm_sex <- function(x) {

  y <- tolower(norm_text(x))

  out <- rep(NA_character_, length(y))

  out[
    y %in% c(
      "female",
      "f",
      "woman"
    )
  ] <- "F"

  out[
    y %in% c(
      "male",
      "m",
      "man"
    )
  ] <- "M"

  out
}

norm_condition <- function(x) {

  # Idempotent normalization:
  # STEP13A already freezes this field as "HFpEF" / "CTRL".
  # STEP13B may therefore receive either the original labels
  # ("hfpef", "control") or the already-standardized labels
  # ("HFpEF", "CTRL").  Both must map correctly.
  y <- tolower(
    norm_text(
      x
    )
  )

  out <- rep(
    NA_character_,
    length(y)
  )

  out[
    !is.na(y) &
      (
        y %in% c(
          "hfpef",
          "heart failure with preserved ejection fraction"
        ) |
        grepl(
          "hfpef|heart.?failure.?with.?preserved|preserved.?ejection",
          y
        )
      )
  ] <- "HFpEF"

  out[
    !is.na(y) &
      (
        y %in% c(
          "ctrl",
          "control",
          "controls",
          "normal",
          "non-failing",
          "non_failing",
          "nonfailing",
          "nf"
        ) |
        grepl(
          "^ctrl$|control|normal|non.?fail|^nf$|non.?hf",
          y
        )
      )
  ] <- "CTRL"

  out
}

# ------------------------------------------------------------------------------
# Gene annotation
# ------------------------------------------------------------------------------

build_gene_manifest <- function(VAR, n_genes) {

  if (nrow(VAR) != n_genes) {
    stop(
      "VAR row count does not equal expected gene count."
    )
  }

  index_col <- if (".index" %in% names(VAR)) {
    ".index"
  } else {
    names(VAR)[1]
  }

  gene_id <- norm_text(
    VAR[[index_col]]
  )

  symbol_candidates <- unique(
    c(
      intersect(
        c(
          "gene_symbol",
          "gene_name",
          "symbol",
          "feature_name",
          "name"
        ),
        names(VAR)
      ),
      names(VAR)[
        grepl(
          "symbol|gene.*name",
          names(VAR),
          ignore.case = TRUE
        )
      ]
    )
  )

  gene_symbol <- rep(
    NA_character_,
    n_genes
  )

  if (length(symbol_candidates)) {

    score <- vapply(
      symbol_candidates,
      function(nm) {

        x <- norm_text(
          VAR[[nm]]
        )

        valid <- !is.na(x) &
          nzchar(x)

        if (!any(valid)) {
          return(-Inf)
        }

        non_ensg <- mean(
          !grepl(
            "^ENSG[0-9]+",
            x[valid]
          )
        )

        unique_fraction <-
          data.table::uniqueN(
            x[valid]
          ) /
          sum(valid)

        non_ensg +
          0.2 *
          unique_fraction
      },
      numeric(1)
    )

    best <- symbol_candidates[
      which.max(score)
    ]

    gene_symbol <- norm_text(
      VAR[[best]]
    )
  }

  # If no useful symbol column is present, use index if it already looks like
  # gene symbols.
  if (
    all(
      is.na(gene_symbol) |
        !nzchar(gene_symbol)
    )
  ) {

    frac_ensg <- mean(
      grepl(
        "^ENSG[0-9]+",
        gene_id
      )
    )

    if (
      is.finite(frac_ensg) &&
      frac_ensg < 0.5
    ) {

      gene_symbol <- gene_id

    } else {

      ensembl_clean <- sub(
        "\\.[0-9]+$",
        "",
        gene_id
      )

      map <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(
          ensembl_clean
        ),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )

      gene_symbol <- unname(
        map[
          ensembl_clean
        ]
      )
    }
  }

  gene_key <- gene_id

  if (
    anyDuplicated(
      gene_key
    )
  ) {
    gene_key <- make.unique(
      gene_key
    )
  }

  data.table(
    gene_index = seq_len(n_genes),
    gene_key = gene_key,
    gene_id = gene_id,
    gene_symbol = gene_symbol
  )
}

# ------------------------------------------------------------------------------
# Read a row block from AnnData CSR group:
#   group/data
#   group/indices
#   group/indptr
#
# AnnData matrix orientation here is cells x genes.
# ------------------------------------------------------------------------------
read_csr_block <- function(
  h5file,
  group_path,
  indptr,
  row_start,
  row_end,
  n_genes
) {

  if (
    row_start < 1L ||
    row_end >
      (
        length(indptr) -
          1L
      ) ||
    row_start >
      row_end
  ) {
    stop("Invalid CSR row block.")
  }

  offset0 <- as.numeric(
    indptr[row_start]
  )

  offset1 <- as.numeric(
    indptr[row_end + 1L]
  )

  n_nz <- as.integer(
    offset1 -
      offset0
  )

  n_rows <- row_end -
    row_start +
    1L

  if (n_nz == 0L) {

    return(
      Matrix::Matrix(
        0,
        nrow = n_rows,
        ncol = n_genes,
        sparse = TRUE
      )
    )
  }

  vals <- rhdf5::h5read(
    h5file,
    paste0(
      group_path,
      "/data"
    ),
    start = offset0 + 1,
    count = n_nz
  )

  idx <- rhdf5::h5read(
    h5file,
    paste0(
      group_path,
      "/indices"
    ),
    start = offset0 + 1,
    count = n_nz
  )

  row_nnz <- diff(
    as.numeric(
      indptr[
        row_start:(
          row_end +
            1L
        )
      ]
    )
  )

  ii <- rep.int(
    seq_len(
      n_rows
    ),
    row_nnz
  )

  jj <- as.integer(
    idx
  ) +
    1L

  Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = as.numeric(vals),
    dims = c(
      n_rows,
      n_genes
    ),
    giveCsparse = TRUE
  )
}

# ------------------------------------------------------------------------------
# Equal-prevalence PPV for target-vs-nontarget cell classification.
# ------------------------------------------------------------------------------
equal_prevalence_ppv <- function(
  target_positive_fraction,
  nontarget_positive_fraction
) {

  den <- target_positive_fraction +
    nontarget_positive_fraction

  out <- rep(
    0,
    length(den)
  )

  ok <- is.finite(den) &
    den >
      0

  out[ok] <-
    target_positive_fraction[ok] /
      den[ok]

  out
}

# ------------------------------------------------------------------------------
# limma-voom model
# ------------------------------------------------------------------------------
run_voom_model <- function(
  counts_gene_by_donor,
  sample_meta,
  design_formula,
  coef_name
) {

  design <- model.matrix(
    design_formula,
    data = sample_meta
  )

  if (
    qr(design)$rank !=
      ncol(design)
  ) {
    stop(
      "Design matrix is not full rank for: ",
      paste(
        deparse(
          design_formula
        ),
        collapse = ""
      )
    )
  }

  if (!coef_name %in%
      colnames(design)) {
    stop(
      "Expected coefficient not found: ",
      coef_name
    )
  }

  y <- edgeR::DGEList(
    counts =
      counts_gene_by_donor
  )

  keep <- edgeR::filterByExpr(
    y,
    design = design
  )

  if (sum(keep) < 1000L) {
    stop(
      "Too few genes passed filterByExpr."
    )
  }

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )

  v <- limma::voom(
    y,
    design,
    plot = FALSE
  )

  fit <- limma::lmFit(
    v,
    design
  )

  fit <- limma::eBayes(
    fit,
    robust = TRUE
  )

  tt <- limma::topTable(
    fit,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )

  list(
    result =
      as.data.table(
        tt,
        keep.rownames = "gene_key"
      ),
    keep = keep,
    y = y,
    voom = v,
    design = design,
    fit = fit
  )
}

# ==============================================================================
# 4. HARD INPUT CHECK
# ==============================================================================

required_inputs <- c(
  H5AD_FILE,
  STEP13A_READINESS,
  COUNT_SOURCE_LOCK,
  META_FILE,
  VAR_FILE
)

if (!all(file.exists(required_inputs))) {
  stop(
    paste0(
      "Missing STEP13 input(s):\n",
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

READY <- fread(
  STEP13A_READINESS
)

if (!all(
  bool_pass(
    READY$pass
  )
)) {
  print(READY)
  stop(
    "STEP13A is not fully PASS."
  )
}

LOCK <- fread(
  COUNT_SOURCE_LOCK
)

if (
  nrow(LOCK) != 1L ||
  !bool_pass(
    LOCK$locked
  ) ||
  LOCK$selected_count_matrix_path[1] !=
    "layers/raw_count_cellranger"
) {
  stop(
    "STEP13A did not uniquely lock layers/raw_count_cellranger."
  )
}

COUNT_GROUP <- "layers/raw_count_cellranger"

# ==============================================================================
# 5. LOAD FROZEN METADATA / VAR
# ==============================================================================

META <- fread(
  META_FILE
)

VAR <- fread(
  VAR_FILE
)

if (
  nrow(META) !=
    EXPECTED_CELLS
) {
  stop(
    "Expected 48,866 metadata rows."
  )
}

if (
  nrow(VAR) !=
    EXPECTED_GENES
) {
  stop(
    "Expected 36,601 VAR rows."
  )
}

required_meta_cols <- c(
  "donor",
  "condition",
  "cell_type",
  "sex",
  "biosample_id"
)

missing_meta <- setdiff(
  required_meta_cols,
  names(META)
)

if (length(missing_meta)) {
  stop(
    "Locked metadata lacks: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

META[
  ,
  donor := norm_text(donor)
]

META[
  ,
  condition :=
    factor(
      norm_condition(
        condition
      ),
      levels = c(
        "CTRL",
        "HFpEF"
      )
    )
]

META[
  ,
  sex_norm :=
    norm_sex(
      sex
    )
]

META[
  ,
  cell_type :=
    norm_text(
      cell_type
    )
]

META[
  ,
  pool :=
    norm_text(
      biosample_id
    )
]

if (anyNA(
  META$condition
)) {

  CONDITION_AUDIT <- unique(
    data.table(
      condition_input =
        as.character(
          fread(
            META_FILE,
            select = "condition"
          )$condition
        ),
      condition_normalized =
        as.character(
          META$condition
        )
    )
  )

  fwrite(
    CONDITION_AUDIT,
    file.path(
      QC_DIR,
      "STEP13B_condition_normalization_audit.csv"
    )
  )

  stop(
    paste0(
      "Condition normalization failed. ",
      "Inspect STEP13B_condition_normalization_audit.csv."
    )
  )
}

if (anyNA(
  META$sex_norm
)) {
  stop(
    "Sex normalization failed."
  )
}

if (
  uniqueN(
    META$donor
  ) !=
    EXPECTED_DONORS
) {
  stop(
    "Expected 43 donors."
  )
}

GENE <- build_gene_manifest(
  VAR,
  EXPECTED_GENES
)

fwrite(
  GENE,
  file.path(
    QC_DIR,
    "STEP13B_gene_manifest.csv"
  )
)

# ==============================================================================
# 6. DONOR METADATA LOCK FOR CARDIOMYOCYTE ANALYSIS
# ==============================================================================

CM_META <- META[
  cell_type ==
    PRIMARY_CELLTYPE
]

if (
  uniqueN(
    CM_META$donor
  ) !=
    EXPECTED_DONORS
) {
  stop(
    "Cardiomyocyte is not represented in all 43 donors."
  )
}

DONOR_META_RAW <- unique(
  CM_META[
    ,
    .(
      donor,
      condition,
      sex_norm,
      pool
    )
  ]
)

if (
  nrow(
    DONOR_META_RAW
  ) !=
    EXPECTED_DONORS
) {
  stop(
    "At least one donor maps to multiple condition/sex/pool values."
  )
}

DONOR_META <- copy(
  DONOR_META_RAW
)

setorder(
  DONOR_META,
  donor
)

DONOR_META[
  ,
  condition :=
    factor(
      as.character(
        condition
      ),
      levels = c(
        "CTRL",
        "HFpEF"
      )
    )
]

DONOR_META[
  ,
  sex :=
    factor(
      sex_norm,
      levels = c(
        "F",
        "M"
      )
    )
]

DONOR_META[
  ,
  pool :=
    factor(
      pool
    )
]

DONOR_META[
  ,
  sex_norm := NULL
]

if (
  sum(
    DONOR_META$condition ==
      "HFpEF"
  ) !=
    EXPECTED_HFPEF
) {
  stop(
    "Expected 19 HFpEF donors."
  )
}

if (
  sum(
    DONOR_META$condition ==
      "CTRL"
  ) !=
    EXPECTED_CTRL
) {
  stop(
    "Expected 24 CTRL donors."
  )
}

CM_N <- CM_META[
  ,
  .(
    n_CM_nuclei = .N
  ),
  by = donor
]

DONOR_META <- merge(
  DONOR_META,
  CM_N,
  by = "donor",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  DONOR_META,
  donor
)

if (
  min(
    DONOR_META$n_CM_nuclei
  ) <
    20L
) {
  stop(
    "At least one donor has <20 cardiomyocyte nuclei."
  )
}

fwrite(
  DONOR_META,
  PB_META_FILE
)

# ==============================================================================
# 7. STREAM CELLRANGER RAW COUNTS -> CM PSEUDOBULK + BACKGROUND HEURISTIC
# ==============================================================================

if (
  file.exists(PB_FILE) &&
  file.exists(BKG_FILE) &&
  file.exists(STREAM_QC_FILE)
) {

  cat(
    "\nReusing STEP13B streaming checkpoints...\n"
  )

  PB <- readRDS(
    PB_FILE
  )

  BKG <- fread(
    BKG_FILE
  )

} else {

  cat(
    "\nStreaming CellRanger raw CSR matrix from H5AD...\n",
    "This pass builds CM patient pseudobulk AND the background heuristic.\n",
    sep = ""
  )

  indptr <- rhdf5::h5read(
    H5AD_FILE,
    paste0(
      COUNT_GROUP,
      "/indptr"
    )
  )

  indptr <- as.numeric(
    indptr
  )

  if (
    length(
      indptr
    ) !=
      EXPECTED_CELLS +
        1L
  ) {
    stop(
      "CellRanger CSR indptr length mismatch."
    )
  }

  donors <- DONOR_META$donor

  PB <- matrix(
    0,
    nrow = length(
      donors
    ),
    ncol = EXPECTED_GENES,
    dimnames = list(
      donors,
      GENE$gene_key
    )
  )

  global_umi <- numeric(
    EXPECTED_GENES
  )

  cm_pos0 <- numeric(
    EXPECTED_GENES
  )

  cm_pos1 <- numeric(
    EXPECTED_GENES
  )

  noncm_pos0 <- numeric(
    EXPECTED_GENES
  )

  noncm_pos1 <- numeric(
    EXPECTED_GENES
  )

  n_cm_total <- sum(
    META$cell_type ==
      PRIMARY_CELLTYPE
  )

  n_noncm_total <- EXPECTED_CELLS -
    n_cm_total

  block_starts <- seq(
    1L,
    EXPECTED_CELLS,
    by = BLOCK_SIZE
  )

  for (
    bb in seq_along(
      block_starts
    )
  ) {

    r1 <- block_starts[
      bb
    ]

    r2 <- min(
      EXPECTED_CELLS,
      r1 +
        BLOCK_SIZE -
        1L
    )

    B <- read_csr_block(
      h5file = H5AD_FILE,
      group_path = COUNT_GROUP,
      indptr = indptr,
      row_start = r1,
      row_end = r2,
      n_genes = EXPECTED_GENES
    )

    block_meta <- META[
      r1:r2
    ]

    # Global gene abundance.
    global_umi <- global_umi +
      as.numeric(
        Matrix::colSums(
          B
        )
      )

    all0 <- as.numeric(
      Matrix::colSums(
        B >
          0
      )
    )

    all1 <- as.numeric(
      Matrix::colSums(
        B >
          1
      )
    )

    cm_local <- which(
      block_meta$cell_type ==
        PRIMARY_CELLTYPE
    )

    if (length(
      cm_local
    )) {

      BCM <- B[
        cm_local,
        ,
        drop = FALSE
      ]

      cm0_block <- as.numeric(
        Matrix::colSums(
          BCM >
            0
        )
      )

      cm1_block <- as.numeric(
        Matrix::colSums(
          BCM >
            1
        )
      )

      cm_pos0 <- cm_pos0 +
        cm0_block

      cm_pos1 <- cm_pos1 +
        cm1_block

      noncm_pos0 <- noncm_pos0 +
        (
          all0 -
            cm0_block
        )

      noncm_pos1 <- noncm_pos1 +
        (
          all1 -
            cm1_block
        )

      donor_local <- match(
        block_meta$donor[
          cm_local
        ],
        donors
      )

      if (anyNA(
        donor_local
      )) {
        stop(
          "Could not match CM cell donor during H5AD streaming."
        )
      }

      M <- Matrix::sparseMatrix(
        i = donor_local,
        j = seq_along(
          donor_local
        ),
        x = 1,
        dims = c(
          length(
            donors
          ),
          length(
            donor_local
          )
        )
      )

      PB_add <- M %*%
        BCM

      PB <- PB +
        as.matrix(
          PB_add
        )

    } else {

      noncm_pos0 <- noncm_pos0 +
        all0

      noncm_pos1 <- noncm_pos1 +
        all1
    }

    if (
      bb %% 10L ==
        0L ||
      bb ==
        length(
          block_starts
        )
    ) {

      cat(
        "  block ",
        bb,
        "/",
        length(
          block_starts
        ),
        " completed\n",
        sep = ""
      )
    }

    rm(
      B
    )

    invisible(
      gc(
        verbose = FALSE
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Background heuristic
  # ---------------------------------------------------------------------------

  total_umi_all <- sum(
    global_umi
  )

  gene_umi_fraction <- global_umi /
    total_umi_all

  bkg_prob <- rank(
    gene_umi_fraction,
    ties.method = "average",
    na.last = "keep"
  ) /
    length(
      gene_umi_fraction
    )

  target_frac0 <- cm_pos0 /
    n_cm_total

  target_frac1 <- cm_pos1 /
    n_cm_total

  nontarget_frac0 <- noncm_pos0 /
    n_noncm_total

  nontarget_frac1 <- noncm_pos1 /
    n_noncm_total

  PPV0 <- equal_prevalence_ppv(
    target_frac0,
    nontarget_frac0
  )

  PPV1 <- equal_prevalence_ppv(
    target_frac1,
    nontarget_frac1
  )

  nontarget_prob <- 1 -
    rowMeans(
      cbind(
        PPV0,
        PPV1
      )
    )

  heuristic <- bkg_prob *
    nontarget_prob

  BKG <- copy(
    GENE
  )

  BKG[
    ,
    `:=`(
      global_UMI =
        global_umi,
      global_UMI_fraction =
        gene_umi_fraction,
      bkg_prob =
        bkg_prob,
      CM_expr_gt0_fraction =
        target_frac0,
      CM_expr_gt1_fraction =
        target_frac1,
      nonCM_expr_gt0_fraction =
        nontarget_frac0,
      nonCM_expr_gt1_fraction =
        nontarget_frac1,
      PPV0_equal_prevalence =
        PPV0,
      PPV1_equal_prevalence =
        PPV1,
      nontarget_prob =
        nontarget_prob,
      background_heuristic =
        heuristic,
      low_background =
        heuristic <=
          BACKGROUND_HEURISTIC_MAX
    )
  ]

  saveRDS(
    PB,
    PB_FILE,
    compress = TRUE
  )

  fwrite(
    BKG,
    BKG_FILE,
    compress = "gzip"
  )

  STREAM_QC <- data.table(
    metric = c(
      "n_cells_total",
      "n_CM_cells",
      "n_nonCM_cells",
      "n_genes",
      "n_donors",
      "CellRanger_total_UMI",
      "CM_pseudobulk_total_UMI",
      "genes_background_heuristic_le_0.40",
      "genes_background_heuristic_gt_0.40"
    ),
    value = c(
      EXPECTED_CELLS,
      n_cm_total,
      n_noncm_total,
      EXPECTED_GENES,
      EXPECTED_DONORS,
      total_umi_all,
      sum(
        PB
      ),
      sum(
        BKG$low_background,
        na.rm = TRUE
      ),
      sum(
        !BKG$low_background,
        na.rm = TRUE
      )
    )
  )

  fwrite(
    STREAM_QC,
    STREAM_QC_FILE
  )
}

# ==============================================================================
# 8. ORIENT PSEUDOBULK AS GENES x DONORS
# ==============================================================================

if (
  nrow(
    PB
  ) !=
    EXPECTED_DONORS ||
  ncol(
    PB
  ) !=
    EXPECTED_GENES
) {
  stop(
    "CM pseudobulk dimensions mismatch."
  )
}

if (!identical(
  rownames(
    PB
  ),
  DONOR_META$donor
)) {

  idx <- match(
    DONOR_META$donor,
    rownames(
      PB
    )
  )

  if (anyNA(
    idx
  )) {
    stop(
      "Could not align CM pseudobulk donors to metadata."
    )
  }

  PB <- PB[
    idx,
    ,
    drop = FALSE
  ]
}

COUNT_GD <- t(
  PB
)

rownames(
  COUNT_GD
) <- GENE$gene_key

colnames(
  COUNT_GD
) <- DONOR_META$donor

# ==============================================================================
# 9. PRIMARY LIMMA-VOOM MODEL: sex + pool + condition
# ==============================================================================

PRIMARY_FORMULA <- ~ sex + pool + condition

PRIMARY <- run_voom_model(
  counts_gene_by_donor = COUNT_GD,
  sample_meta = DONOR_META,
  design_formula = PRIMARY_FORMULA,
  coef_name = "conditionHFpEF"
)

RES <- PRIMARY$result

# Add gene annotation.
gidx <- match(
  RES$gene_key,
  GENE$gene_key
)

if (anyNA(
  gidx
)) {
  stop(
    "Primary DE results cannot be mapped to gene manifest."
  )
}

RES[
  ,
  `:=`(
    gene_id =
      GENE$gene_id[
        gidx
      ],
    gene_symbol =
      GENE$gene_symbol[
        gidx
      ]
  )
]

# Add background heuristic.
bidx <- match(
  RES$gene_key,
  BKG$gene_key
)

if (anyNA(
  bidx
)) {
  stop(
    "Primary DE results cannot be mapped to background heuristic."
  )
}

RES[
  ,
  background_heuristic :=
    BKG$background_heuristic[
      bidx
    ]
]

RES[
  ,
  low_background :=
    background_heuristic <=
      BACKGROUND_HEURISTIC_MAX
]

RES[
  ,
  FDR05 :=
    adj.P.Val <
      DE_FDR
]

RES[
  ,
  FDR05_low_background :=
    adj.P.Val <
      DE_FDR &
    low_background
]

RES[
  ,
  direction :=
    fifelse(
      logFC >
        0,
      "HFpEF_up",
      fifelse(
        logFC <
          0,
        "HFpEF_down",
        "No_change"
      )
    )
]

setorder(
  RES,
  adj.P.Val,
  P.Value
)

fwrite(
  RES,
  ALL_DE_FILE,
  compress = "gzip"
)

fwrite(
  RES[
    FDR05 ==
      TRUE
  ],
  FDR_DE_FILE
)

fwrite(
  RES[
    FDR05_low_background ==
      TRUE
  ],
  FDR_BKG_DE_FILE
)

# ==============================================================================
# 10. SENSITIVITY MODEL: sex + condition
# ==============================================================================

SENS_FORMULA <- ~ sex + condition

SENS <- run_voom_model(
  counts_gene_by_donor = COUNT_GD,
  sample_meta = DONOR_META,
  design_formula = SENS_FORMULA,
  coef_name = "conditionHFpEF"
)

SRES <- SENS$result

sidx <- match(
  SRES$gene_key,
  GENE$gene_key
)

SRES[
  ,
  `:=`(
    gene_id =
      GENE$gene_id[
        sidx
      ],
    gene_symbol =
      GENE$gene_symbol[
        sidx
      ]
  )
]

fwrite(
  SRES,
  SENS_DE_FILE,
  compress = "gzip"
)

COMMON <- merge(
  RES[
    ,
    .(
      gene_key,
      primary_logFC =
        logFC,
      primary_P =
        P.Value,
      primary_FDR =
        adj.P.Val
    )
  ],
  SRES[
    ,
    .(
      gene_key,
      sensitivity_logFC =
        logFC,
      sensitivity_P =
        P.Value,
      sensitivity_FDR =
        adj.P.Val
    )
  ],
  by = "gene_key"
)

LOGFC_SPEARMAN <- suppressWarnings(
  cor(
    COMMON$primary_logFC,
    COMMON$sensitivity_logFC,
    method = "spearman",
    use = "complete.obs"
  )
)

DIRECTION_AGREEMENT <- mean(
  sign(
    COMMON$primary_logFC
  ) ==
    sign(
      COMMON$sensitivity_logFC
    )
)

fwrite(
  COMMON,
  file.path(
    SENS_DIR,
    "STEP13B_primary_vs_sensitivity_gene_comparison.csv.gz"
  ),
  compress = "gzip"
)

# ==============================================================================
# 11. PUBLICATION BENCHMARK — DESCRIPTIVE ONLY
# ==============================================================================

BENCHMARK <- data.table(
  gene_symbol = c(
    "FKBP5",
    "MTFR1",
    "HMOX2",
    "NMRAL1",
    "EGLN1",
    "PLPP3",
    "MYOF",
    "PGR",
    "HSD11B1",
    "MID1",
    "EFCAB11",
    "ZMAT1",
    "ZNF20",
    "ZFP2",
    "SCN2B"
  ),
  published_direction = c(
    "down",
    "down",
    "down",
    "down",
    "down",
    "down",
    "up",
    "up",
    "up",
    "up",
    "up",
    "up",
    "up",
    "up",
    "up"
  )
)

BENCH <- merge(
  BENCHMARK,
  RES[
    !is.na(
      gene_symbol
    ),
    .(
      gene_symbol,
      reanalysis_logFC =
        logFC,
      reanalysis_P =
        P.Value,
      reanalysis_FDR =
        adj.P.Val,
      background_heuristic,
      FDR05_low_background
    )
  ],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

BENCH[
  ,
  reanalysis_direction :=
    fifelse(
      reanalysis_logFC >
        0,
      "up",
      fifelse(
        reanalysis_logFC <
          0,
        "down",
        NA_character_
      )
    )
]

BENCH[
  ,
  direction_concordant :=
    reanalysis_direction ==
      published_direction
]

fwrite(
  BENCH,
  file.path(
    QC_DIR,
    "STEP13B_published_CM_gene_direction_benchmark.csv"
  )
)

# ==============================================================================
# 12. SAMPLE / MODEL QC
# ==============================================================================

SAMPLE_QC <- copy(
  DONOR_META
)

SAMPLE_QC[
  ,
  raw_library_size :=
    as.numeric(
      colSums(
        COUNT_GD
      )
    )
]

SAMPLE_QC[
  ,
  TMM_norm_factor :=
    PRIMARY$y$samples$norm.factors
]

SAMPLE_QC[
  ,
  effective_library_size :=
    PRIMARY$y$samples$lib.size *
      PRIMARY$y$samples$norm.factors
]

fwrite(
  SAMPLE_QC,
  file.path(
    QC_DIR,
    "STEP13B_Cardiomyocyte_sample_QC.csv"
  )
)

DESIGN_AUDIT <- data.table(
  model = c(
    "primary",
    "sensitivity"
  ),
  formula = c(
    "~ sex + pool + condition",
    "~ sex + condition"
  ),
  n_samples = c(
    nrow(
      DONOR_META
    ),
    nrow(
      DONOR_META
    )
  ),
  n_coefficients = c(
    ncol(
      PRIMARY$design
    ),
    ncol(
      SENS$design
    )
  ),
  rank = c(
    qr(
      PRIMARY$design
    )$rank,
    qr(
      SENS$design
    )$rank
  ),
  full_rank = c(
    qr(
      PRIMARY$design
    )$rank ==
      ncol(
        PRIMARY$design
      ),
    qr(
      SENS$design
    )$rank ==
      ncol(
        SENS$design
      )
  )
)

fwrite(
  DESIGN_AUDIT,
  file.path(
    QC_DIR,
    "STEP13B_design_matrix_QC.csv"
  )
)

SUMMARY <- data.table(
  metric = c(
    "CM_donors",
    "HFpEF_donors",
    "CTRL_donors",
    "CM_nuclei_total",
    "CM_nuclei_min_per_donor",
    "CM_nuclei_median_per_donor",
    "genes_input",
    "genes_tested_primary",
    "FDR05_primary",
    "FDR05_low_background",
    "HFpEF_up_FDR05_low_background",
    "HFpEF_down_FDR05_low_background",
    "primary_vs_sensitivity_logFC_Spearman",
    "primary_vs_sensitivity_direction_agreement",
    "published_benchmark_direction_concordance"
  ),
  value = c(
    nrow(
      DONOR_META
    ),
    sum(
      DONOR_META$condition ==
        "HFpEF"
    ),
    sum(
      DONOR_META$condition ==
        "CTRL"
    ),
    sum(
      DONOR_META$n_CM_nuclei
    ),
    min(
      DONOR_META$n_CM_nuclei
    ),
    median(
      DONOR_META$n_CM_nuclei
    ),
    nrow(
      COUNT_GD
    ),
    nrow(
      RES
    ),
    sum(
      RES$FDR05
    ),
    sum(
      RES$FDR05_low_background
    ),
    sum(
      RES$FDR05_low_background &
        RES$logFC >
          0
    ),
    sum(
      RES$FDR05_low_background &
        RES$logFC <
          0
    ),
    LOGFC_SPEARMAN,
    DIRECTION_AGREEMENT,
    mean(
      BENCH$direction_concordant,
      na.rm = TRUE
    )
  )
)

fwrite(
  SUMMARY,
  file.path(
    QC_DIR,
    "STEP13B_DE_summary.csv"
  )
)

# ==============================================================================
# 13. FIGURES — SINGLE PANEL, NO LEGEND INSIDE
# ==============================================================================

# Volcano uses raw P on y-axis, consistent with the SCP3342 publication figure.
PLOT <- copy(
  RES
)

PLOT[
  ,
  plot_group :=
    fifelse(
      FDR05_low_background &
        logFC >
          0,
      "HFpEF_up",
      fifelse(
        FDR05_low_background &
          logFC <
            0,
        "HFpEF_down",
        "Other"
      )
    )
]

LABEL_UP <- PLOT[
  FDR05_low_background &
    logFC >
      0 &
    !is.na(
      gene_symbol
    ) &
    nzchar(
      gene_symbol
    )
][
  order(
    P.Value
  )
][
  1:min(
    .N,
    10L
  )
]

LABEL_DOWN <- PLOT[
  FDR05_low_background &
    logFC <
      0 &
    !is.na(
      gene_symbol
    ) &
    nzchar(
      gene_symbol
    )
][
  order(
    P.Value
  )
][
  1:min(
    .N,
    10L
  )
]

LABEL <- rbindlist(
  list(
    LABEL_UP,
    LABEL_DOWN
  ),
  fill = TRUE
)

p_vol <- ggplot(
  PLOT,
  aes(
    x = logFC,
    y = -log10(
      pmax(
        P.Value,
        1e-300
      )
    ),
    colour = plot_group
  )
) +
  geom_point(
    size = 0.55,
    alpha = 0.58
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_hline(
    yintercept =
      -log10(
        DE_FDR
      ),
    linetype = 2,
    linewidth = 0.45,
    colour = "grey55"
  ) +
  scale_colour_manual(
    values = c(
      "Other" = "grey78",
      "HFpEF_up" = "#B34D4D",
      "HFpEF_down" = "#557DA1"
    )
  ) +
  ggrepel::geom_text_repel(
    data = LABEL,
    aes(
      label = gene_symbol
    ),
    size = 2.7,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  labs(
    x = "log2 fold change (HFpEF vs CTRL)",
    y = expression(
      -log[10](P)
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
      colour = "black",
      linewidth = 0.55
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
    "Figure_STEP13B_CM_HFpEF_vs_CTRL_volcano.tiff"
  ),
  p_vol,
  width = 5.8,
  height = 5.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13B_CM_HFpEF_vs_CTRL_volcano.pdf"
  ),
  p_vol,
  width = 5.8,
  height = 5.0,
  units = "in"
)

# MDS from primary voom/DGE object.
mds <- limma::plotMDS(
  PRIMARY$y,
  top = 500,
  plot = FALSE
)

MDS_DT <- data.table(
  donor =
    DONOR_META$donor,
  condition =
    as.character(
      DONOR_META$condition
    ),
  sex =
    as.character(
      DONOR_META$sex
    ),
  pool =
    as.character(
      DONOR_META$pool
    ),
  MDS1 =
    mds$x,
  MDS2 =
    mds$y
)

fwrite(
  MDS_DT,
  file.path(
    SOURCE_DIR,
    "Figure_STEP13B_CM_MDS_source_data.csv"
  )
)

p_mds <- ggplot(
  MDS_DT,
  aes(
    x = MDS1,
    y = MDS2,
    shape = condition
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  labs(
    x = "MDS1",
    y = "MDS2"
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
    "Figure_STEP13B_CM_MDS_QC.tiff"
  ),
  p_mds,
  width = 5.2,
  height = 4.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================================================================
# 14. READINESS
# ==============================================================================

READINESS <- data.table(
  check = c(
    "STEP13A_fully_passed",
    "CellRanger_raw_count_source_locked",
    "Cardiomyocyte_present_in_43_donors",
    "19_HFpEF_donors",
    "24_CTRL_donors",
    "all_CM_donors_have_ge20_nuclei",
    "sex_available_for_all_donors",
    "pool_available_for_all_donors",
    "primary_design_sex_pool_condition_full_rank",
    "sensitivity_design_sex_condition_full_rank",
    "genes_tested_primary_gt1000",
    "all_primary_Pvalues_finite",
    "all_primary_FDR_finite",
    "background_heuristic_computed",
    "primary_DE_result_written",
    "primary_vs_sensitivity_comparison_written",
    "STEP13B_complete"
  ),
  pass = c(
    all(
      bool_pass(
        READY$pass
      )
    ),
    LOCK$selected_count_matrix_path[1] ==
      "layers/raw_count_cellranger",
    nrow(
      DONOR_META
    ) ==
      EXPECTED_DONORS,
    sum(
      DONOR_META$condition ==
        "HFpEF"
    ) ==
      EXPECTED_HFPEF,
    sum(
      DONOR_META$condition ==
        "CTRL"
    ) ==
      EXPECTED_CTRL,
    min(
      DONOR_META$n_CM_nuclei
    ) >=
      20L,
    all(
      !is.na(
        DONOR_META$sex
      )
    ),
    all(
      !is.na(
        DONOR_META$pool
      )
    ),
    qr(
      PRIMARY$design
    )$rank ==
      ncol(
        PRIMARY$design
      ),
    qr(
      SENS$design
    )$rank ==
      ncol(
        SENS$design
      ),
    nrow(
      RES
    ) >
      1000L,
    all(
      is.finite(
        RES$P.Value
      )
    ),
    all(
      is.finite(
        RES$adj.P.Val
      )
    ),
    nrow(
      BKG
    ) ==
      EXPECTED_GENES &&
      all(
        is.finite(
          BKG$background_heuristic
        )
      ),
    file.exists(
      ALL_DE_FILE
    ),
    file.exists(
      file.path(
        SENS_DIR,
        "STEP13B_primary_vs_sensitivity_gene_comparison.csv.gz"
      )
    ),
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP13B_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "dataset",
    "tissue",
    "primary_cell_type",
    "statistical_unit",
    "raw_count_source",
    "pseudobulk_definition",
    "DE_method",
    "primary_model",
    "sensitivity_model",
    "primary_contrast",
    "multiple_testing",
    "background_heuristic_threshold",
    "CellBender_raw_count_status",
    "candidate_gene_filtering_STEP13B",
    "next_step"
  ),
  value = c(
    "SCP3342",
    "Human ventricular septal myocardium",
    "Cardiomyocyte",
    "Biological donor/patient",
    "layers/raw_count_cellranger",
    "CellRanger raw counts summed within donor x Cardiomyocyte",
    "limma-voom",
    "~ sex + pool + condition",
    "~ sex + condition",
    "HFpEF vs CTRL",
    "Benjamini-Hochberg within Cardiomyocyte",
    as.character(
      BACKGROUND_HEURISTIC_MAX
    ),
    "Released H5AD layer is continuous/non-integer and is not summed as raw counts; no false dual-matrix reproduction claim",
    "None; transcriptome-wide analysis blind to upstream candidate genes",
    "STEP13C cross-disease integration: frozen genetics + AF-CM + HFpEF-CM"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP13B_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP13B_sessionInfo.txt"
  )
)

# ==============================================================================
# 15. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP13B COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

print(
  SUMMARY
)

cat(
  "\nPublished cardiomyocyte benchmark directions:\n"
)
print(
  BENCH
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
    "\nSTEP13B PASSED.\n",
    "Do not select new genetics candidates from the HFpEF transcriptome.\n",
    "Next: STEP13C will intersect the already-frozen genetics/regulatory ",
    "candidate set with frozen AF-CM and HFpEF-CM disease-tissue results.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP13B failed one or more technical checks.\n",
    "Do not start STEP13C yet.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  DE_DIR,
  "\n",
  SENS_DIR,
  "\n",
  sep = ""
)
