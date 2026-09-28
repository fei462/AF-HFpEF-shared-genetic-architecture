# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP12B_V2_FIX_LIMMA_MDS_GSE255612_PATIENT_LEVEL_DE(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP12B V2 — GSE255612 AF vs CTRL PATIENT-LEVEL PSEUDOBULK DE
#             PRIMARY: CARDIOMYOCYTES
#             SECONDARY: MACROPHAGES
#
# Project: AF–HFpEF–BMI–OSA shared genetics and atrial cellular architecture
#
# STUDY-ROLE OF THIS STEP
#   This is the disease-tissue validation layer.
#
#   Upstream evidence has already been frozen:
#     GWAS / cross-trait LDSC
#       -> conjFDR + LAVA
#       -> MR/MVMR
#       -> SMR/HEIDI + sQTL-SMR/HEIDI + coloc/SuSiE
#       -> GSE238242 S-LDSC
#       -> AF-priority SCAVENGE
#       -> Cardiomyocytes (CM) prioritized
#
#   STEP12B asks an independent question:
#     In real AF left atrial tissue, do cardiomyocytes show reproducible
#     patient-level transcriptional alterations?
#
# IMPORTANT ANTI-CIRCULARITY RULE
#   STEP12B is intentionally BLIND to the upstream candidate-gene list.
#   We do NOT select genes because they were significant upstream.
#   We first freeze all patient-level CM DE results.
#   STEP12C will then intersect the frozen STEP12B results with the already-
#   frozen genetics/regulatory candidate genes.
#
# PRIMARY ANALYSIS
#   Author-defined Cardiomyocytes
#   Statistical unit = biological donor/patient
#   n = 34 donors (18 AF, 16 CTRL)
#   edgeR quasi-likelihood pseudobulk model
#   primary design = ~ sex + condition, if sex is identifiable and full rank
#   fallback only if sex unavailable/non-estimable: ~ condition
#
# SECONDARY CONTEXT ANALYSIS
#   Author-defined Macrophages
#   Same patient-level pipeline
#   This is secondary context, NOT allowed to replace the CM primary analysis.
#
# NO CELL-LEVEL WILCOXON TESTS.
# NO RECLUSTERING.
# NO REANNOTATION.
# NO CANDIDATE-GENE-DRIVEN FILTERING.
#
# OUTPUT
#   D:/A/data/STEP12_GSE255612_AF_snRNA/STEP12B_AF_vs_CTRL_DE/
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = max(3600, getOption("timeout")))
set.seed(9527)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

STEP12_ROOT <- file.path(
  DATA_ROOT,
  "STEP12_GSE255612_AF_snRNA"
)

STEP12A_QC <- file.path(
  STEP12_ROOT,
  "00_QC",
  "STEP12A_readiness.csv"
)

PB_COUNTS_FILE <- file.path(
  STEP12_ROOT,
  "02_PSEUDOBULK",
  "GSE255612_patient_by_celltype_pseudobulk_counts.rds"
)

PB_META_FILE <- file.path(
  STEP12_ROOT,
  "02_PSEUDOBULK",
  "GSE255612_patient_by_celltype_pseudobulk_metadata.csv"
)

META_LOCK_FILE <- file.path(
  STEP12_ROOT,
  "01_OBJECTS",
  "GSE255612_author_metadata_LOCKED.csv.gz"
)

GENE_MANIFEST_FILE <- file.path(
  STEP12_ROOT,
  "01_OBJECTS",
  "GSE255612_gene_manifest_LOCKED.csv"
)

PRIMARY_CELLTYPE <- "Cardiomyocytes"
SECONDARY_CELLTYPE <- "Macrophages"

EXPECTED_DONORS <- 34L
EXPECTED_AF <- 18L
EXPECTED_CTRL <- 16L

# Gene-level interpretation thresholds.
# Primary genome-wide DE significance is FDR < 0.05.
# We additionally export a larger-effect subset at |log2FC| >= 0.5,
# but do NOT require an FC threshold for downstream candidate validation.
DE_FDR <- 0.05
LARGE_EFFECT_ABS_LOGFC <- 0.5

# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

OUT_ROOT <- file.path(
  STEP12_ROOT,
  "STEP12B_AF_vs_CTRL_DE"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
PRIMARY_DIR <- file.path(OUT_ROOT, "01_PRIMARY_CM")
SECONDARY_DIR <- file.path(OUT_ROOT, "02_SECONDARY_MACROPHAGE")
SOURCE_DIR <- file.path(OUT_ROOT, "03_SOURCE_DATA")
FIG_DIR <- file.path(OUT_ROOT, "04_FIGURES")

for (d in c(
  OUT_ROOT,
  QC_DIR,
  PRIMARY_DIR,
  SECONDARY_DIR,
  SOURCE_DIR,
  FIG_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

cran_pkgs <- c(
  "data.table",
  "ggplot2",
  "ggrepel",
  "Matrix"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Could not install/load CRAN package: ", p)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("edgeR", "limma")

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
    "Required Bioconductor package(s) unavailable: ",
    paste(still_missing_bioc, collapse = ", ")
  )
}

library(data.table)
library(ggplot2)
library(ggrepel)
library(Matrix)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

bool_pass <- function(x) {
  x %in% c(TRUE, "TRUE", 1, "1")
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_sex <- function(x) {

  y <- tolower(trimws(as.character(x)))

  out <- rep(NA_character_, length(y))

  out[
    y %in% c(
      "f", "female", "woman", "women"
    )
  ] <- "F"

  out[
    y %in% c(
      "m", "male", "man", "men"
    )
  ] <- "M"

  out
}

make_gene_annotation <- function(gene_manifest, count_rownames) {

  gm <- copy(gene_manifest)

  required <- c(
    "gene_key",
    "gene_id",
    "gene_id_clean",
    "gene_symbol"
  )

  missing <- setdiff(required, names(gm))

  if (length(missing)) {
    stop(
      "Gene manifest lacks columns: ",
      paste(missing, collapse = ", ")
    )
  }

  idx <- match(
    count_rownames,
    gm$gene_key
  )

  if (anyNA(idx)) {
    stop(
      "Could not map all pseudobulk matrix rownames to frozen gene manifest."
    )
  }

  gm[idx]
}

# ------------------------------------------------------------------------------
# Patient-level pseudobulk DE function
# ------------------------------------------------------------------------------

run_pseudobulk_DE <- function(
  target_celltype,
  analysis_label,
  out_dir,
  pb_counts,
  pb_meta,
  donor_covariates,
  gene_anno
) {

  cat(
    "\n============================================================\n",
    analysis_label,
    " : ",
    target_celltype,
    "\n============================================================\n",
    sep = ""
  )

  META <- copy(
    pb_meta[
      cell_type == target_celltype &
        eligible_primary_DE == TRUE
    ]
  )

  if (!nrow(META)) {
    stop("No eligible pseudobulk samples for ", target_celltype)
  }

  # Match pseudobulk matrix columns.
  idx_col <- match(
    META$pb_id,
    colnames(pb_counts)
  )

  if (anyNA(idx_col)) {
    stop(
      "Could not align ",
      target_celltype,
      " metadata to pseudobulk count matrix."
    )
  }

  COUNT <- pb_counts[
    ,
    idx_col,
    drop = FALSE
  ]

  colnames(COUNT) <- META$pb_id

  # Add donor-level covariates.
  idx_cov <- match(
    META$donor,
    donor_covariates$donor
  )

  if (anyNA(idx_cov)) {
    stop("Could not map donor covariates for ", target_celltype)
  }

  META[
    ,
    sex := donor_covariates$sex[idx_cov]
  ]

  # Freeze condition levels: CTRL is reference, AF is tested.
  META[
    ,
    condition :=
      factor(
        condition,
        levels = c("CTRL", "AF")
      )
  ]

  # Donor/sample counts.
  if (uniqueN(META$donor) != EXPECTED_DONORS) {
    stop(
      target_celltype,
      ": expected 34 eligible donors; observed ",
      uniqueN(META$donor)
    )
  }

  if (sum(META$condition == "AF") != EXPECTED_AF) {
    stop(
      target_celltype,
      ": expected 18 AF donors."
    )
  }

  if (sum(META$condition == "CTRL") != EXPECTED_CTRL) {
    stop(
      target_celltype,
      ": expected 16 CTRL donors."
    )
  }

  # ---------------------------------------------------------------------------
  # Decide whether sex can be used as a covariate.
  # We prefer ~ sex + condition.
  # It is used only if:
  #   - no missing sex,
  #   - two sex levels,
  #   - each condition contains both sexes,
  #   - model matrix is full rank.
  # ---------------------------------------------------------------------------

  use_sex <- FALSE

  if (
    all(!is.na(META$sex)) &&
    uniqueN(META$sex) == 2L
  ) {

    META[
      ,
      sex := factor(
        sex,
        levels = c("F", "M")
      )
    ]

    SEX_TAB <- table(
      META$condition,
      META$sex
    )

    candidate_design <- model.matrix(
      ~ sex + condition,
      data = META
    )

    if (
      all(SEX_TAB > 0) &&
      qr(candidate_design)$rank ==
        ncol(candidate_design)
    ) {
      use_sex <- TRUE
    }
  }

  if (use_sex) {

    DESIGN <- model.matrix(
      ~ sex + condition,
      data = META
    )

    DESIGN_FORMULA <- "~ sex + condition"

  } else {

    DESIGN <- model.matrix(
      ~ condition,
      data = META
    )

    DESIGN_FORMULA <- "~ condition"
  }

  if (
    qr(DESIGN)$rank !=
      ncol(DESIGN)
  ) {
    stop(
      target_celltype,
      ": DE design matrix is not full rank."
    )
  }

  if (!"conditionAF" %in% colnames(DESIGN)) {
    stop(
      target_celltype,
      ": expected conditionAF coefficient not found."
    )
  }

  # ---------------------------------------------------------------------------
  # edgeR quasi-likelihood DE
  # ---------------------------------------------------------------------------

  Y <- edgeR::DGEList(
    counts = COUNT,
    samples = as.data.frame(META)
  )

  KEEP <- edgeR::filterByExpr(
    Y,
    design = DESIGN
  )

  if (sum(KEEP) < 1000L) {
    stop(
      target_celltype,
      ": fewer than 1000 genes passed filterByExpr."
    )
  }

  YF <- Y[
    KEEP,
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
    coef = "conditionAF"
  )

  TAB <- edgeR::topTags(
    QLF,
    n = Inf,
    sort.by = "PValue"
  )$table

  RES <- as.data.table(
    TAB,
    keep.rownames = "gene_key"
  )

  idx_gene <- match(
    RES$gene_key,
    gene_anno$gene_key
  )

  if (anyNA(idx_gene)) {
    stop(
      target_celltype,
      ": result genes cannot be mapped to frozen gene manifest."
    )
  }

  RES[
    ,
    `:=`(
      gene_id =
        gene_anno$gene_id[idx_gene],
      gene_id_clean =
        gene_anno$gene_id_clean[idx_gene],
      gene_symbol =
        gene_anno$gene_symbol[idx_gene]
    )
  ]

  RES[
    ,
    direction :=
      fifelse(
        logFC > 0,
        "AF_up",
        fifelse(
          logFC < 0,
          "AF_down",
          "No_change"
        )
      )
  ]

  RES[
    ,
    FDR05 :=
      FDR < DE_FDR
  ]

  RES[
    ,
    FDR05_large_effect :=
      FDR < DE_FDR &
      abs(logFC) >=
        LARGE_EFFECT_ABS_LOGFC
  ]

  # Preserve ranked output.
  setorder(
    RES,
    FDR,
    PValue
  )

  # ---------------------------------------------------------------------------
  # Donor-level normalized expression for downstream STEP12C.
  # This allows plotting only the PRE-SPECIFIED genes after integration.
  # ---------------------------------------------------------------------------

  LOGCPM <- edgeR::cpm(
    YF,
    log = TRUE,
    prior.count = 2
  )

  LOGCPM_DT <- as.data.table(
    LOGCPM,
    keep.rownames = "gene_key"
  )

  # ---------------------------------------------------------------------------
  # Sample QC
  # ---------------------------------------------------------------------------

  SAMPLE_QC <- copy(META)

  SAMPLE_QC[
    ,
    raw_library_size :=
      as.numeric(
        Matrix::colSums(
          COUNT
        )
      )
  ]

  SAMPLE_QC[
    ,
    TMM_norm_factor :=
      YF$samples$norm.factors
  ]

  SAMPLE_QC[
    ,
    effective_library_size :=
      YF$samples$lib.size *
        YF$samples$norm.factors
  ]

  # ---------------------------------------------------------------------------
  # MDS using the limma plotMDS generic with the DGEList method.
  # Only for QC; no clustering is used for inference.
  # ---------------------------------------------------------------------------

  mds <- limma::plotMDS(
    YF,
    top = 500,
    plot = FALSE
  )

  MDS_DT <- data.table(
    pb_id = META$pb_id,
    donor = META$donor,
    condition = as.character(META$condition),
    sex = as.character(META$sex),
    MDS1 = mds$x,
    MDS2 = mds$y
  )

  # ---------------------------------------------------------------------------
  # DE summary
  # ---------------------------------------------------------------------------

  N_FDR <- sum(
    RES$FDR05
  )

  N_UP <- sum(
    RES$FDR05 &
      RES$logFC > 0
  )

  N_DOWN <- sum(
    RES$FDR05 &
      RES$logFC < 0
  )

  N_LARGE <- sum(
    RES$FDR05_large_effect
  )

  SUMMARY <- data.table(
    metric = c(
      "cell_type",
      "analysis_role",
      "eligible_donors",
      "AF_donors",
      "CTRL_donors",
      "design_formula",
      "sex_covariate_used",
      "genes_input",
      "genes_tested",
      "FDR05_DEGs",
      "FDR05_AF_up",
      "FDR05_AF_down",
      "FDR05_abslogFC_ge_0.5",
      "common_dispersion"
    ),
    value = c(
      target_celltype,
      analysis_label,
      uniqueN(META$donor),
      sum(META$condition == "AF"),
      sum(META$condition == "CTRL"),
      DESIGN_FORMULA,
      use_sex,
      nrow(COUNT),
      nrow(RES),
      N_FDR,
      N_UP,
      N_DOWN,
      N_LARGE,
      YF$common.dispersion
    )
  )

  # ---------------------------------------------------------------------------
  # Write data
  # ---------------------------------------------------------------------------

  prefix <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    target_celltype
  )

  fwrite(
    RES,
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_AF_vs_CTRL_all_genes.csv.gz"
      )
    ),
    compress = "gzip"
  )

  fwrite(
    RES[
      FDR05 == TRUE
    ],
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_AF_vs_CTRL_FDR05_DEG.csv"
      )
    )
  )

  fwrite(
    RES[
      FDR05_large_effect == TRUE
    ],
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_AF_vs_CTRL_FDR05_abslogFC05_DEG.csv"
      )
    )
  )

  fwrite(
    SAMPLE_QC,
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_sample_QC.csv"
      )
    )
  )

  fwrite(
    MDS_DT,
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_MDS_source_data.csv"
      )
    )
  )

  fwrite(
    SUMMARY,
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_DE_summary.csv"
      )
    )
  )

  saveRDS(
    LOGCPM,
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_filtered_logCPM_patient_matrix.rds"
      )
    ),
    compress = TRUE
  )

  fwrite(
    data.table(
      gene_key = rownames(LOGCPM)
    ),
    file.path(
      out_dir,
      paste0(
        "STEP12B_",
        prefix,
        "_filtered_logCPM_gene_order.csv"
      )
    )
  )

  # ---------------------------------------------------------------------------
  # Figures — single figure per file, no legend inside
  # ---------------------------------------------------------------------------

  # Volcano:
  # grey = not FDR significant
  # red = AF-up FDR significant
  # blue = AF-down FDR significant

  PLOT <- copy(
    RES
  )

  PLOT[
    ,
    volcano_group :=
      fifelse(
        FDR < DE_FDR &
          logFC > 0,
        "AF_up",
        fifelse(
          FDR < DE_FDR &
            logFC < 0,
          "AF_down",
          "NS"
        )
      )
  ]

  # Label only strongest transcriptome-wide DE genes.
  LABEL <- PLOT[
    FDR < DE_FDR &
      !is.na(gene_symbol) &
      nzchar(gene_symbol)
  ][
    order(
      FDR,
      -abs(logFC)
    )
  ][
    1:min(.N, 15L)
  ]

  p_vol <- ggplot(
    PLOT,
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
      alpha = 0.62
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
      colour = "grey45"
    ) +
    scale_colour_manual(
      values = c(
        "NS" = "grey78",
        "AF_up" = "#B34D4D",
        "AF_down" = "#557DA1"
      )
    ) +
    ggrepel::geom_text_repel(
      data = LABEL,
      aes(
        label = gene_symbol
      ),
      size = 2.7,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    labs(
      x = "log2 fold change (AF vs CTRL)",
      y = expression(
        -log[10](FDR)
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
        8, 10, 8, 8
      )
    )

  ggsave(
    file.path(
      FIG_DIR,
      paste0(
        "Figure_STEP12B_",
        prefix,
        "_AF_vs_CTRL_volcano.tiff"
      )
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
      paste0(
        "Figure_STEP12B_",
        prefix,
        "_AF_vs_CTRL_volcano.pdf"
      )
    ),
    p_vol,
    width = 5.8,
    height = 5.0,
    units = "in"
  )

  # MDS.
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
      paste0(
        "Figure_STEP12B_",
        prefix,
        "_MDS_QC.tiff"
      )
    ),
    p_mds,
    width = 5.2,
    height = 4.6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  # Return frozen objects.
  list(
    result = RES,
    summary = SUMMARY,
    sample_qc = SAMPLE_QC,
    logCPM = LOGCPM,
    design = DESIGN,
    design_formula = DESIGN_FORMULA,
    sex_used = use_sex,
    keep = KEEP
  )
}

# ==============================================================================
# 4. HARD INPUT / STEP12A CHECK
# ==============================================================================

required_inputs <- c(
  STEP12A_QC,
  PB_COUNTS_FILE,
  PB_META_FILE,
  META_LOCK_FILE,
  GENE_MANIFEST_FILE
)

if (!all(file.exists(required_inputs))) {
  stop(
    paste0(
      "Missing STEP12A input(s):\n",
      paste(
        required_inputs[
          !file.exists(required_inputs)
        ],
        collapse = "\n"
      )
    )
  )
}

READY12A <- fread(
  STEP12A_QC
)

if (!all(
  bool_pass(
    READY12A$pass
  )
)) {
  print(READY12A)
  stop(
    "STEP12A is not fully PASS. Do not run STEP12B."
  )
}

# ==============================================================================
# 5. LOAD PSEUDOBULK + FROZEN AUTHOR METADATA
# ==============================================================================

cat("\nLoading patient×cell-type pseudobulk counts...\n")

PB_COUNTS <- readRDS(
  PB_COUNTS_FILE
)

PB_COUNTS <- methods::as(
  PB_COUNTS,
  "dgCMatrix"
)

PB_META <- fread(
  PB_META_FILE
)

META_LOCK <- fread(
  META_LOCK_FILE
)

GENE_MANIFEST <- fread(
  GENE_MANIFEST_FILE
)

required_pb_cols <- c(
  "pb_id",
  "donor",
  "condition",
  "cell_type",
  "n_cells",
  "library_size",
  "eligible_primary_DE"
)

missing_pb <- setdiff(
  required_pb_cols,
  names(PB_META)
)

if (length(missing_pb)) {
  stop(
    "Pseudobulk metadata lacks columns: ",
    paste(missing_pb, collapse = ", ")
  )
}

if (
  ncol(PB_COUNTS) !=
    nrow(PB_META)
) {
  stop(
    "Pseudobulk count matrix/sample metadata dimension mismatch."
  )
}

idx_pb <- match(
  colnames(PB_COUNTS),
  PB_META$pb_id
)

if (anyNA(idx_pb)) {
  stop(
    "Could not align pseudobulk count matrix to metadata."
  )
}

PB_META <- PB_META[
  idx_pb
]

if (!identical(
  colnames(PB_COUNTS),
  PB_META$pb_id
)) {
  stop(
    "Pseudobulk count matrix and metadata order mismatch."
  )
}

# ==============================================================================
# 6. DONOR-LEVEL COVARIATE LOCK
# ==============================================================================

if (!"sex" %in% names(META_LOCK)) {
  stop(
    "Frozen author metadata lacks sex column."
  )
}

DONOR_COV <- META_LOCK[
  ,
  .(
    condition_n =
      uniqueN(condition),
    condition =
      unique(condition)[1],
    sex_n =
      uniqueN(
        normalize_sex(sex)[
          !is.na(
            normalize_sex(sex)
          )
        ]
      ),
    sex =
      unique(
        normalize_sex(sex)[
          !is.na(
            normalize_sex(sex)
          )
        ]
      )[1]
  ),
  by = donor
]

if (any(
  DONOR_COV$condition_n != 1L
)) {
  stop(
    "At least one donor maps to multiple AF/CTRL labels."
  )
}

if (nrow(DONOR_COV) != EXPECTED_DONORS) {
  stop(
    "Expected 34 donors in frozen author metadata."
  )
}

if (
  sum(DONOR_COV$condition == "AF") !=
    EXPECTED_AF ||
  sum(DONOR_COV$condition == "CTRL") !=
    EXPECTED_CTRL
) {
  stop(
    "Donor AF/CTRL counts differ from 18/16."
  )
}

fwrite(
  DONOR_COV,
  file.path(
    QC_DIR,
    "STEP12B_donor_covariates_LOCKED.csv"
  )
)

# ==============================================================================
# 7. GENE ANNOTATION LOCK
# ==============================================================================

GENE_ANNO <- make_gene_annotation(
  GENE_MANIFEST,
  rownames(PB_COUNTS)
)

# ==============================================================================
# 8. PRIMARY ANALYSIS — CARDIOMYOCYTES
# ==============================================================================

CM <- run_pseudobulk_DE(
  target_celltype = PRIMARY_CELLTYPE,
  analysis_label = "PRIMARY",
  out_dir = PRIMARY_DIR,
  pb_counts = PB_COUNTS,
  pb_meta = PB_META,
  donor_covariates = DONOR_COV,
  gene_anno = GENE_ANNO
)

# ==============================================================================
# 9. PRE-SPECIFIED SECONDARY CONTEXT — MACROPHAGES
# ==============================================================================

MAC <- run_pseudobulk_DE(
  target_celltype = SECONDARY_CELLTYPE,
  analysis_label = "SECONDARY_CONTEXT",
  out_dir = SECONDARY_DIR,
  pb_counts = PB_COUNTS,
  pb_meta = PB_META,
  donor_covariates = DONOR_COV,
  gene_anno = GENE_ANNO
)

# ==============================================================================
# 10. PRIMARY-vs-SECONDARY SUMMARY
# ==============================================================================

extract_metric <- function(summary_dt, metric_name) {

  x <- summary_dt[
    metric == metric_name,
    value
  ]

  if (!length(x)) {
    return(NA_character_)
  }

  as.character(x[1])
}

SUMMARY_COMPARE <- data.table(
  cell_type = c(
    PRIMARY_CELLTYPE,
    SECONDARY_CELLTYPE
  ),
  role = c(
    "PRIMARY",
    "SECONDARY_CONTEXT"
  ),
  donors = c(
    extract_metric(
      CM$summary,
      "eligible_donors"
    ),
    extract_metric(
      MAC$summary,
      "eligible_donors"
    )
  ),
  design = c(
    CM$design_formula,
    MAC$design_formula
  ),
  genes_tested = c(
    extract_metric(
      CM$summary,
      "genes_tested"
    ),
    extract_metric(
      MAC$summary,
      "genes_tested"
    )
  ),
  FDR05_DEGs = c(
    extract_metric(
      CM$summary,
      "FDR05_DEGs"
    ),
    extract_metric(
      MAC$summary,
      "FDR05_DEGs"
    )
  ),
  FDR05_AF_up = c(
    extract_metric(
      CM$summary,
      "FDR05_AF_up"
    ),
    extract_metric(
      MAC$summary,
      "FDR05_AF_up"
    )
  ),
  FDR05_AF_down = c(
    extract_metric(
      CM$summary,
      "FDR05_AF_down"
    ),
    extract_metric(
      MAC$summary,
      "FDR05_AF_down"
    )
  )
)

fwrite(
  SUMMARY_COMPARE,
  file.path(
    QC_DIR,
    "STEP12B_primary_secondary_DE_summary.csv"
  )
)

# ==============================================================================
# 11. TECHNICAL READINESS
# ==============================================================================

cm_res <- CM$result
mac_res <- MAC$result

READINESS <- data.table(
  check = c(
    "STEP12A_fully_passed",
    "34_donors_locked",
    "18_AF_donors_locked",
    "16_CTRL_donors_locked",
    "Cardiomyocytes_34_donors",
    "Macrophages_34_donors",
    "CM_design_full_rank",
    "Macrophage_design_full_rank",
    "CM_genes_tested_gt_1000",
    "Macrophage_genes_tested_gt_1000",
    "CM_all_Pvalues_finite",
    "CM_all_FDR_finite",
    "Macrophage_all_Pvalues_finite",
    "Macrophage_all_FDR_finite",
    "CM_primary_result_written",
    "Macrophage_secondary_result_written",
    "candidate_gene_list_not_used_in_DE",
    "STEP12B_complete"
  ),
  pass = c(
    all(
      bool_pass(
        READY12A$pass
      )
    ),
    nrow(DONOR_COV) ==
      EXPECTED_DONORS,
    sum(
      DONOR_COV$condition ==
        "AF"
    ) ==
      EXPECTED_AF,
    sum(
      DONOR_COV$condition ==
        "CTRL"
    ) ==
      EXPECTED_CTRL,
    uniqueN(
      CM$sample_qc$donor
    ) ==
      EXPECTED_DONORS,
    uniqueN(
      MAC$sample_qc$donor
    ) ==
      EXPECTED_DONORS,
    qr(
      CM$design
    )$rank ==
      ncol(
        CM$design
      ),
    qr(
      MAC$design
    )$rank ==
      ncol(
        MAC$design
      ),
    nrow(
      cm_res
    ) >
      1000L,
    nrow(
      mac_res
    ) >
      1000L,
    all(
      is.finite(
        cm_res$PValue
      )
    ),
    all(
      is.finite(
        cm_res$FDR
      )
    ),
    all(
      is.finite(
        mac_res$PValue
      )
    ),
    all(
      is.finite(
        mac_res$FDR
      )
    ),
    file.exists(
      file.path(
        PRIMARY_DIR,
        "STEP12B_Cardiomyocytes_AF_vs_CTRL_all_genes.csv.gz"
      )
    ),
    file.exists(
      file.path(
        SECONDARY_DIR,
        "STEP12B_Macrophages_AF_vs_CTRL_all_genes.csv.gz"
      )
    ),
    TRUE,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP12B_readiness.csv"
  )
)

# ==============================================================================
# 12. METHOD PROVENANCE
# ==============================================================================

METHOD <- data.table(
  field = c(
    "dataset",
    "tissue",
    "statistical_unit",
    "primary_cell_type",
    "primary_cell_type_rationale",
    "secondary_context_cell_type",
    "DE_method",
    "normalization",
    "primary_model",
    "contrast",
    "genome_wide_DE_threshold",
    "large_effect_export_threshold",
    "candidate_gene_selection_in_STEP12B",
    "candidate_gene_integration_step",
    "reclustering",
    "reannotation",
    "cell_level_DE"
  ),
  value = c(
    "GSE255612 / SCP2489",
    "Human left atrium",
    "Biological donor/patient",
    "Cardiomyocytes",
    "Prespecified by independent upstream AF S-LDSC + SCAVENGE cellular localization",
    "Macrophages",
    "edgeR quasi-likelihood negative-binomial GLM",
    "TMM",
    "Prefer ~ sex + condition; fallback ~ condition only if sex non-estimable",
    "AF vs CTRL",
    "FDR < 0.05",
    "FDR < 0.05 and |log2FC| >= 0.5",
    "None; transcriptome-wide DE performed without upstream candidate filtering",
    "STEP12C after STEP12B is frozen",
    "No",
    "No",
    "No"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP12B_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP12B_sessionInfo.txt"
  )
)

# ==============================================================================
# 13. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP12B COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "PRIMARY — Cardiomyocytes\n"
)
print(
  CM$summary
)

cat(
  "\nSECONDARY CONTEXT — Macrophages\n"
)
print(
  MAC$summary
)

cat(
  "\nTechnical readiness:\n"
)
print(
  READINESS
)

if (all(
  READINESS$pass
)) {

  cat(
    "\nSTEP12B PASSED.\n",
    "DO NOT interpret candidate genes yet.\n",
    "Next: STEP12C will integrate the frozen CM DE table with the ",
    "pre-existing genetics/regulatory evidence matrix.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP12B has failed one or more readiness checks.\n",
    "Do not start STEP12C yet.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  PRIMARY_DIR,
  "\n",
  SECONDARY_DIR,
  "\n",
  sep = ""
)
