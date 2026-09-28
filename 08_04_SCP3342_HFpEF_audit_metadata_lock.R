# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP13A_V3_FIX_AUTHOR_CELLTYPE_HFPEF_H5AD_AUDIT.R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP13A V3 — HUMAN HFpEF snRNA-seq (SCP3342)
#           H5AD STRUCTURE AUDIT + METADATA LOCK + COUNT-LAYER AUDIT
#
# Project main line:
#   Shared genetics
#     -> regulatory loci / genes
#     -> GSE238242 S-LDSC + SCAVENGE: AF genetic risk localizes to CM
#     -> GSE255612 AF disease tissue: patient-level CM transcriptome
#     -> SCP3342 HFpEF myocardium: independent disease-tissue validation
#
# DATASET
#   Broad Single Cell Portal SCP3342
#   Human ventricular septal myocardium snRNA-seq
#   Published final processed dataset:
#      48,866 nuclei
#      36,601 genes
#      43 subjects = 19 HFpEF + 24 non-failing controls
#      14 author-defined cell types
#
# PURPOSE OF STEP13A
#   1) Audit the already-downloaded H5AD WITHOUT loading the full expression
#      matrix into RAM.
#   2) Read and freeze author obs/var metadata.
#   3) Identify donor, HFpEF/control status, cell type and sex fields.
#   4) Verify 48,866 nuclei / 36,601 genes / 19 HFpEF / 24 CTRL / 14 cell types.
#   5) Audit every possible expression-count location:
#         X
#         raw/X
#         layers/*
#      and determine whether values look like raw integer counts or normalized
#      expression.
#   6) Quantify donor x cell-type sample availability for STEP13B.
#
# IMPORTANT
#   - NO reclustering.
#   - NO reannotation.
#   - NO differential expression.
#   - NO candidate-gene filtering.
#   - NO full H5AD matrix import in this step.
#
# STEP13B will only start after the raw-count source is unambiguously identified.
#
# INPUT EXPECTED
#   D:/A/data/HFpEF_snRNAseq_single_cell_portal_10.14.2025.h5ad
#
# OUTPUT
#   D:/A/data/STEP13_HFpEF_SCP3342_snRNA/00_QC/
#   D:/A/data/STEP13_HFpEF_SCP3342_snRNA/01_METADATA/
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

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP13_HFpEF_SCP3342_snRNA"
)

QC_DIR <- file.path(OUT_ROOT, "00_QC")
META_DIR <- file.path(OUT_ROOT, "01_METADATA")
FIG_DIR <- file.path(OUT_ROOT, "02_FIGURES")

for (d in c(OUT_ROOT, QC_DIR, META_DIR, FIG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

EXPECTED_CELLS <- 48866L
EXPECTED_GENES <- 36601L
EXPECTED_DONORS <- 43L
EXPECTED_HFPEF <- 19L
EXPECTED_CTRL <- 24L
EXPECTED_CELLTYPES <- 14L

MIN_CELLS_PER_DONOR_CELLTYPE <- 20L

# Number of stored nonzero values sampled from a matrix candidate when checking
# whether it behaves like integer raw counts.
N_VALUE_SAMPLE <- 100000L

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

pkgs <- c(
  "data.table",
  "hdf5r",
  "ggplot2"
)

for (p in pkgs) {
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
library(hdf5r)
library(ggplot2)

# ==============================================================================
# 2. HELPERS
# ==============================================================================

norm_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "nan", "NaN")] <- NA_character_
  x
}

safe_attr <- function(obj, name, default = NA) {
  tryCatch(
    obj$attr_open(name)$read(),
    error = function(e) default
  )
}

safe_read_dataset <- function(obj) {
  tryCatch(
    obj[],
    error = function(e) NULL
  )
}

decode_categorical <- function(grp) {

  nm <- grp$ls(recursive = FALSE)$name

  if (!all(c("codes", "categories") %in% nm)) {
    return(NULL)
  }

  codes <- safe_read_dataset(grp[["codes"]])
  cats <- safe_read_dataset(grp[["categories"]])

  if (is.null(codes) || is.null(cats)) {
    return(NULL)
  }

  codes <- as.integer(codes)
  cats <- norm_text(cats)

  out <- rep(NA_character_, length(codes))

  ok <- is.finite(codes) & codes >= 0L & (codes + 1L) <= length(cats)

  out[ok] <- cats[codes[ok] + 1L]

  out
}

read_dataframe_group <- function(grp, expected_n = NULL) {

  listing <- grp$ls(recursive = FALSE)

  # AnnData dataframe index column is stored in the _index attribute.
  idx_name <- safe_attr(grp, "_index", default = "_index")
  idx_name <- as.character(idx_name)[1]

  result <- list()

  # Read index first.
  if (idx_name %in% listing$name) {
    idx_obj <- grp[[idx_name]]
    idx_val <- safe_read_dataset(idx_obj)
    if (!is.null(idx_val)) {
      result[[".index"]] <- norm_text(idx_val)
    }
  }

  for (nm in listing$name) {

    if (identical(nm, idx_name)) {
      next
    }

    obj <- grp[[nm]]

    val <- NULL

    if (inherits(obj, "H5D")) {

      tmp <- safe_read_dataset(obj)

      # Only retain one-dimensional metadata columns.
      if (!is.null(tmp) && length(dim(tmp)) <= 1L) {
        val <- tmp
      }

    } else if (inherits(obj, "H5Group")) {

      val <- decode_categorical(obj)
    }

    if (!is.null(val)) {

      val <- norm_text(val)

      if (is.null(expected_n) || length(val) == expected_n) {
        result[[nm]] <- val
      }
    }
  }

  if (!length(result)) {
    stop("Could not decode H5AD dataframe group.")
  }

  DT <- as.data.table(result)

  # Guarantee unique column names.
  setnames(DT, make.unique(names(DT)))

  DT
}

profile_columns <- function(DT) {

  rbindlist(
    lapply(
      names(DT),
      function(nm) {

        x <- norm_text(DT[[nm]])
        x2 <- x[!is.na(x)]

        vals <- tolower(unique(x2))

        data.table(
          column = nm,
          n_unique = uniqueN(x2),
          example = paste(head(unique(x2), 8), collapse = " | "),
          has_HFpEF = any(grepl("hfpef|preserved.?ejection", vals)),
          has_control = any(
            grepl(
              "non.?fail|control|healthy|normal|^nf$|donor",
              vals
            )
          ),
          has_CM = any(
            grepl(
              "cardiomyocyte|cardiac.?muscle|^cm$",
              vals
            )
          ),
          has_sex = any(vals %in% c("m", "f", "male", "female"))
        )
      }
    ),
    fill = TRUE
  )
}

pick_column <- function(PROF, patterns, structural) {

  scores <- vapply(
    PROF$column,
    function(nm) {
      z <- tolower(nm)
      max(
        c(
          0,
          vapply(
            seq_along(patterns),
            function(i) {
              if (grepl(patterns[i], z)) {
                length(patterns) - i + 1
              } else {
                0
              }
            },
            numeric(1)
          )
        )
      )
    },
    numeric(1)
  )

  idx <- which(structural & scores > 0)

  if (!length(idx)) {
    return(NA_character_)
  }

  PROF$column[idx[which.max(scores[idx])]]
}

normalize_condition <- function(x) {

  y <- tolower(norm_text(x))

  out <- rep(NA_character_, length(y))

  out[
    !is.na(y) &
      grepl(
        "hfpef|heart.?failure.?with.?preserved|preserved.?ejection",
        y
      )
  ] <- "HFpEF"

  out[
    !is.na(y) &
      grepl(
        "non.?fail|control|healthy|normal|^nf$|non.?hf",
        y
      )
  ] <- "CTRL"

  out
}

normalize_sex <- function(x) {

  y <- tolower(norm_text(x))

  out <- rep(NA_character_, length(y))

  out[y %in% c("f", "female", "woman")] <- "F"
  out[y %in% c("m", "male", "man")] <- "M"

  out
}

# ------------------------------------------------------------------------------
# Return dimensions/encoding for a possible matrix object.
# Supports dense HDF5 dataset and AnnData sparse CSR/CSC groups.
# ------------------------------------------------------------------------------

matrix_candidate_info <- function(h5, path) {

  exists <- tryCatch(
    h5$exists(path),
    error = function(e) FALSE
  )

  if (!exists) {
    return(NULL)
  }

  obj <- h5[[path]]

  encoding <- safe_attr(obj, "encoding-type", default = NA)
  encoding <- as.character(encoding)[1]

  shape <- safe_attr(obj, "shape", default = NULL)

  if (inherits(obj, "H5D")) {

    dims <- obj$dims

    # Sample stored values without loading entire dense matrix.
    n1 <- min(as.integer(dims[1]), 250L)
    n2 <- if (length(dims) >= 2L) min(as.integer(dims[2]), 250L) else 1L

    vals <- tryCatch(
      {
        if (length(dims) >= 2L) {
          as.numeric(obj[1:n1, 1:n2])
        } else {
          as.numeric(obj[1:min(length(obj[]), N_VALUE_SAMPLE)])
        }
      },
      error = function(e) numeric(0)
    )

    return(
      list(
        path = path,
        object_type = "dense_dataset",
        encoding = encoding,
        nrow = if (length(dims) >= 1L) as.integer(dims[1]) else NA_integer_,
        ncol = if (length(dims) >= 2L) as.integer(dims[2]) else NA_integer_,
        sample_values = vals
      )
    )
  }

  if (inherits(obj, "H5Group")) {

    members <- obj$ls(recursive = FALSE)$name

    if (is.null(shape) || length(shape) < 2L) {
      shape <- c(NA_integer_, NA_integer_)
    }

    vals <- numeric(0)

    if ("data" %in% members) {

      ds <- obj[["data"]]
      n <- min(as.integer(ds$dims[1]), N_VALUE_SAMPLE)

      if (n > 0L) {
        vals <- tryCatch(
          as.numeric(ds[1:n]),
          error = function(e) numeric(0)
        )
      }
    }

    return(
      list(
        path = path,
        object_type = "sparse_group",
        encoding = encoding,
        nrow = as.integer(shape[1]),
        ncol = as.integer(shape[2]),
        sample_values = vals
      )
    )
  }

  NULL
}

summarize_matrix_candidate <- function(x) {

  vals <- x$sample_values
  vals <- vals[is.finite(vals)]

  integer_fraction <- if (length(vals)) {
    mean(abs(vals - round(vals)) < 1e-8)
  } else {
    NA_real_
  }

  data.table(
    matrix_path = x$path,
    object_type = x$object_type,
    encoding_type = x$encoding,
    nrow = x$nrow,
    ncol = x$ncol,
    sampled_values = length(vals),
    sample_min = if (length(vals)) min(vals) else NA_real_,
    sample_median = if (length(vals)) median(vals) else NA_real_,
    sample_max = if (length(vals)) max(vals) else NA_real_,
    sampled_integer_fraction = integer_fraction,
    looks_like_raw_counts =
      length(vals) > 0L &&
      integer_fraction > 0.999 &&
      min(vals) >= 0
  )
}

# ==============================================================================
# 3. HARD FILE CHECK
# ==============================================================================

if (!file.exists(H5AD_FILE)) {
  stop(
    "HFpEF H5AD file not found:\n",
    H5AD_FILE
  )
}

FILE_QC <- data.table(
  file = basename(H5AD_FILE),
  path = H5AD_FILE,
  size_GB = file.info(H5AD_FILE)$size / 1024^3
)

fwrite(
  FILE_QC,
  file.path(QC_DIR, "STEP13A_input_file_QC.csv")
)

# ==============================================================================
# 4. OPEN H5AD — METADATA ONLY
# ==============================================================================

cat("\nOpening H5AD in read-only mode...\n")

h5 <- H5File$new(
  H5AD_FILE,
  mode = "r"
)

on.exit(
  try(h5$close_all(), silent = TRUE),
  add = TRUE
)

ROOT_LIST <- h5$ls(recursive = FALSE)

fwrite(
  as.data.table(ROOT_LIST),
  file.path(QC_DIR, "STEP13A_H5AD_root_structure.csv")
)

if (!h5$exists("obs")) {
  stop("H5AD lacks /obs group.")
}

if (!h5$exists("var")) {
  stop("H5AD lacks /var group.")
}

# ==============================================================================
# 5. READ OBS / VAR WITHOUT READING X
# ==============================================================================

cat("Reading obs metadata only...\n")

OBS <- read_dataframe_group(
  h5[["obs"]],
  expected_n = EXPECTED_CELLS
)

cat("Reading var metadata only...\n")

VAR <- read_dataframe_group(
  h5[["var"]],
  expected_n = EXPECTED_GENES
)

# If exact expected_n caused metadata fields to be omitted because file dimensions
# differ, fail explicitly rather than silently changing study identity.
if (nrow(OBS) != EXPECTED_CELLS) {
  stop(
    "Expected 48,866 obs rows; decoded ",
    nrow(OBS),
    "."
  )
}

if (nrow(VAR) != EXPECTED_GENES) {
  stop(
    "Expected 36,601 var rows; decoded ",
    nrow(VAR),
    "."
  )
}

fwrite(
  OBS,
  file.path(
    META_DIR,
    "SCP3342_obs_author_metadata_RAW.csv.gz"
  ),
  compress = "gzip"
)

fwrite(
  VAR,
  file.path(
    META_DIR,
    "SCP3342_var_gene_metadata_RAW.csv.gz"
  ),
  compress = "gzip"
)

# ==============================================================================
# 6. IDENTIFY METADATA COLUMNS
# ==============================================================================

PROF <- profile_columns(OBS)

fwrite(
  PROF,
  file.path(
    QC_DIR,
    "STEP13A_obs_column_profile.csv"
  )
)

# Cell type.
# IMPORTANT:
# SCP3342 contains ontology-normalized "cell_type" fields with 13 unique
# ontology classes, while the publication uses 14 AUTHOR-DEFINED cell types.
# We must preserve the author-defined 14-class annotation for all downstream
# patient-level pseudobulk analyses.
#
# Preferred order:
#   1) original_cell_type          (human-readable author label)
#   2) original_cell_type_label    (numbered author label)
#   3) another column with exactly 14 unique values AND a CM label
# Only if these are unavailable do we fall back to generic cell-type columns.

CELLTYPE_COL <- NA_character_

preferred_author_ct <- c(
  "original_cell_type",
  "original_cell_type_label"
)

for (nm in preferred_author_ct) {

  if (nm %in% PROF$column) {

    rr <- PROF[
      column == nm
    ]

    if (
      nrow(rr) == 1L &&
      rr$n_unique == EXPECTED_CELLTYPES &&
      isTRUE(rr$has_CM)
    ) {

      CELLTYPE_COL <- nm
      break
    }
  }
}

if (is.na(CELLTYPE_COL)) {

  idx_author <- which(
    PROF$n_unique == EXPECTED_CELLTYPES &
      PROF$has_CM
  )

  if (length(idx_author) == 1L) {
    CELLTYPE_COL <- PROF$column[idx_author]
  }
}

if (is.na(CELLTYPE_COL)) {

  CELLTYPE_COL <- pick_column(
    PROF,
    patterns = c(
      "^original_cell_type$",
      "^original_cell_type_label$",
      "^cell_type$",
      "^celltype$",
      "cell_type",
      "celltype",
      "annotation",
      "cell.*type",
      "cluster"
    ),
    structural =
      PROF$n_unique >= 8L &
      PROF$n_unique <= 30L
  )
}

# Hard guard: do not silently accept an ontology-collapsed annotation when a
# 14-class author annotation is expected.
if (!is.na(CELLTYPE_COL)) {

  n_selected_ct <- PROF[
    column == CELLTYPE_COL,
    n_unique
  ][1]

  if (
    is.na(n_selected_ct) ||
    n_selected_ct != EXPECTED_CELLTYPES
  ) {

    candidates14 <- PROF[
      n_unique == EXPECTED_CELLTYPES &
        has_CM == TRUE
    ]

    fwrite(
      candidates14,
      file.path(
        QC_DIR,
        "STEP13A_author_celltype_14class_candidates.csv"
      )
    )

    stop(
      paste0(
        "Selected cell-type column '",
        CELLTYPE_COL,
        "' has ",
        n_selected_ct,
        " unique classes, but SCP3342 requires 14 author-defined classes.\n",
        "The script will not use ontology-collapsed cell types.\n",
        "Inspect STEP13A_author_celltype_14class_candidates.csv."
      )
    )
  }
}

# Condition.
CONDITION_COL <- pick_column(
  PROF,
  patterns = c(
    "^condition$",
    "^diagnosis$",
    "disease",
    "condition",
    "diagnosis",
    "phenotype",
    "status",
    "group"
  ),
  structural =
    PROF$n_unique >= 2L &
    PROF$n_unique <= 10L &
    PROF$has_HFpEF &
    PROF$has_control
)

if (is.na(CONDITION_COL)) {

  idx <- which(
    PROF$n_unique >= 2L &
    PROF$n_unique <= 10L &
    PROF$has_HFpEF &
    PROF$has_control
  )

  if (length(idx) == 1L) {
    CONDITION_COL <- PROF$column[idx]
  }
}

# Donor.
DONOR_COL <- pick_column(
  PROF,
  patterns = c(
    "^donor$",
    "^sample$",
    "subject",
    "patient",
    "individual",
    "donor",
    "sample_id",
    "sample"
  ),
  structural =
    PROF$n_unique >= 35L &
    PROF$n_unique <= 50L
)

if (is.na(DONOR_COL)) {

  idx <- which(
    PROF$n_unique >= 35L &
    PROF$n_unique <= 50L
  )

  if (length(idx)) {

    d <- abs(
      PROF$n_unique[idx] -
        EXPECTED_DONORS
    )

    best <- idx[d == min(d)]

    if (length(best) == 1L) {
      DONOR_COL <- PROF$column[best]
    }
  }
}

# Sex is useful later as a covariate but is not required to identify the dataset.
SEX_COL <- pick_column(
  PROF,
  patterns = c(
    "^sex$",
    "biological_sex",
    "gender",
    "sex"
  ),
  structural =
    PROF$n_unique >= 2L &
    PROF$n_unique <= 4L &
    PROF$has_sex
)

DETECTION <- data.table(
  role = c(
    "donor",
    "condition",
    "cell_type",
    "sex"
  ),
  column = c(
    DONOR_COL,
    CONDITION_COL,
    CELLTYPE_COL,
    SEX_COL
  )
)

fwrite(
  DETECTION,
  file.path(
    QC_DIR,
    "STEP13A_detected_metadata_columns.csv"
  )
)

if (any(is.na(DETECTION$column[1:3]))) {

  print(DETECTION)

  stop(
    paste0(
      "Could not unambiguously identify donor/condition/cell type.\n",
      "Upload STEP13A_obs_column_profile.csv and ",
      "STEP13A_detected_metadata_columns.csv."
    )
  )
}

# ==============================================================================
# 7. LOCK AUTHOR METADATA
# ==============================================================================

META <- copy(OBS)

META[
  ,
  donor :=
    norm_text(
      get(DONOR_COL)
    )
]

META[
  ,
  condition_original :=
    norm_text(
      get(CONDITION_COL)
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
      get(CELLTYPE_COL)
    )
]

if (!is.na(SEX_COL)) {

  META[
    ,
    sex_original :=
      norm_text(
        get(SEX_COL)
      )
  ]

  META[
    ,
    sex :=
      normalize_sex(
        sex_original
      )
  ]

} else {

  META[
    ,
    sex_original :=
      NA_character_
  ]

  META[
    ,
    sex :=
      NA_character_
  ]
}

# Stable cell ID = AnnData obs index.
if (".index" %in% names(META)) {

  META[
    ,
    cell_id :=
      norm_text(
        .index
      )
  ]

} else {

  META[
    ,
    cell_id :=
      paste0(
        "cell_",
        seq_len(.N)
      )
  ]
}

if (anyNA(META$condition)) {

  bad <- unique(
    META[
      is.na(condition),
      condition_original
    ]
  )

  stop(
    paste0(
      "Unknown HFpEF/control labels remain after normalization:\n",
      paste(bad, collapse = " | ")
    )
  )
}

if (anyNA(META$donor) || anyNA(META$cell_type)) {
  stop("Missing donor or cell-type labels in locked metadata.")
}

if (anyDuplicated(META$cell_id)) {
  stop("Duplicated AnnData cell IDs.")
}

DONOR_CONDITION <- unique(
  META[
    ,
    .(
      donor,
      condition
    )
  ]
)

DONOR_CHECK <- DONOR_CONDITION[
  ,
  .N,
  by = donor
]

if (any(DONOR_CHECK$N != 1L)) {
  stop("At least one donor maps to multiple HFpEF/control labels.")
}

N_DONORS <- uniqueN(META$donor)
N_HFPEF <- uniqueN(META[condition == "HFpEF", donor])
N_CTRL <- uniqueN(META[condition == "CTRL", donor])
N_CT <- uniqueN(META$cell_type)

fwrite(
  META,
  file.path(
    META_DIR,
    "SCP3342_author_metadata_LOCKED.csv.gz"
  ),
  compress = "gzip"
)

# ==============================================================================
# 8. CELL-TYPE / DONOR READINESS
# ==============================================================================

CT_SUM <- META[
  ,
  .(
    n_cells = .N,
    n_donors = uniqueN(donor),
    n_HFpEF_cells = sum(condition == "HFpEF"),
    n_CTRL_cells = sum(condition == "CTRL"),
    n_HFpEF_donors = uniqueN(donor[condition == "HFpEF"]),
    n_CTRL_donors = uniqueN(donor[condition == "CTRL"])
  ),
  by = cell_type
]

setorder(
  CT_SUM,
  -n_cells
)

fwrite(
  CT_SUM,
  file.path(
    META_DIR,
    "SCP3342_author_celltype_summary.csv"
  )
)

DONOR_CT <- META[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    donor,
    condition,
    cell_type
  )
]

DONOR_CT[
  ,
  eligible_primary_DE :=
    n_cells >=
      MIN_CELLS_PER_DONOR_CELLTYPE
]

CT_READY <- DONOR_CT[
  ,
  .(
    n_donors_all =
      as.integer(
        uniqueN(
          donor
        )
      ),
    n_HFpEF_all =
      as.integer(
        uniqueN(
          donor[
            condition ==
              "HFpEF"
          ]
        )
      ),
    n_CTRL_all =
      as.integer(
        uniqueN(
          donor[
            condition ==
              "CTRL"
          ]
        )
      ),
    n_eligible_HFpEF =
      as.integer(
        uniqueN(
          donor[
            condition ==
              "HFpEF" &
              eligible_primary_DE
          ]
        )
      ),
    n_eligible_CTRL =
      as.integer(
        uniqueN(
          donor[
            condition ==
              "CTRL" &
              eligible_primary_DE
          ]
        )
      ),
    median_cells_per_donor =
      as.numeric(
        median(
          as.numeric(
            n_cells
          )
        )
      ),
    min_cells_per_donor =
      as.numeric(
        min(
          as.numeric(
            n_cells
          )
        )
      )
  ),
  by = cell_type
]

CT_READY[
  ,
  ready_for_STEP13B :=
    n_eligible_HFpEF >=
      5L &
    n_eligible_CTRL >=
      5L
]

setorder(
  CT_READY,
  -ready_for_STEP13B,
  -n_eligible_HFpEF,
  -n_eligible_CTRL
)

fwrite(
  CT_READY,
  file.path(
    QC_DIR,
    "STEP13A_celltype_patient_level_readiness.csv"
  )
)

# Cardiomyocyte label.
CM_LABELS <- CT_READY[
  grepl(
    "cardiomyocyte|cardiac.?muscle|^cm$",
    cell_type,
    ignore.case = TRUE
  ),
  cell_type
]

fwrite(
  data.table(
    target = "Cardiomyocyte",
    detected_author_label =
      if (length(CM_LABELS)) {
        paste(CM_LABELS, collapse = " | ")
      } else {
        NA_character_
      }
  ),
  file.path(
    QC_DIR,
    "STEP13A_primary_celltype_label.csv"
  )
)

# ==============================================================================
# 9. AUDIT X / RAW/X / LAYERS — DO NOT LOAD FULL MATRIX
# ==============================================================================

candidate_paths <- c("X")

if (h5$exists("raw") && h5[["raw"]]$exists("X")) {
  candidate_paths <- c(candidate_paths, "raw/X")
}

if (h5$exists("layers")) {

  layer_names <- h5[["layers"]]$ls(recursive = FALSE)$name

  if (length(layer_names)) {
    candidate_paths <- c(
      candidate_paths,
      paste0("layers/", layer_names)
    )
  }
}

MATRIX_INFO <- rbindlist(
  lapply(
    candidate_paths,
    function(path) {

      info <- matrix_candidate_info(
        h5,
        path
      )

      if (is.null(info)) {
        return(NULL)
      }

      summarize_matrix_candidate(
        info
      )
    }
  ),
  fill = TRUE
)

if (!nrow(MATRIX_INFO)) {
  stop("No expression matrix candidate could be audited.")
}

# AnnData matrix orientation is cells x genes.
MATRIX_INFO[
  ,
  dimensions_match_obs_var :=
    nrow ==
      EXPECTED_CELLS &
    ncol ==
      EXPECTED_GENES
]

# Raw can have a different var count; preserve that fact rather than rejecting it.
MATRIX_INFO[
  ,
  plausible_count_source :=
    looks_like_raw_counts &
    nrow ==
      EXPECTED_CELLS
]

fwrite(
  MATRIX_INFO,
  file.path(
    QC_DIR,
    "STEP13A_expression_matrix_candidate_audit.csv"
  )
)

COUNT_CANDIDATES <- MATRIX_INFO[
  plausible_count_source ==
    TRUE
]

# We only lock automatically when exactly one candidate looks like raw counts.
COUNT_SOURCE_LOCKED <- nrow(COUNT_CANDIDATES) == 1L

COUNT_SOURCE <- if (COUNT_SOURCE_LOCKED) {
  COUNT_CANDIDATES$matrix_path[1]
} else {
  NA_character_
}

COUNT_SOURCE_TABLE <- data.table(
  locked = COUNT_SOURCE_LOCKED,
  selected_count_matrix_path = COUNT_SOURCE,
  n_plausible_count_candidates = nrow(COUNT_CANDIDATES)
)

fwrite(
  COUNT_SOURCE_TABLE,
  file.path(
    QC_DIR,
    "STEP13A_count_source_lock.csv"
  )
)

# ==============================================================================
# 10. DATASET STRUCTURE QC
# ==============================================================================

STRUCTURE_QC <- data.table(
  metric = c(
    "n_cells",
    "n_genes",
    "n_donors",
    "n_HFpEF_donors",
    "n_CTRL_donors",
    "n_author_celltypes",
    "cardiomyocyte_labels_detected",
    "count_source_uniquely_locked"
  ),
  value = c(
    nrow(META),
    nrow(VAR),
    N_DONORS,
    N_HFPEF,
    N_CTRL,
    N_CT,
    length(CM_LABELS),
    COUNT_SOURCE_LOCKED
  )
)

fwrite(
  STRUCTURE_QC,
  file.path(
    QC_DIR,
    "STEP13A_dataset_structure_QC.csv"
  )
)

# ==============================================================================
# 11. FIGURE — AUTHOR CELL-TYPE COMPOSITION
# ==============================================================================

PLOT <- copy(
  CT_SUM
)

PLOT[
  ,
  cell_type_plot :=
    factor(
      cell_type,
      levels = rev(cell_type)
    )
]

p <- ggplot(
  PLOT,
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
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.55)
  )

ggsave(
  file.path(
    FIG_DIR,
    "Figure_STEP13A_SCP3342_author_celltype_counts.tiff"
  ),
  p,
  width = 6.2,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================================================================
# 12. READINESS
# ==============================================================================

CM_UNIQUE <- length(CM_LABELS) == 1L

CM_READY <- FALSE

if (CM_UNIQUE) {

  tmp <- CT_READY[
    cell_type ==
      CM_LABELS[1],
    ready_for_STEP13B
  ]

  if (length(tmp) == 1L) {
    CM_READY <- isTRUE(tmp)
  }
}

READINESS <- data.table(
  check = c(
    "H5AD_file_present",
    "48866_nuclei_present",
    "36601_genes_present",
    "43_donors_present",
    "19_HFpEF_donors_present",
    "24_CTRL_donors_present",
    "14_author_celltypes_present",
    "donor_column_locked",
    "condition_column_locked",
    "celltype_column_locked",
    "condition_normalization_complete",
    "cardiomyocyte_label_unique",
    "cardiomyocyte_ready_for_patient_level_DE",
    "expression_matrix_candidates_audited",
    "raw_count_source_uniquely_locked",
    "STEP13A_complete"
  ),
  pass = c(
    file.exists(H5AD_FILE),
    nrow(META) == EXPECTED_CELLS,
    nrow(VAR) == EXPECTED_GENES,
    N_DONORS == EXPECTED_DONORS,
    N_HFPEF == EXPECTED_HFPEF,
    N_CTRL == EXPECTED_CTRL,
    N_CT == EXPECTED_CELLTYPES,
    !is.na(DONOR_COL),
    !is.na(CONDITION_COL),
    !is.na(CELLTYPE_COL),
    all(!is.na(META$condition)),
    CM_UNIQUE,
    CM_READY,
    nrow(MATRIX_INFO) >= 1L,
    COUNT_SOURCE_LOCKED,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP13A_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "dataset",
    "portal",
    "tissue",
    "assay",
    "published_final_nuclei",
    "published_final_HFpEF",
    "published_final_controls",
    "published_celltypes",
    "metadata_policy",
    "celltype_policy",
    "matrix_policy",
    "statistical_unit_next_step",
    "primary_celltype_next_step",
    "candidate_gene_filtering_STEP13A",
    "differential_expression_STEP13A"
  ),
  value = c(
    "SCP3342",
    "Broad Single Cell Portal",
    "Human ventricular septal myocardium",
    "snRNA-seq",
    as.character(EXPECTED_CELLS),
    as.character(EXPECTED_HFPEF),
    as.character(EXPECTED_CTRL),
    as.character(EXPECTED_CELLTYPES),
    "Preserve author metadata; no reclustering or reannotation",
    "Use the 14-class author-defined original_cell_type annotation; do not use the 13-class ontology-collapsed cell_type field",
    "Audit H5AD X/raw/X/layers without full matrix import; lock raw-count source before pseudobulk",
    "Biological donor/patient",
    "Cardiomyocyte, prespecified by upstream AF cellular localization and cross-disease study design",
    "None",
    "Not performed"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP13A_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP13A_sessionInfo.txt"
  )
)

# ==============================================================================
# 13. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP13A HFpEF SCP3342 AUDIT COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat("Detected metadata columns:\n")
print(DETECTION)

cat("\nDataset structure:\n")
print(STRUCTURE_QC)

cat("\nExpression matrix candidates:\n")
print(MATRIX_INFO)

cat("\nCount-source lock:\n")
print(COUNT_SOURCE_TABLE)

cat("\nCell-type patient-level readiness:\n")
print(CT_READY)

cat("\nTechnical readiness:\n")
print(READINESS)

if (all(READINESS$pass)) {

  cat(
    "\nSTEP13A PASSED.\n",
    "Next: STEP13B will build donor×cell-type pseudobulk directly from the ",
    "locked raw-count H5AD matrix and run HFpEF vs CTRL patient-level DE.\n",
    sep = ""
  )

} else {

  cat(
    "\nSTEP13A has one or more failed checks.\n",
    "Do NOT start STEP13B until the metadata/count-source audit is resolved.\n",
    sep = ""
  )
}

cat(
  "\nUPLOAD AFTER COMPLETION:\n",
  QC_DIR,
  "\n",
  sep = ""
)
