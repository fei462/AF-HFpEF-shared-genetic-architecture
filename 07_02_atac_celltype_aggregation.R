# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11B_GSE238242_ATAC_CELLTYPE_AGGREGATION_V2.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11B V2 — GSE238242 paired-barcode lock + cell-type ATAC aggregation
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# WHY THIS STEP
#   STEP11A passed. The 7 RNA and 7 ATAC matrices are paired by donor,
#   and the author metadata contains complete donor/Rhythm/cell_type labels.
#
#   For the next S-LDSC stage, rebuilding/reclustering a full Seurat/Signac
#   object is unnecessary and memory-inefficient. Instead, STEP11B:
#     1) locks the author-provided labels;
#     2) proves RNA/ATAC/metadata barcodes pair exactly;
#     3) proves the RNA gene and ATAC peak universes are consistent;
#     4) streams the large ATAC TSV files in small batches;
#     5) creates peak x cell-type pseudobulk detection/count summaries
#        for ALL donors (primary) and AF/SR strata (sensitivity);
#     6) preserves the cell manifest for later AF-priority SCAVENGE.
#
# IMPORTANT
#   - NO re-clustering
#   - NO re-annotation
#   - NO biological filtering
#   - NO peak threshold is chosen for S-LDSC yet
#   - NO liftOver yet
#   - NO SCAVENGE yet
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"

STEP11_ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

QC_DIR <- file.path(
  STEP11_ROOT,
  "00_QC"
)

EXTRACT_DIR <- file.path(
  STEP11_ROOT,
  "01_EXTRACTED"
)

OBJECT_DIR <- file.path(
  STEP11_ROOT,
  "02_OBJECTS"
)

SLDSC_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC"
)

AGG_DIR <- file.path(
  SLDSC_DIR,
  "00_ATAC_AGGREGATES"
)

for (d in c(
  QC_DIR,
  OBJECT_DIR,
  SLDSC_DIR,
  AGG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

META_GZ <- file.path(
  DATA_ROOT,
  "GSE238242_snAF.metadata.tsv.gz"
)

FILE_MANIFEST_CSV <- file.path(
  QC_DIR,
  "STEP11A_extracted_file_manifest.csv"
)

# Streaming batch size.
# 500 is conservative for ~2,400 cells/sample.
BATCH_N_PEAKS <- 500L

# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)

# ============================================================
# 2. HARD INPUT CHECK
# ============================================================

required <- c(
  metadata = META_GZ,
  manifest = FILE_MANIFEST_CSV
)

missing <- required[
  !file.exists(required)
]

if (length(missing)) {
  stop(
    paste0(
      "Missing STEP11B input(s):\n",
      paste(
        names(missing),
        missing,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

# ============================================================
# 3. LOAD AUTHOR METADATA + STEP11A MANIFEST
# ============================================================

META <- fread(
  META_GZ,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  showProgress = FALSE
)

MANIFEST <- fread(
  FILE_MANIFEST_CSV,
  showProgress = FALSE
)

required_meta_cols <- c(
  "V1",
  "cell_type",
  "sample",
  "sex",
  "Rhythm"
)

if (!all(required_meta_cols %in% names(META))) {
  stop(
    paste0(
      "Expected author metadata columns are missing.\n",
      "Required: ",
      paste(required_meta_cols, collapse = ", "),
      "\nObserved: ",
      paste(names(META), collapse = ", ")
    )
  )
}

if (anyDuplicated(META$V1)) {
  stop(
    "Author metadata barcode column V1 contains duplicated barcodes."
  )
}

if (nrow(META) != 11986L) {
  warning(
    "Expected 11,986 author-annotated cells, observed ",
    nrow(META),
    ". The script will continue but this must be reviewed."
  )
}

# Freeze clearer names without altering source file.
META_LOCKED <- copy(
  META
)

setnames(
  META_LOCKED,
  c(
    "V1",
    "sample",
    "Rhythm"
  ),
  c(
    "barcode",
    "donor",
    "rhythm"
  )
)

META_LOCKED[
  ,
  cell_type :=
    as.character(
      cell_type
    )
]

META_LOCKED[
  ,
  donor :=
    as.character(
      donor
    )
]

META_LOCKED[
  ,
  rhythm :=
    as.character(
      rhythm
    )
]

META_LOCKED[
  ,
  barcode :=
    as.character(
      barcode
    )
]

saveRDS(
  META_LOCKED,
  file.path(
    OBJECT_DIR,
    "GSE238242_author_metadata_LOCKED.rds"
  ),
  compress = TRUE
)

# ============================================================
# 4. EXPECTED DESIGN
# ============================================================

EXPECTED_DONORS <- c(
  "CF69",
  "CF77",
  "CF89",
  "CF91",
  "CF93",
  "CF97",
  "CF102"
)

EXPECTED_RHYTHM <- data.table(
  donor = EXPECTED_DONORS,
  expected_rhythm = c(
    "SR",
    "SR",
    "SR",
    "SR",
    "AF",
    "AF",
    "AF"
  )
)

observed_design <- unique(
  META_LOCKED[
    ,
    .(
      donor,
      rhythm,
      sex
    )
  ]
)

DESIGN_QC <- merge(
  EXPECTED_RHYTHM,
  observed_design,
  by = "donor",
  all.x = TRUE
)

DESIGN_QC[
  ,
  rhythm_match :=
    expected_rhythm == rhythm
]

if (
  nrow(DESIGN_QC) != 7L ||
  !all(DESIGN_QC$rhythm_match)
) {
  stop(
    "Donor/rhythm design does not match the frozen STEP11A design."
  )
}

fwrite(
  DESIGN_QC,
  file.path(
    QC_DIR,
    "STEP11B_donor_design_QC.csv"
  )
)

# ============================================================
# 5. HELPERS
# ============================================================

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
      "Empty matrix file: ",
      path
    )
  }

  out <- strsplit(
    line,
    "\t",
    fixed = TRUE
  )[[1]]

  clean_quotes(
    out
  )
}

read_feature_names_gz <- function(path) {

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  # skip header
  invisible(
    readLines(
      con,
      n = 1L,
      warn = FALSE
    )
  )

  ans <- character()

  repeat {

    lines <- readLines(
      con,
      n = 5000L,
      warn = FALSE
    )

    if (!length(lines)) {
      break
    }

    first_fields <- sub(
      "\t.*$",
      "",
      lines
    )

    first_fields <- clean_quotes(
      first_fields
    )

    ans <- c(
      ans,
      first_fields
    )
  }

  ans
}

get_one_file <- function(
  donor,
  modality
) {

  donor_value <- donor
  modality_value <- modality

  hit <- MANIFEST[
    donor == donor_value &
      modality == modality_value,
    full_path
  ]

  if (length(hit) != 1L) {
    stop(
      "Expected exactly one ",
      modality,
      " file for ",
      donor,
      "; found ",
      length(hit)
    )
  }

  hit
}

# ============================================================
# 6. DONOR-LEVEL RNA/ATAC/METADATA BARCODE PAIRING
# ============================================================

PAIR_QC_LIST <- list()
CELL_MANIFEST_LIST <- list()

for (d in EXPECTED_DONORS) {

  rna_file <- get_one_file(
    d,
    "RNA"
  )

  atac_file <- get_one_file(
    d,
    "ATAC"
  )

  rna_bc <- read_header_barcodes(
    rna_file
  )

  atac_bc <- read_header_barcodes(
    atac_file
  )

  meta_d <- META_LOCKED[
    donor == d
  ]

  meta_bc <- meta_d$barcode

  rna_atac_exact_order <- identical(
    rna_bc,
    atac_bc
  )

  rna_meta_set_equal <- setequal(
    rna_bc,
    meta_bc
  )

  atac_meta_set_equal <- setequal(
    atac_bc,
    meta_bc
  )

  all_unique <- (
    !anyDuplicated(rna_bc) &&
    !anyDuplicated(atac_bc) &&
    !anyDuplicated(meta_bc)
  )

  PAIR_QC_LIST[[
    length(PAIR_QC_LIST) + 1L
  ]] <- data.table(
    donor = d,
    rhythm = unique(meta_d$rhythm),
    n_metadata = length(meta_bc),
    n_RNA_barcodes = length(rna_bc),
    n_ATAC_barcodes = length(atac_bc),
    RNA_ATAC_exact_same_order = rna_atac_exact_order,
    RNA_metadata_set_equal = rna_meta_set_equal,
    ATAC_metadata_set_equal = atac_meta_set_equal,
    all_barcodes_unique = all_unique
  )

  if (
    !rna_atac_exact_order ||
    !rna_meta_set_equal ||
    !atac_meta_set_equal ||
    !all_unique
  ) {
    stop(
      "Barcode pairing failed for donor ",
      d,
      ". Inspect STEP11B_barcode_pairing_QC.csv after rerun."
    )
  }

  # Order metadata exactly as the matrix columns.
  idx <- match(
    atac_bc,
    meta_d$barcode
  )

  meta_ordered <- meta_d[
    idx
  ]

  meta_ordered[
    ,
    matrix_column_index :=
      seq_len(
        .N
      )
  ]

  CELL_MANIFEST_LIST[[
    length(CELL_MANIFEST_LIST) + 1L
  ]] <- meta_ordered
}

PAIR_QC <- rbindlist(
  PAIR_QC_LIST,
  fill = TRUE
)

CELL_MANIFEST <- rbindlist(
  CELL_MANIFEST_LIST,
  fill = TRUE
)

fwrite(
  PAIR_QC,
  file.path(
    QC_DIR,
    "STEP11B_barcode_pairing_QC.csv"
  )
)

fwrite(
  CELL_MANIFEST,
  file.path(
    OBJECT_DIR,
    "STEP11B_cell_manifest.csv"
  )
)

if (!all(
  PAIR_QC$RNA_ATAC_exact_same_order &
  PAIR_QC$RNA_metadata_set_equal &
  PAIR_QC$ATAC_metadata_set_equal &
  PAIR_QC$all_barcodes_unique
)) {
  stop(
    "One or more barcode-pairing checks failed."
  )
}

cat(
  "\nBarcode pairing: PASS for all 7 donors.\n"
)

# ============================================================
# 7. RNA GENE UNIVERSE CONSISTENCY
# ============================================================

RNA_GENE_QC_LIST <- list()
RNA_GENE_REFERENCE <- NULL

for (d in EXPECTED_DONORS) {

  rna_file <- get_one_file(
    d,
    "RNA"
  )

  genes <- read_feature_names_gz(
    rna_file
  )

  if (is.null(RNA_GENE_REFERENCE)) {
    RNA_GENE_REFERENCE <- genes
  }

  identical_to_reference <- identical(
    genes,
    RNA_GENE_REFERENCE
  )

  RNA_GENE_QC_LIST[[
    length(RNA_GENE_QC_LIST) + 1L
  ]] <- data.table(
    donor = d,
    n_genes = length(genes),
    identical_gene_order_to_reference = identical_to_reference
  )

  if (!identical_to_reference) {
    stop(
      "RNA gene universe/order differs for donor ",
      d
    )
  }
}

RNA_GENE_QC <- rbindlist(
  RNA_GENE_QC_LIST
)

fwrite(
  RNA_GENE_QC,
  file.path(
    QC_DIR,
    "STEP11B_RNA_gene_universe_QC.csv"
  )
)

saveRDS(
  RNA_GENE_REFERENCE,
  file.path(
    OBJECT_DIR,
    "GSE238242_RNA_gene_universe.rds"
  ),
  compress = TRUE
)

cat(
  "RNA gene universe: PASS (",
  length(RNA_GENE_REFERENCE),
  " genes).\n",
  sep = ""
)

# ============================================================
# 8. ATAC PEAK UNIVERSE FROM FIRST DONOR
# ============================================================

FIRST_ATAC <- get_one_file(
  EXPECTED_DONORS[1],
  "ATAC"
)

PEAKS <- read_feature_names_gz(
  FIRST_ATAC
)

N_PEAKS <- length(
  PEAKS
)

if (N_PEAKS < 100000L) {
  stop(
    "Unexpectedly small ATAC peak universe: ",
    N_PEAKS
  )
}

if (anyDuplicated(PEAKS)) {
  stop(
    "ATAC peak universe contains duplicated peak IDs."
  )
}

saveRDS(
  PEAKS,
  file.path(
    OBJECT_DIR,
    "GSE238242_ATAC_peak_universe_GRCh38.rds"
  ),
  compress = TRUE
)

cat(
  "ATAC peak universe loaded: ",
  format(N_PEAKS, big.mark = ","),
  " peaks.\n",
  sep = ""
)

# ============================================================
# 9. PARSE PEAK COORDINATES
# ============================================================

peak_parts <- tstrsplit(
  PEAKS,
  "-",
  fixed = TRUE
)

if (length(peak_parts) != 3L) {
  stop(
    "ATAC peak IDs are not consistently formatted as chr-start-end."
  )
}

PEAK_COORDS <- data.table(
  peak_index = seq_along(PEAKS),
  peak_id = PEAKS,
  chr = peak_parts[[1]],
  start = suppressWarnings(
    as.integer(
      peak_parts[[2]]
    )
  ),
  end = suppressWarnings(
    as.integer(
      peak_parts[[3]]
    )
  )
)

PEAK_COORDS[
  ,
  valid_coordinate :=
    !is.na(start) &
    !is.na(end) &
    end >= start &
    grepl(
      "^chr([0-9]+|X|Y)$",
      chr
    )
]

PEAK_COORDS[
  ,
  autosomal :=
    grepl(
      "^chr([1-9]|1[0-9]|2[0-2])$",
      chr
    )
]

fwrite(
  PEAK_COORDS,
  file.path(
    OBJECT_DIR,
    "GSE238242_ATAC_peak_coordinates_GRCh38.csv"
  )
)

PEAK_COORD_QC <- PEAK_COORDS[
  ,
  .(
    n_peaks = .N,
    n_valid = sum(valid_coordinate),
    n_autosomal = sum(autosomal),
    n_non_autosomal = sum(!autosomal),
    n_invalid = sum(!valid_coordinate)
  )
]

fwrite(
  PEAK_COORD_QC,
  file.path(
    QC_DIR,
    "STEP11B_peak_coordinate_QC.csv"
  )
)

if (PEAK_COORD_QC$n_invalid > 0L) {
  stop(
    "Invalid ATAC peak coordinates were detected."
  )
}

# ============================================================
# 10. CELL TYPES + CELL COUNTS
# ============================================================

CELL_TYPES <- sort(
  unique(
    META_LOCKED$cell_type
  )
)

if (length(CELL_TYPES) != 12L) {
  warning(
    "Expected 12 author cell types; observed ",
    length(CELL_TYPES)
  )
}

CELL_COUNTS_ALL <- META_LOCKED[
  ,
  .N,
  by = cell_type
]

CELL_COUNTS_ALL[
  ,
  stratum := "ALL"
]

CELL_COUNTS_AF <- META_LOCKED[
  rhythm == "AF",
  .N,
  by = cell_type
]

CELL_COUNTS_AF[
  ,
  stratum := "AF"
]

CELL_COUNTS_SR <- META_LOCKED[
  rhythm == "SR",
  .N,
  by = cell_type
]

CELL_COUNTS_SR[
  ,
  stratum := "SR"
]

CELLTYPE_COUNTS <- rbindlist(
  list(
    CELL_COUNTS_ALL,
    CELL_COUNTS_AF,
    CELL_COUNTS_SR
  ),
  fill = TRUE
)

setcolorder(
  CELLTYPE_COUNTS,
  c(
    "stratum",
    "cell_type",
    "N"
  )
)

fwrite(
  CELLTYPE_COUNTS,
  file.path(
    QC_DIR,
    "STEP11B_celltype_cell_counts_ALL_AF_SR.csv"
  )
)

# ============================================================
# 11. PREALLOCATE AGGREGATE MATRICES
# ============================================================

make_zero_matrix <- function() {

  m <- matrix(
    0,
    nrow = N_PEAKS,
    ncol = length(CELL_TYPES)
  )

  colnames(
    m
  ) <- CELL_TYPES

  m
}

DET_ALL <- make_zero_matrix()
CNT_ALL <- make_zero_matrix()

DET_AF <- make_zero_matrix()
CNT_AF <- make_zero_matrix()

DET_SR <- make_zero_matrix()
CNT_SR <- make_zero_matrix()

# ============================================================
# 12. STREAM EACH ATAC FILE AND AGGREGATE BY AUTHOR CELL TYPE
# ============================================================

ATAC_STREAM_QC_LIST <- list()

for (d in EXPECTED_DONORS) {

  cat(
    "\n============================================================\n",
    "Streaming ATAC donor: ",
    d,
    "\n",
    "============================================================\n",
    sep = ""
  )

  atac_file <- get_one_file(
    d,
    "ATAC"
  )

  barcodes <- read_header_barcodes(
    atac_file
  )

  meta_d <- META_LOCKED[
    donor == d
  ]

  idx_meta <- match(
    barcodes,
    meta_d$barcode
  )

  if (anyNA(idx_meta)) {
    stop(
      "Metadata alignment failed during ATAC streaming for ",
      d
    )
  }

  meta_ordered <- meta_d[
    idx_meta
  ]

  # cell x cell-type indicator
  indicator <- matrix(
    0,
    nrow = length(barcodes),
    ncol = length(CELL_TYPES)
  )

  colnames(
    indicator
  ) <- CELL_TYPES

  ct_index <- match(
    meta_ordered$cell_type,
    CELL_TYPES
  )

  indicator[
    cbind(
      seq_along(ct_index),
      ct_index
    )
  ] <- 1

  rhythm_d <- unique(
    meta_ordered$rhythm
  )

  if (length(rhythm_d) != 1L) {
    stop(
      "Donor ",
      d,
      " has more than one Rhythm label."
    )
  }

  con <- gzfile(
    atac_file,
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

  row_start <- 1L
  batch_id <- 0L
  n_rows_seen <- 0L
  peak_order_match <- TRUE

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

    row_end <- row_start + n_batch - 1L

    if (row_end > N_PEAKS) {
      close(con)
      stop(
        "ATAC file has more rows than expected for donor ",
        d
      )
    }

    peak_batch <- as.character(
      DT[[1]]
    )

    expected_peak_batch <- PEAKS[
      row_start:row_end
    ]

    if (!identical(
      peak_batch,
      expected_peak_batch
    )) {
      peak_order_match <- FALSE
      close(con)
      stop(
        "ATAC peak order mismatch for donor ",
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

    if (ncol(X) != length(barcodes)) {
      close(con)
      stop(
        "ATAC cell-column count mismatch for donor ",
        d,
        " at batch ",
        batch_id
      )
    }

    DET_BLOCK <- (X > 0) %*% indicator
    CNT_BLOCK <- X %*% indicator

    DET_ALL[
      row_start:row_end,
    ] <- DET_ALL[
      row_start:row_end,
    ] + DET_BLOCK

    CNT_ALL[
      row_start:row_end,
    ] <- CNT_ALL[
      row_start:row_end,
    ] + CNT_BLOCK

    if (rhythm_d == "AF") {

      DET_AF[
        row_start:row_end,
      ] <- DET_AF[
        row_start:row_end,
      ] + DET_BLOCK

      CNT_AF[
        row_start:row_end,
      ] <- CNT_AF[
        row_start:row_end,
      ] + CNT_BLOCK

    } else if (rhythm_d == "SR") {

      DET_SR[
        row_start:row_end,
      ] <- DET_SR[
        row_start:row_end,
      ] + DET_BLOCK

      CNT_SR[
        row_start:row_end,
      ] <- CNT_SR[
        row_start:row_end,
      ] + CNT_BLOCK
    }

    n_rows_seen <- n_rows_seen + n_batch
    row_start <- row_end + 1L

    rm(
      DT,
      X,
      DET_BLOCK,
      CNT_BLOCK
    )

    if (
      batch_id %% 25L == 0L
    ) {

      cat(
        "  processed peaks: ",
        format(
          n_rows_seen,
          big.mark = ","
        ),
        " / ",
        format(
          N_PEAKS,
          big.mark = ","
        ),
        "\n",
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

  if (n_rows_seen != N_PEAKS) {
    stop(
      "ATAC peak count mismatch for donor ",
      d,
      ": observed ",
      n_rows_seen,
      ", expected ",
      N_PEAKS
    )
  }

  ATAC_STREAM_QC_LIST[[
    length(ATAC_STREAM_QC_LIST) + 1L
  ]] <- data.table(
    donor = d,
    rhythm = rhythm_d,
    n_cells = length(barcodes),
    n_peaks = n_rows_seen,
    peak_order_matches_reference = peak_order_match,
    n_batches = batch_id
  )

  cat(
    "Completed ",
    d,
    ": ",
    format(n_rows_seen, big.mark = ","),
    " peaks x ",
    length(barcodes),
    " cells.\n",
    sep = ""
  )
}

ATAC_STREAM_QC <- rbindlist(
  ATAC_STREAM_QC_LIST
)

fwrite(
  ATAC_STREAM_QC,
  file.path(
    QC_DIR,
    "STEP11B_ATAC_streaming_QC.csv"
  )
)

# ============================================================
# 13. SAVE AGGREGATES
# ============================================================

make_agg_object <- function(
  det,
  cnt,
  stratum
) {

  stratum_value <- stratum

  n_cells_vec <- setNames(
    rep(
      0L,
      length(CELL_TYPES)
    ),
    CELL_TYPES
  )

  tmp <- CELLTYPE_COUNTS[
    stratum == stratum_value
  ]

  n_cells_vec[
    tmp$cell_type
  ] <- tmp$N

  list(
    source = "GSE238242",
    tissue = "Left atrial appendage",
    genome_build = "GRCh38",
    author_cell_types = CELL_TYPES,
    stratum = stratum,
    peaks = PEAKS,
    n_cells_by_celltype = n_cells_vec,
    accessible_cell_counts = det,
    total_ATAC_counts = cnt
  )
}

AGG_ALL <- make_agg_object(
  DET_ALL,
  CNT_ALL,
  "ALL"
)

AGG_AF <- make_agg_object(
  DET_AF,
  CNT_AF,
  "AF"
)

AGG_SR <- make_agg_object(
  DET_SR,
  CNT_SR,
  "SR"
)

saveRDS(
  AGG_ALL,
  file.path(
    AGG_DIR,
    "GSE238242_ATAC_celltype_aggregate_ALL_GRCh38.rds"
  ),
  compress = TRUE
)

saveRDS(
  AGG_AF,
  file.path(
    AGG_DIR,
    "GSE238242_ATAC_celltype_aggregate_AF_GRCh38.rds"
  ),
  compress = TRUE
)

saveRDS(
  AGG_SR,
  file.path(
    AGG_DIR,
    "GSE238242_ATAC_celltype_aggregate_SR_GRCh38.rds"
  ),
  compress = TRUE
)

# ============================================================
# 14. PEAK-PREVALENCE QC AT MULTIPLE THRESHOLDS
#     QC ONLY — DOES NOT CHOOSE FINAL S-LDSC ANNOTATION.
# ============================================================

summarize_prevalence <- function(
  agg,
  stratum
) {

  stratum_value <- stratum
  out <- list()

  for (ct in CELL_TYPES) {

    n_cells <- agg$n_cells_by_celltype[
      ct
    ]

    det <- agg$accessible_cell_counts[
      ,
      ct
    ]

    if (
      is.na(n_cells) ||
      n_cells <= 0L
    ) {

      out[[
        length(out) + 1L
      ]] <- data.table(
        stratum = stratum_value,
        cell_type = ct,
        n_cells = 0L,
        threshold_rule = c(
          ">=1_cell",
          ">=max(5,1%)",
          ">=max(5,5%)",
          ">=max(5,10%)"
        ),
        required_cells = NA_integer_,
        n_peaks_passing = NA_integer_
      )

      next
    }

    rules <- data.table(
      threshold_rule = c(
        ">=1_cell",
        ">=max(5,1%)",
        ">=max(5,5%)",
        ">=max(5,10%)"
      ),
      required_cells = c(
        1L,
        max(
          5L,
          ceiling(
            0.01 * n_cells
          )
        ),
        max(
          5L,
          ceiling(
            0.05 * n_cells
          )
        ),
        max(
          5L,
          ceiling(
            0.10 * n_cells
          )
        )
      )
    )

    rules[
      ,
      n_peaks_passing :=
        vapply(
          required_cells,
          function(k) {
            sum(
              det >= k
            )
          },
          integer(1)
        )
    ]

    rules[
      ,
      `:=`(
        stratum = stratum_value,
        cell_type = ct,
        n_cells = n_cells
      )
    ]

    setcolorder(
      rules,
      c(
        "stratum",
        "cell_type",
        "n_cells",
        "threshold_rule",
        "required_cells",
        "n_peaks_passing"
      )
    )

    out[[
      length(out) + 1L
    ]] <- rules
  }

  rbindlist(
    out,
    fill = TRUE
  )
}

PREV_QC <- rbindlist(
  list(
    summarize_prevalence(
      AGG_ALL,
      "ALL"
    ),
    summarize_prevalence(
      AGG_AF,
      "AF"
    ),
    summarize_prevalence(
      AGG_SR,
      "SR"
    )
  ),
  fill = TRUE
)

fwrite(
  PREV_QC,
  file.path(
    QC_DIR,
    "STEP11B_peak_prevalence_threshold_QC.csv"
  )
)

# ============================================================
# 15. GLOBAL AGGREGATE CONSISTENCY CHECKS
# ============================================================

# Detection counts can never exceed cell-type cell counts.
check_det_limits <- function(
  agg
) {

  ok <- TRUE

  for (ct in CELL_TYPES) {

    n_cells <- agg$n_cells_by_celltype[
      ct
    ]

    if (
      is.na(n_cells) ||
      n_cells == 0L
    ) {
      next
    }

    if (
      any(
        agg$accessible_cell_counts[
          ,
          ct
        ] > n_cells
      )
    ) {
      ok <- FALSE
      break
    }
  }

  ok
}

CONSISTENCY_QC <- data.table(
  check = c(
    "ALL_detection_counts_within_cell_counts",
    "AF_detection_counts_within_cell_counts",
    "SR_detection_counts_within_cell_counts",
    "ALL_equals_AF_plus_SR_detection",
    "ALL_equals_AF_plus_SR_total_counts"
  ),
  pass = c(
    check_det_limits(AGG_ALL),
    check_det_limits(AGG_AF),
    check_det_limits(AGG_SR),
    isTRUE(
      all.equal(
        DET_ALL,
        DET_AF + DET_SR,
        tolerance = 0
      )
    ),
    isTRUE(
      all.equal(
        CNT_ALL,
        CNT_AF + CNT_SR,
        tolerance = 0
      )
    )
  )
)

fwrite(
  CONSISTENCY_QC,
  file.path(
    QC_DIR,
    "STEP11B_aggregate_consistency_QC.csv"
  )
)

if (!all(CONSISTENCY_QC$pass)) {
  stop(
    "Aggregate consistency QC failed. Inspect STEP11B_aggregate_consistency_QC.csv."
  )
}

# ============================================================
# 16. FINAL READINESS
# ============================================================

READINESS <- data.table(
  check = c(
    "7_donors_present",
    "11986_author_annotated_cells",
    "RNA_ATAC_barcodes_exactly_paired_all_donors",
    "metadata_barcodes_match_matrices_all_donors",
    "RNA_gene_universe_identical_all_donors",
    "ATAC_peak_universe_identical_all_donors",
    "ATAC_peak_coordinates_valid",
    "12_author_cell_types_present",
    "ALL_AF_SR_aggregates_consistent",
    "S_LDSC_peak_threshold_not_yet_selected"
  ),
  pass = c(
    uniqueN(META_LOCKED$donor) == 7L,
    nrow(META_LOCKED) == 11986L,
    all(PAIR_QC$RNA_ATAC_exact_same_order),
    all(
      PAIR_QC$RNA_metadata_set_equal &
      PAIR_QC$ATAC_metadata_set_equal
    ),
    all(
      RNA_GENE_QC$identical_gene_order_to_reference
    ),
    all(
      ATAC_STREAM_QC$peak_order_matches_reference
    ),
    PEAK_COORD_QC$n_invalid == 0L,
    length(CELL_TYPES) == 12L,
    all(CONSISTENCY_QC$pass),
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11B_readiness.csv"
  )
)

# ============================================================
# 17. METHOD / PROVENANCE
# ============================================================

METHOD <- data.table(
  field = c(
    "dataset",
    "tissue",
    "source_build",
    "n_donors",
    "SR_donors",
    "AF_donors",
    "n_cells",
    "n_author_cell_types",
    "n_RNA_genes",
    "n_ATAC_peaks",
    "primary_S_LDSC_stratum",
    "AF_SR_strata_role",
    "reclustering",
    "reannotation",
    "peak_threshold_status"
  ),
  value = c(
    "GSE238242",
    "Left atrial appendage",
    "GRCh38",
    "7",
    "CF69,CF77,CF89,CF91",
    "CF93,CF97,CF102",
    as.character(nrow(META_LOCKED)),
    as.character(length(CELL_TYPES)),
    as.character(length(RNA_GENE_REFERENCE)),
    as.character(N_PEAKS),
    "ALL donors pooled within author cell type",
    "Sensitivity only",
    "No",
    "No; author cell_type labels frozen",
    "Not selected in STEP11B"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11B_method_provenance.csv"
  )
)

# ============================================================
# 18. README / WHAT TO SEND BACK
# ============================================================

README <- c(
  "STEP11B — GSE238242 paired-barcode lock + cell-type ATAC aggregation",
  "",
  "Primary purpose:",
  "Prepare author-cell-type ATAC aggregates for S-LDSC without re-clustering.",
  "",
  "Primary S-LDSC annotation source:",
  "ALL 7 donors pooled within each author-defined cell type.",
  "",
  "AF-only and SR-only aggregates are sensitivity resources only.",
  "",
  "No peak prevalence threshold is frozen in STEP11B.",
  "STEP11C will inspect these QC distributions, construct binary cell-type annotations,",
  "perform GRCh38 -> GRCh37 liftover, and generate LDSC annotation files.",
  "",
  "Please send back:",
  "00_QC/STEP11B_readiness.csv",
  "00_QC/STEP11B_barcode_pairing_QC.csv",
  "00_QC/STEP11B_RNA_gene_universe_QC.csv",
  "00_QC/STEP11B_ATAC_streaming_QC.csv",
  "00_QC/STEP11B_peak_coordinate_QC.csv",
  "00_QC/STEP11B_celltype_cell_counts_ALL_AF_SR.csv",
  "00_QC/STEP11B_peak_prevalence_threshold_QC.csv",
  "00_QC/STEP11B_aggregate_consistency_QC.csv",
  "00_QC/STEP11B_method_provenance.csv",
  "",
  "Also include the 03_SLDSC/00_ATAC_AGGREGATES folder if practical."
)

writeLines(
  README,
  file.path(
    STEP11_ROOT,
    "README_STEP11B.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11B_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP11B COMPLETE\n",
  "============================================================\n",
  sep = ""
)

print(
  READINESS
)

cat(
  "\nAuthor cell types:\n"
)

print(
  CELLTYPE_COUNTS[
    stratum == "ALL"
  ][
    order(
      -N
    )
  ]
)

cat(
  "\nSend back the QC files listed in README_STEP11B.txt.\n",
  "Do not start STEP11C until STEP11B QC is reviewed.\n",
  sep = ""
)
