# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11D_GSE238242_LIFTOVER_GRCh38_TO_hg19_FINAL.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11D — GSE238242 ATAC annotation liftOver
# GRCh38 -> GRCh37/hg19
#
# INPUT
#   D:/A/data/STEP11_GSE238242/03_SLDSC/01_ANNOTATION_GRCh38/
#     GSE238242_celltype_binary_annotations_GRCh38.rds
#
# CHAIN
#   D:/A/data/hg38ToHg19.over.chain
#
# OUTPUT
#   D:/A/data/STEP11_GSE238242/03_SLDSC/02_LIFTOVER_GRCh37/
#
# DESIGN
#   - Lift the shared 212,084-peak universe ONCE.
#   - Keep only uniquely mapped intervals for primary analysis.
#   - Keep chr1-22 only for S-LDSC.
#   - Apply the frozen STEP11C binary annotations after mapping.
#   - Export both raw uniquely-mapped intervals and reduced/merged BEDs.
#
# IMPORTANT
#   - No S-LDSC is run in this step.
#   - No threshold is changed.
#   - No AF/SR redefinition is performed.
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

INPUT_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "01_ANNOTATION_GRCh38"
)

OUTPUT_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "02_LIFTOVER_GRCh37"
)

QC_DIR <- file.path(
  STEP11_ROOT,
  "00_QC"
)

ANNOT_RDS <- file.path(
  INPUT_DIR,
  "GSE238242_celltype_binary_annotations_GRCh38.rds"
)

CHAIN_FILE <- file.path(
  DATA_ROOT,
  "hg38ToHg19.over.chain"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "GenomicRanges",
  "IRanges",
  "rtracklayer"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall them first, then rerun STEP11D."
    )
  )
}

library(data.table)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)

# ============================================================
# 2. HARD INPUT CHECK
# ============================================================

required_files <- c(
  annotation_object = ANNOT_RDS,
  liftover_chain = CHAIN_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files)) {
  stop(
    paste0(
      "Missing STEP11D input(s):\n",
      paste(
        names(missing_files),
        missing_files,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

# ============================================================
# 3. LOAD FROZEN STEP11C ANNOTATION OBJECT
# ============================================================

ANNOT_OBJECT <- readRDS(
  ANNOT_RDS
)

required_names <- c(
  "genome_build",
  "cell_types",
  "peak_ids",
  "binary_annotation",
  "cell_numbers"
)

if (!all(required_names %in% names(ANNOT_OBJECT))) {
  stop(
    "STEP11C annotation RDS does not contain the expected fields."
  )
}

if (!identical(
  as.character(ANNOT_OBJECT$genome_build),
  "GRCh38"
)) {
  stop(
    "Input annotation object is not labelled GRCh38."
  )
}

PEAKS <- as.character(
  ANNOT_OBJECT$peak_ids
)

ANNOT <- ANNOT_OBJECT$binary_annotation

CELL_TYPES <- as.character(
  ANNOT_OBJECT$cell_types
)

if (nrow(ANNOT) != length(PEAKS)) {
  stop(
    "Peak count and annotation-row count are inconsistent."
  )
}

if (ncol(ANNOT) != length(CELL_TYPES)) {
  stop(
    "Cell-type count and annotation-column count are inconsistent."
  )
}

if (is.null(colnames(ANNOT))) {
  colnames(ANNOT) <- CELL_TYPES
}

if (!all(CELL_TYPES %in% colnames(ANNOT))) {
  stop(
    "Annotation matrix does not contain all expected cell-type columns."
  )
}

# ============================================================
# 4. PARSE GRCh38 PEAK IDS
#    Expected format: chr-start-end
# ============================================================

parts <- tstrsplit(
  PEAKS,
  "-",
  fixed = TRUE
)

if (length(parts) != 3L) {
  stop(
    "Peak IDs are not consistently formatted as chr-start-end."
  )
}

PEAK_DT <- data.table(
  original_index = seq_along(PEAKS),
  peak_id_GRCh38 = PEAKS,
  chr38 = parts[[1]],
  start38 = suppressWarnings(
    as.integer(parts[[2]])
  ),
  end38 = suppressWarnings(
    as.integer(parts[[3]])
  )
)

PEAK_DT[
  ,
  valid38 :=
    !is.na(start38) &
    !is.na(end38) &
    start38 >= 1L &
    end38 >= start38 &
    grepl(
      "^chr([0-9]+|X|Y|M)$",
      chr38
    )
]

if (!all(PEAK_DT$valid38)) {
  stop(
    "Invalid GRCh38 peak coordinates detected: ",
    sum(!PEAK_DT$valid38)
  )
}

# ============================================================
# 5. IMPORT VERIFIED CHAIN
# ============================================================

CHAIN <- rtracklayer::import.chain(
  CHAIN_FILE
)

cat(
  "\nChain imported successfully.\n",
  "Chain length: ",
  length(CHAIN),
  "\n",
  sep = ""
)

# ============================================================
# 6. LIFT THE SHARED PEAK UNIVERSE ONCE
# ============================================================

GR38 <- GRanges(
  seqnames = PEAK_DT$chr38,
  ranges = IRanges(
    start = PEAK_DT$start38,
    end = PEAK_DT$end38
  )
)

names(GR38) <- as.character(
  PEAK_DT$original_index
)

cat(
  "\nRunning liftOver for ",
  format(length(GR38), big.mark = ","),
  " shared ATAC peaks...\n",
  sep = ""
)

LIFTED <- rtracklayer::liftOver(
  GR38,
  CHAIN
)

N_MAP <- lengths(
  LIFTED
)

if (length(N_MAP) != length(PEAKS)) {
  stop(
    "liftOver output length does not match input peak universe."
  )
}

MAP_CLASS <- fifelse(
  N_MAP == 0L,
  "unmapped",
  fifelse(
    N_MAP == 1L,
    "unique",
    "multimapped"
  )
)

UNIQUE_INDEX <- which(
  N_MAP == 1L
)

if (!length(UNIQUE_INDEX)) {
  stop(
    "No peaks mapped uniquely to hg19."
  )
}

GR19_UNIQUE <- unlist(
  LIFTED[
    UNIQUE_INDEX
  ],
  use.names = FALSE
)

if (length(GR19_UNIQUE) != length(UNIQUE_INDEX)) {
  stop(
    "Unexpected mismatch between unique input peaks and lifted ranges."
  )
}

MAP_DT <- data.table(
  original_index = UNIQUE_INDEX,
  peak_id_GRCh38 = PEAKS[UNIQUE_INDEX],
  chr38 = PEAK_DT$chr38[UNIQUE_INDEX],
  start38 = PEAK_DT$start38[UNIQUE_INDEX],
  end38 = PEAK_DT$end38[UNIQUE_INDEX],
  chr19 = as.character(
    seqnames(GR19_UNIQUE)
  ),
  start19 = start(
    GR19_UNIQUE
  ),
  end19 = end(
    GR19_UNIQUE
  )
)

MAP_DT[
  ,
  autosomal19 :=
    grepl(
      "^chr([1-9]|1[0-9]|2[0-2])$",
      chr19
    )
]

# Save all unique mappings (including X/Y/M if present).
fwrite(
  MAP_DT,
  file.path(
    OUTPUT_DIR,
    "STEP11D_peak_universe_unique_liftover_GRCh38_to_hg19.csv.gz"
  )
)

# Save the map status for every input peak.
MAP_STATUS <- data.table(
  original_index = seq_along(PEAKS),
  peak_id_GRCh38 = PEAKS,
  n_hg19_mappings = N_MAP,
  mapping_class = MAP_CLASS
)

fwrite(
  MAP_STATUS,
  file.path(
    OUTPUT_DIR,
    "STEP11D_peak_universe_mapping_status.csv.gz"
  )
)

# ============================================================
# 7. GLOBAL LIFTOVER QC
# ============================================================

GLOBAL_QC <- data.table(
  metric = c(
    "input_peaks",
    "unmapped_peaks",
    "unique_mapped_peaks",
    "multimapped_peaks",
    "unique_mapping_rate",
    "unique_autosomal_hg19_peaks",
    "unique_autosomal_rate_among_input"
  ),
  value = c(
    length(PEAKS),
    sum(N_MAP == 0L),
    sum(N_MAP == 1L),
    sum(N_MAP > 1L),
    mean(N_MAP == 1L),
    sum(MAP_DT$autosomal19),
    sum(MAP_DT$autosomal19) / length(PEAKS)
  )
)

fwrite(
  GLOBAL_QC,
  file.path(
    OUTPUT_DIR,
    "STEP11D_global_liftover_QC.csv"
  )
)

# ============================================================
# 8. APPLY EACH CELL-TYPE ANNOTATION
# ============================================================

CELLTYPE_QC_LIST <- list()

for (ct in CELL_TYPES) {

  cat(
    "\nProcessing cell type: ",
    ct,
    "\n",
    sep = ""
  )

  anno_vec <- as.integer(
    ANNOT[
      ,
      ct
    ]
  )

  if (!all(anno_vec %in% c(0L, 1L))) {
    stop(
      "Non-binary annotation values detected for ",
      ct
    )
  }

  input_annot_idx <- which(
    anno_vec == 1L
  )

  # Retain only uniquely mapped + autosomal peaks.
  ct_map <- MAP_DT[
    original_index %in% input_annot_idx &
      autosomal19 == TRUE
  ]

  # Raw uniquely-mapped intervals, preserving the original peak identity.
  RAW_OUT <- data.table(
    chr = ct_map$chr19,
    start_1based = ct_map$start19,
    end_1based = ct_map$end19,
    peak_id_GRCh38 = ct_map$peak_id_GRCh38,
    original_index = ct_map$original_index,
    annotation = 1L
  )

  fwrite(
    RAW_OUT,
    file.path(
      OUTPUT_DIR,
      paste0(
        "GSE238242_",
        ct,
        "_ATAC_hg19_unique_autosomal.tsv.gz"
      )
    ),
    sep = "\t"
  )

  # Reduce overlapping intervals.
  if (nrow(ct_map) > 0L) {

    ct_gr <- GRanges(
      seqnames = ct_map$chr19,
      ranges = IRanges(
        start = ct_map$start19,
        end = ct_map$end19
      )
    )

    ct_gr_reduced <- GenomicRanges::reduce(
      ct_gr,
      ignore.strand = TRUE
    )

    MERGED_OUT <- data.table(
      chr = as.character(
        seqnames(ct_gr_reduced)
      ),
      start_1based = start(
        ct_gr_reduced
      ),
      end_1based = end(
        ct_gr_reduced
      ),
      annotation = 1L
    )

    # BED version for downstream interval tools:
    # BED start is 0-based; BED end is end-exclusive.
    BED_OUT <- data.table(
      chr = MERGED_OUT$chr,
      start_0based = pmax(
        0L,
        MERGED_OUT$start_1based - 1L
      ),
      end_0based_exclusive = MERGED_OUT$end_1based,
      name = paste0(
        "GSE238242_",
        ct,
        "_ATAC"
      )
    )

  } else {

    MERGED_OUT <- data.table(
      chr = character(),
      start_1based = integer(),
      end_1based = integer(),
      annotation = integer()
    )

    BED_OUT <- data.table(
      chr = character(),
      start_0based = integer(),
      end_0based_exclusive = integer(),
      name = character()
    )
  }

  fwrite(
    MERGED_OUT,
    file.path(
      OUTPUT_DIR,
      paste0(
        "GSE238242_",
        ct,
        "_ATAC_hg19_merged.tsv.gz"
      )
    ),
    sep = "\t"
  )

  fwrite(
    BED_OUT,
    file.path(
      OUTPUT_DIR,
      paste0(
        "GSE238242_",
        ct,
        "_ATAC_hg19_merged.bed.gz"
      )
    ),
    sep = "\t",
    col.names = FALSE
  )

  input_n <- length(
    input_annot_idx
  )

  unique_n <- sum(
    N_MAP[
      input_annot_idx
    ] == 1L
  )

  multi_n <- sum(
    N_MAP[
      input_annot_idx
    ] > 1L
  )

  unmapped_n <- sum(
    N_MAP[
      input_annot_idx
    ] == 0L
  )

  autosomal_n <- nrow(
    ct_map
  )

  CELLTYPE_QC_LIST[[
    length(CELLTYPE_QC_LIST) + 1L
  ]] <- data.table(
    cell_type = ct,
    input_annotated_peaks_GRCh38 = input_n,
    unmapped = unmapped_n,
    multimapped = multi_n,
    unique_mapped = unique_n,
    unique_mapping_rate = if (input_n > 0L) unique_n / input_n else NA_real_,
    unique_autosomal_hg19 = autosomal_n,
    autosomal_retention_rate = if (input_n > 0L) autosomal_n / input_n else NA_real_,
    merged_hg19_intervals = nrow(
      MERGED_OUT
    )
  )
}

CELLTYPE_QC <- rbindlist(
  CELLTYPE_QC_LIST,
  fill = TRUE
)

fwrite(
  CELLTYPE_QC,
  file.path(
    OUTPUT_DIR,
    "STEP11D_celltype_liftover_QC.csv"
  )
)

# ============================================================
# 9. READINESS RULES
# ============================================================

# Conservative QC thresholds:
#   - >=90% unique mapping for every cell type
#   - >=90% autosomal retention for every cell type
# We do not fail automatically for lower values; instead readiness records it.

READINESS <- data.table(
  check = c(
    "chain_imported",
    "shared_peak_universe_lifted",
    "all_12_celltypes_exported",
    "all_celltypes_unique_mapping_rate_ge_0.90",
    "all_celltypes_autosomal_retention_rate_ge_0.90",
    "hg19_BED_files_created"
  ),
  pass = c(
    length(CHAIN) > 0L,
    nrow(MAP_STATUS) == length(PEAKS),
    nrow(CELLTYPE_QC) == length(CELL_TYPES),
    all(
      CELLTYPE_QC$unique_mapping_rate >= 0.90,
      na.rm = TRUE
    ),
    all(
      CELLTYPE_QC$autosomal_retention_rate >= 0.90,
      na.rm = TRUE
    ),
    all(
      file.exists(
        file.path(
          OUTPUT_DIR,
          paste0(
            "GSE238242_",
            CELL_TYPES,
            "_ATAC_hg19_merged.bed.gz"
          )
        )
      )
    )
  )
)

fwrite(
  READINESS,
  file.path(
    OUTPUT_DIR,
    "STEP11D_readiness.csv"
  )
)

# ============================================================
# 10. PROVENANCE
# ============================================================

PROVENANCE <- data.table(
  field = c(
    "dataset",
    "tissue",
    "source_build",
    "target_build",
    "chain_file",
    "source_annotation_threshold",
    "liftover_rule",
    "multimapping_rule",
    "chromosomes_used_for_S_LDSC",
    "interval_export"
  ),
  value = c(
    "GSE238242",
    "Left atrial appendage",
    "GRCh38",
    "GRCh37/hg19",
    normalizePath(
      CHAIN_FILE,
      winslash = "/",
      mustWork = TRUE
    ),
    as.character(
      ANNOT_OBJECT$primary_threshold
    ),
    "lift shared peak universe once; apply cell-type annotations after mapping",
    "retain uniquely mapped peaks only",
    "chr1-22",
    "1-based TSV plus 0-based BED"
  )
)

fwrite(
  PROVENANCE,
  file.path(
    OUTPUT_DIR,
    "STEP11D_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTPUT_DIR,
    "STEP11D_sessionInfo.txt"
  )
)

# ============================================================
# 11. CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================================\n",
  "STEP11D COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Global liftOver QC:\n"
)

print(
  GLOBAL_QC
)

cat(
  "\nCell-type liftOver QC:\n"
)

print(
  CELLTYPE_QC
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

cat(
  "\nPlease send back these files:\n",
  "  STEP11D_global_liftover_QC.csv\n",
  "  STEP11D_celltype_liftover_QC.csv\n",
  "  STEP11D_readiness.csv\n",
  "  STEP11D_method_provenance.csv\n",
  "\nOutput folder:\n",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
