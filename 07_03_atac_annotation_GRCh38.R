# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11C_GSE238242_ATAC_ANNOTATION_GRCh38.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11C — GSE238242 cell-type ATAC annotation construction
# Purpose:
#   Build GRCh38 cell-type-specific binary ATAC annotations
#   from STEP11B aggregates.
#
# Strategy:
#   Primary: ALL donors pooled.
#   Threshold: peak accessible in >=10% of cells of a cell type.
#   Sensitivity thresholds (5%,1%) are exported but NOT primary.
#
# Output:
#   - peak annotation tables
#   - cell-type peak summary
#   - GRCh38 annotation objects for later liftover/LDSC
# ============================================================

rm(list=ls())
options(stringsAsFactors=FALSE, scipen=999)

if(!requireNamespace("data.table", quietly=TRUE)){
  install.packages("data.table",
                   repos="https://cloud.r-project.org")
}

library(data.table)

DATA_ROOT <- "D:/A/data"
STEP11_ROOT <- file.path(DATA_ROOT,"STEP11_GSE238242")

AGG_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "00_ATAC_AGGREGATES"
)

OUT_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "01_ANNOTATION_GRCh38"
)

QC_DIR <- file.path(
  STEP11_ROOT,
  "00_QC"
)

dir.create(
  OUT_DIR,
  recursive=TRUE,
  showWarnings=FALSE
)

# ============================================================
# 1. Load ALL donor aggregate
# ============================================================

agg_file <- file.path(
  AGG_DIR,
  "GSE238242_ATAC_celltype_aggregate_ALL_GRCh38.rds"
)

if(!file.exists(agg_file)){
  stop("Missing aggregate file: ", agg_file)
}

AGG <- readRDS(agg_file)

PEAKS <- AGG$peaks
CELL_TYPES <- AGG$author_cell_types

DET <- AGG$accessible_cell_counts

N_CELLS <- AGG$n_cells_by_celltype

# ============================================================
# 2. Construct binary annotations
# ============================================================

PRIMARY_THRESHOLD <- 0.10

ANNOT <- matrix(
  0L,
  nrow=nrow(DET),
  ncol=ncol(DET)
)

colnames(ANNOT) <- CELL_TYPES
rownames(ANNOT) <- PEAKS

SUMMARY <- list()

for(ct in CELL_TYPES){

  ncell <- N_CELLS[ct]

  cutoff <- ceiling(
    ncell * PRIMARY_THRESHOLD
  )

  # minimum 5 cells
  cutoff <- max(
    cutoff,
    5L
  )

  keep <- DET[,ct] >= cutoff

  ANNOT[,ct] <- as.integer(keep)

  SUMMARY[[length(SUMMARY)+1]] <- data.table(
    cell_type=ct,
    n_cells=ncell,
    threshold=">=10%",
    required_accessible_cells=cutoff,
    n_annotated_peaks=sum(keep)
  )
}

SUMMARY <- rbindlist(SUMMARY)

fwrite(
  SUMMARY,
  file.path(
    OUT_DIR,
    "STEP11C_primary_annotation_summary.csv"
  )
)

# ============================================================
# 3. Export each cell type annotation
# ============================================================

for(ct in CELL_TYPES){

  out <- data.table(
    peak=PEAKS,
    chr=sub("-.*","",PEAKS),
    start=as.integer(sub(
      ".*-|-[^-]*$",
      "",
      PEAKS
    )),
    end=as.integer(sub(
      ".*-",
      "",
      PEAKS
    )),
    annotation=ANNOT[,ct]
  )

  fwrite(
    out,
    file.path(
      OUT_DIR,
      paste0(
        "GSE238242_",
        ct,
        "_ATAC_annotation_GRCh38.tsv.gz"
      )
    ),
    sep="\t"
  )
}

# ============================================================
# 4. Save combined object
# ============================================================

ANNOT_OBJECT <- list(
  dataset="GSE238242",
  genome_build="GRCh38",
  tissue="Left atrial appendage",
  primary_threshold=">=10% accessible cells",
  cell_types=CELL_TYPES,
  peak_ids=PEAKS,
  binary_annotation=ANNOT,
  cell_numbers=N_CELLS
)

saveRDS(
  ANNOT_OBJECT,
  file.path(
    OUT_DIR,
    "GSE238242_celltype_binary_annotations_GRCh38.rds"
  ),
  compress=TRUE
)

# ============================================================
# 5. Export QC matrix
# ============================================================

QC_MATRIX <- data.table(
  peak=PEAKS,
  ANNOT
)

fwrite(
  QC_MATRIX,
  file.path(
    OUT_DIR,
    "GSE238242_annotation_binary_matrix_GRCh38.csv.gz"
  )
)

# ============================================================
# 6. Sensitivity threshold export
# ============================================================

SENS <- list()

for(thr in c(0.01,0.05,0.10)){

  tmp <- matrix(
    0L,
    nrow=nrow(DET),
    ncol=ncol(DET)
  )

  colnames(tmp) <- CELL_TYPES

  for(ct in CELL_TYPES){

    cutoff <- max(
      5L,
      ceiling(
        N_CELLS[ct]*thr
      )
    )

    tmp[,ct] <- as.integer(
      DET[,ct] >= cutoff
    )
  }

  SENS[[paste0(thr*100,"%")]] <- tmp
}

saveRDS(
  SENS,
  file.path(
    OUT_DIR,
    "GSE238242_annotation_sensitivity_thresholds.rds"
  ),
  compress=TRUE
)

cat("\nSTEP11C COMPLETE\n")
print(SUMMARY)
