# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11A_GSE238242_LOCAL_AUDIT.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11A — GSE238242 LAA sn-multiome local audit
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# PURPOSE
#   1) Use the already-downloaded local GSE238242 files
#   2) Extract GSE238242_RAW.tar without modifying source files
#   3) Audit the 14 expected RNA/ATAC sample files
#   4) Audit author-provided cell metadata and original labels
#   5) Freeze the exact donor / rhythm / modality design for STEP11B
#
# IMPORTANT
#   - NO re-clustering
#   - NO cell filtering
#   - NO marker-based re-annotation
#   - NO S-LDSC / SCAVENGE yet
#   - Source files under D:/A/data remain untouched
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"

RAW_TAR <- file.path(
  DATA_ROOT,
  "GSE238242_RAW.tar"
)

META_GZ <- file.path(
  DATA_ROOT,
  "GSE238242_snAF.metadata.tsv.gz"
)

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

EXTRACT_DIR <- file.path(
  OUT_ROOT,
  "01_EXTRACTED"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "00_QC"
)

OBJECT_DIR <- file.path(
  OUT_ROOT,
  "02_OBJECTS"
)

SLDSC_DIR <- file.path(
  OUT_ROOT,
  "03_SLDSC"
)

SCAVENGE_DIR <- file.path(
  OUT_ROOT,
  "04_SCAVENGE"
)

FINAL_DIR <- file.path(
  OUT_ROOT,
  "05_FINAL"
)

for (d in c(
  OUT_ROOT,
  EXTRACT_DIR,
  QC_DIR,
  OBJECT_DIR,
  SLDSC_DIR,
  SCAVENGE_DIR,
  FINAL_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 1. PACKAGE
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

required_files <- c(
  RAW_TAR = RAW_TAR,
  META_GZ = META_GZ
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files)) {
  stop(
    paste0(
      "Missing required local file(s):\n",
      paste(
        names(missing_files),
        missing_files,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

INPUT_QC <- data.table(
  item = names(
    required_files
  ),
  path = normalizePath(
    unname(
      required_files
    ),
    winslash = "/",
    mustWork = TRUE
  ),
  size_MB = as.numeric(
    file.info(
      unname(
        required_files
      )
    )$size
  ) / 1024^2
)

fwrite(
  INPUT_QC,
  file.path(
    QC_DIR,
    "STEP11A_input_file_QC.csv"
  )
)

cat(
  "\nLocal GSE238242 inputs found:\n"
)
print(
  INPUT_QC
)

# ============================================================
# 3. FROZEN GEO DESIGN
# ============================================================

EXPECTED_DESIGN <- data.table(
  GSM = c(
    "GSM7660990",
    "GSM7660991",
    "GSM7660992",
    "GSM7660993",
    "GSM7660994",
    "GSM7660995",
    "GSM7660996",
    "GSM7660997",
    "GSM7660998",
    "GSM7660999",
    "GSM7661000",
    "GSM7661001",
    "GSM7661002",
    "GSM7661003"
  ),
  donor = rep(
    c(
      "CF69",
      "CF77",
      "CF89",
      "CF91",
      "CF93",
      "CF97",
      "CF102"
    ),
    times = 2
  ),
  rhythm = rep(
    c(
      "SR",
      "SR",
      "SR",
      "SR",
      "AF",
      "AF",
      "AF"
    ),
    times = 2
  ),
  modality = c(
    rep(
      "RNA",
      7
    ),
    rep(
      "ATAC",
      7
    )
  ),
  tissue = "LAA",
  genome_build = "GRCh38"
)

fwrite(
  EXPECTED_DESIGN,
  file.path(
    QC_DIR,
    "STEP11A_expected_GEO_design.csv"
  )
)

# ============================================================
# 4. LIST TAR MEMBERS BEFORE EXTRACTION
# ============================================================

tar_members <- tryCatch(
  utils::untar(
    RAW_TAR,
    list = TRUE
  ),
  error = function(e) {
    stop(
      "Could not list GSE238242_RAW.tar: ",
      conditionMessage(e)
    )
  }
)

if (!length(tar_members)) {
  stop(
    "GSE238242_RAW.tar contains no members."
  )
}

TAR_MANIFEST <- data.table(
  archive_member = tar_members
)

TAR_MANIFEST[
  ,
  file_name :=
    basename(
      archive_member
    )
]

TAR_MANIFEST[
  ,
  GSM :=
    fifelse(
      grepl(
        "GSM[0-9]+",
        file_name,
        ignore.case = TRUE
      ),
      sub(
        ".*?(GSM[0-9]+).*",
        "\\1",
        file_name,
        ignore.case = TRUE
      ),
      NA_character_
    )
]

TAR_MANIFEST[
  ,
  donor :=
    fifelse(
      grepl(
        "CF[0-9]+",
        file_name,
        ignore.case = TRUE
      ),
      toupper(
        sub(
          ".*?(CF[0-9]+).*",
          "\\1",
          file_name,
          ignore.case = TRUE
        )
      ),
      NA_character_
    )
]

TAR_MANIFEST[
  ,
  modality :=
    fifelse(
      grepl(
        "RNA",
        file_name,
        ignore.case = TRUE
      ),
      "RNA",
      fifelse(
        grepl(
          "ATAC",
          file_name,
          ignore.case = TRUE
        ),
        "ATAC",
        NA_character_
      )
    )
]

fwrite(
  TAR_MANIFEST,
  file.path(
    QC_DIR,
    "STEP11A_RAW_tar_manifest.csv"
  )
)

# ============================================================
# 5. EXTRACT ARCHIVE
# ============================================================

extract_marker <- file.path(
  EXTRACT_DIR,
  ".STEP11A_EXTRACT_COMPLETE"
)

if (!file.exists(extract_marker)) {

  cat(
    "\nExtracting GSE238242_RAW.tar...\n"
  )

  utils::untar(
    RAW_TAR,
    exdir = EXTRACT_DIR
  )

  writeLines(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    extract_marker
  )

} else {

  cat(
    "\nExisting extraction marker detected; skipping extraction.\n"
  )
}

# ============================================================
# 6. EXTRACTED FILE MANIFEST
# ============================================================

extracted_files <- list.files(
  EXTRACT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

if (!length(extracted_files)) {
  stop(
    "No extracted files were found."
  )
}

FILE_MANIFEST <- data.table(
  full_path = normalizePath(
    extracted_files,
    winslash = "/",
    mustWork = TRUE
  )
)

FILE_MANIFEST[
  ,
  file_name :=
    basename(
      full_path
    )
]

FILE_MANIFEST[
  ,
  size_MB :=
    as.numeric(
      file.info(
        full_path
      )$size
    ) / 1024^2
]

FILE_MANIFEST[
  ,
  GSM :=
    fifelse(
      grepl(
        "GSM[0-9]+",
        file_name,
        ignore.case = TRUE
      ),
      sub(
        ".*?(GSM[0-9]+).*",
        "\\1",
        file_name,
        ignore.case = TRUE
      ),
      NA_character_
    )
]

FILE_MANIFEST[
  ,
  donor :=
    fifelse(
      grepl(
        "CF[0-9]+",
        file_name,
        ignore.case = TRUE
      ),
      toupper(
        sub(
          ".*?(CF[0-9]+).*",
          "\\1",
          file_name,
          ignore.case = TRUE
        )
      ),
      NA_character_
    )
]

FILE_MANIFEST[
  ,
  modality :=
    fifelse(
      grepl(
        "RNA",
        file_name,
        ignore.case = TRUE
      ),
      "RNA",
      fifelse(
        grepl(
          "ATAC",
          file_name,
          ignore.case = TRUE
        ),
        "ATAC",
        "UNKNOWN"
      )
    )
]

FILE_MANIFEST <- merge(
  FILE_MANIFEST,
  EXPECTED_DESIGN[
    ,
    .(
      GSM,
      expected_donor = donor,
      expected_rhythm = rhythm,
      expected_modality = modality
    )
  ],
  by = "GSM",
  all.x = TRUE,
  sort = FALSE
)

FILE_MANIFEST[
  ,
  donor_match :=
    !is.na(expected_donor) &
    donor == expected_donor
]

FILE_MANIFEST[
  ,
  modality_match :=
    !is.na(expected_modality) &
    modality == expected_modality
]

setorder(
  FILE_MANIFEST,
  GSM,
  file_name
)

fwrite(
  FILE_MANIFEST,
  file.path(
    QC_DIR,
    "STEP11A_extracted_file_manifest.csv"
  )
)

# ============================================================
# 7. LIGHTWEIGHT TSV HEADER AUDIT
#    No full matrix loading at this stage.
# ============================================================

read_first_lines <- function(
  path,
  n = 2L
) {

  con <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(con),
    add = TRUE
  )

  readLines(
    con,
    n = n,
    warn = FALSE
  )
}

HEADER_AUDIT_LIST <- vector(
  "list",
  nrow(
    FILE_MANIFEST
  )
)

for (i in seq_len(nrow(FILE_MANIFEST))) {

  f <- FILE_MANIFEST$full_path[i]

  first_lines <- tryCatch(
    read_first_lines(
      f,
      n = 2L
    ),
    error = function(e) {
      character(0)
    }
  )

  header_fields <- if (length(first_lines) >= 1L) {
    length(
      strsplit(
        first_lines[1],
        "\t",
        fixed = TRUE
      )[[1]]
    )
  } else {
    NA_integer_
  }

  row1_fields <- if (length(first_lines) >= 2L) {
    length(
      strsplit(
        first_lines[2],
        "\t",
        fixed = TRUE
      )[[1]]
    )
  } else {
    NA_integer_
  }

  header_preview <- if (length(first_lines) >= 1L) {
    paste0(
      substr(
        first_lines[1],
        1,
        300
      )
    )
  } else {
    NA_character_
  }

  HEADER_AUDIT_LIST[[i]] <- data.table(
    GSM = FILE_MANIFEST$GSM[i],
    donor = FILE_MANIFEST$donor[i],
    modality = FILE_MANIFEST$modality[i],
    file_name = FILE_MANIFEST$file_name[i],
    size_MB = FILE_MANIFEST$size_MB[i],
    header_n_fields = header_fields,
    first_data_n_fields = row1_fields,
    header_preview = header_preview
  )
}

HEADER_AUDIT <- rbindlist(
  HEADER_AUDIT_LIST,
  fill = TRUE
)

fwrite(
  HEADER_AUDIT,
  file.path(
    QC_DIR,
    "STEP11A_matrix_header_audit.csv"
  )
)

# ============================================================
# 8. READ AUTHOR CELL METADATA
# ============================================================

meta_con <- gzfile(
  META_GZ,
  open = "rt"
)

META <- fread(
  meta_con,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  showProgress = FALSE
)

close(
  meta_con
)

if (!nrow(META)) {
  stop(
    "Author metadata contains zero rows."
  )
}

# Preserve original metadata exactly as read.
saveRDS(
  META,
  file.path(
    OBJECT_DIR,
    "GSE238242_author_metadata_original.rds"
  )
)

# ============================================================
# 9. METADATA COLUMN AUDIT
# ============================================================

META_COLUMN_AUDIT <- data.table(
  column_index = seq_along(
    names(META)
  ),
  column_name = names(
    META
  ),
  R_class = vapply(
    META,
    function(x) {
      paste(
        class(x),
        collapse = ";"
      )
    },
    character(1)
  ),
  n_unique = vapply(
    META,
    uniqueN,
    integer(1)
  ),
  n_missing = vapply(
    META,
    function(x) {
      sum(
        is.na(x) |
        trimws(
          as.character(x)
        ) == ""
      )
    },
    integer(1)
  )
)

fwrite(
  META_COLUMN_AUDIT,
  file.path(
    QC_DIR,
    "STEP11A_metadata_column_audit.csv"
  )
)

fwrite(
  head(
    META,
    100L
  ),
  file.path(
    QC_DIR,
    "STEP11A_metadata_head100.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 10. LOW-CARDINALITY LEVEL AUDIT
#     Helps identify author cell-type / cluster / donor columns.
# ============================================================

LOW_CARD_LIST <- list()

for (nm in names(META)) {

  n_u <- uniqueN(
    META[[nm]]
  )

  if (
    is.finite(n_u) &&
    n_u >= 1L &&
    n_u <= 50L
  ) {

    vals <- sort(
      unique(
        as.character(
          META[[nm]]
        )
      )
    )

    LOW_CARD_LIST[[
      length(
        LOW_CARD_LIST
      ) + 1L
    ]] <- data.table(
      column = nm,
      n_unique = n_u,
      levels = paste(
        vals,
        collapse = " | "
      )
    )
  }
}

LOW_CARD <- if (length(LOW_CARD_LIST)) {
  rbindlist(
    LOW_CARD_LIST,
    fill = TRUE
  )
} else {
  data.table(
    column = character(),
    n_unique = integer(),
    levels = character()
  )
}

fwrite(
  LOW_CARD,
  file.path(
    QC_DIR,
    "STEP11A_metadata_low_cardinality_levels.csv"
  )
)

# ============================================================
# 11. AUTO-DETECT IMPORTANT METADATA COLUMNS
# ============================================================

pick_col <- function(
  patterns,
  nms
) {

  for (pat in patterns) {

    hit <- grep(
      pat,
      nms,
      value = TRUE,
      ignore.case = TRUE
    )

    if (length(hit)) {
      return(
        hit[1]
      )
    }
  }

  NA_character_
}

COL_BARCODE <- pick_col(
  c(
    "^barcode$",
    "barcode",
    "^cell$",
    "cell.*id"
  ),
  names(META)
)

COL_DONOR <- pick_col(
  c(
    "^donor$",
    "donor",
    "patient",
    "subject",
    "^sample$",
    "sample.*id",
    "orig.ident"
  ),
  names(META)
)

COL_RHYTHM <- pick_col(
  c(
    "^rhythm$",
    "rhythm",
    "^condition$",
    "disease",
    "diagnosis",
    "^group$",
    "status"
  ),
  names(META)
)

COL_CELLTYPE <- pick_col(
  c(
    "^celltype$",
    "^cell_type$",
    "cell.*type",
    "annotation",
    "cell.*label",
    "cluster.*name",
    "celltype"
  ),
  names(META)
)

AUTO_COLS <- data.table(
  role = c(
    "barcode",
    "donor_or_sample",
    "rhythm_or_condition",
    "celltype_or_annotation"
  ),
  detected_column = c(
    COL_BARCODE,
    COL_DONOR,
    COL_RHYTHM,
    COL_CELLTYPE
  )
)

fwrite(
  AUTO_COLS,
  file.path(
    QC_DIR,
    "STEP11A_metadata_auto_detected_columns.csv"
  )
)

# ============================================================
# 12. COUNTS FROM DETECTED AUTHOR LABELS
# ============================================================

write_count_table <- function(
  colname,
  outfile,
  output_name
) {

  if (
    is.na(colname) ||
    !colname %in% names(META)
  ) {
    return(
      invisible(FALSE)
    )
  }

  x <- META[
    ,
    .N,
    by = colname
  ]

  setnames(
    x,
    colname,
    output_name
  )

  setorder(
    x,
    -N
  )

  fwrite(
    x,
    outfile
  )

  invisible(
    TRUE
  )
}

write_count_table(
  COL_DONOR,
  file.path(
    QC_DIR,
    "STEP11A_donor_counts.csv"
  ),
  "donor_or_sample"
)

write_count_table(
  COL_RHYTHM,
  file.path(
    QC_DIR,
    "STEP11A_rhythm_counts.csv"
  ),
  "rhythm_or_condition"
)

write_count_table(
  COL_CELLTYPE,
  file.path(
    QC_DIR,
    "STEP11A_author_celltype_counts.csv"
  ),
  "author_celltype"
)

if (
  !is.na(COL_DONOR) &&
  !is.na(COL_CELLTYPE) &&
  COL_DONOR %in% names(META) &&
  COL_CELLTYPE %in% names(META)
) {

  DONOR_CELLTYPE <- META[
    ,
    .N,
    by = c(
      COL_DONOR,
      COL_CELLTYPE
    )
  ]

  setnames(
    DONOR_CELLTYPE,
    c(
      COL_DONOR,
      COL_CELLTYPE
    ),
    c(
      "donor_or_sample",
      "author_celltype"
    )
  )

  fwrite(
    DONOR_CELLTYPE,
    file.path(
      QC_DIR,
      "STEP11A_donor_by_author_celltype.csv"
    )
  )
}

# ============================================================
# 13. EXPECTED 14-GSM COMPLETENESS
# ============================================================

GSM_CHECK <- EXPECTED_DESIGN[
  ,
  .(
    GSM,
    donor,
    rhythm,
    modality
  )
]

GSM_CHECK[
  ,
  n_extracted_files :=
    vapply(
      GSM,
      function(g) {
        sum(
          FILE_MANIFEST$GSM == g,
          na.rm = TRUE
        )
      },
      integer(1)
    )
]

GSM_CHECK[
  ,
  present :=
    n_extracted_files >= 1L
]

fwrite(
  GSM_CHECK,
  file.path(
    QC_DIR,
    "STEP11A_expected_14GSM_completeness.csv"
  )
)

# ============================================================
# 14. FINAL STEP11A READINESS
# ============================================================

READINESS <- data.table(
  check = c(
    "local_RAW_tar_exists",
    "local_metadata_exists",
    "tar_has_members",
    "extracted_files_exist",
    "all_14_expected_GSM_present",
    "7_unique_expected_donors",
    "4_SR_plus_3_AF_design",
    "author_metadata_readable",
    "author_metadata_nonempty",
    "author_celltype_column_auto_detected"
  ),
  pass = c(
    file.exists(RAW_TAR),
    file.exists(META_GZ),
    length(tar_members) > 0L,
    nrow(FILE_MANIFEST) > 0L,
    all(GSM_CHECK$present),
    uniqueN(EXPECTED_DESIGN$donor) == 7L,
    (
      uniqueN(
        EXPECTED_DESIGN[
          rhythm == "SR",
          donor
        ]
      ) == 4L &&
      uniqueN(
        EXPECTED_DESIGN[
          rhythm == "AF",
          donor
        ]
      ) == 3L
    ),
    TRUE,
    nrow(META) > 0L,
    !is.na(COL_CELLTYPE)
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11A_readiness.csv"
  )
)

# ============================================================
# 15. SESSION INFO + README
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11A_sessionInfo.txt"
  )
)

README <- c(
  "STEP11A — GSE238242 LAA sn-multiome local audit",
  "",
  "Source files were already present under D:/A/data.",
  "No data were downloaded by this script.",
  "",
  "Frozen expected GEO design:",
  "SR donors: CF69, CF77, CF89, CF91",
  "AF donors: CF93, CF97, CF102",
  "7 paired RNA + 7 paired ATAC libraries",
  "Tissue: left atrial appendage (LAA)",
  "Source genome build: GRCh38/hg38",
  "",
  "Important:",
  "STEP11A performs no biological filtering, re-clustering, or re-annotation.",
  "Author-provided cell labels are audited and preserved.",
  "",
  "Send back these files first:",
  "00_QC/STEP11A_readiness.csv",
  "00_QC/STEP11A_extracted_file_manifest.csv",
  "00_QC/STEP11A_matrix_header_audit.csv",
  "00_QC/STEP11A_metadata_column_audit.csv",
  "00_QC/STEP11A_metadata_low_cardinality_levels.csv",
  "00_QC/STEP11A_metadata_auto_detected_columns.csv",
  "00_QC/STEP11A_expected_14GSM_completeness.csv",
  "00_QC/STEP11A_author_celltype_counts.csv (if generated)",
  "00_QC/STEP11A_donor_by_author_celltype.csv (if generated)"
)

writeLines(
  README,
  file.path(
    OUT_ROOT,
    "README_STEP11A.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP11A COMPLETE\n",
  "Output root:\n",
  OUT_ROOT,
  "\n\nReadiness:\n",
  sep = ""
)

print(
  READINESS
)

cat(
  "\nPlease send back the QC files listed in README_STEP11A.txt.\n",
  "Do not start STEP11B until the file structure and author metadata are checked.\n",
  "============================================================\n",
  sep = ""
)
