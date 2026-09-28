# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP10B_V3_2_FINAL_LOCUS_CENTRIC_INTEGRATION.R
# See repository README.md for execution order and external dependencies.


# ============================================================
# STEP10B_V3_2_FINAL_LOCUS_CENTRIC
# Final locus-centric multi-layer evidence integration
# AF / HFpEF / BMI / OSA
#
# SCIENTIFIC PURPOSE
# ------------------
# Fix the central problem of gene-centric integration:
# conjFDR, LAVA, coloc and SuSiE are LOCUS / TRAIT-PAIR evidence,
# not gene-level evidence. Therefore this script:
#
#   1) defines loci from FINAL FUMA GenomicRiskLoci.txt files;
#   2) validates each FUMA locus against FINAL conjFDR loci;
#   3) maps LAVA by genomic interval + trait pair;
#   4) maps coloc by the four pre-specified target regions;
#   5) uses STEP8C_V3 as the FINAL formal SuSiE result for all
#      eight targeted tests;
#   6) adds STEP8C_V4 only as reliability / SER audit annotation;
#   7) maps FUMA genes to their actual loci;
#   8) maps FINAL eQTL-SMR and FULL-GENOME sQTL-SMR results
#      to genes inside their actual loci;
#   9) distinguishes shared regulatory evidence across BOTH
#      traits from one-trait-only regulatory evidence;
#  10) creates journal-ready standalone figures and source data.
#
# IMPORTANT
# ---------
# This script does NOT rerun any upstream analysis.
# It integrates only explicitly selected FINAL files.
#
# NO arbitrary "causal score" is used.
# Candidate ordering is lexicographic and evidence-structure based.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "ggplot2",
  "stringr",
  "tidyr",
  "dplyr"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(data.table)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)

# ============================================================
# 1. FROZEN PATHS
# ============================================================

ROOT <- "D:/A/data/STEP10A/STEP10A_INPUT_FINAL"

CONJ_DIR  <- file.path(ROOT, "01_conjFDR")
FUMA_DIR  <- file.path(ROOT, "02_FUMA")
LAVA_DIR  <- file.path(ROOT, "03_LAVA")
COLOC_DIR <- file.path(ROOT, "04_coloc")
SUSIE_DIR <- file.path(ROOT, "05_SuSiE")
EQTL_DIR  <- file.path(ROOT, "06_eQTL_SMR")
SQTL_DIR  <- file.path(ROOT, "07_sQTL_SMR")

OUT_ROOT <- "D:/A/data/STEP10B/STEP10B_V3_FINAL_LOCUS_CENTRIC_RESULTS"

OUT_DIRS <- c(
  "00_QC",
  "01_SOURCE_DATA",
  "02_MAIN_TABLES",
  "03_SUPPLEMENTARY_TABLES",
  "04_MAIN_FIGURES",
  "05_SUPPLEMENTARY_FIGURES",
  "06_LOGS"
)

for (d in OUT_DIRS) {
  dir.create(
    file.path(OUT_ROOT, d),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 2. FROZEN TRAITS / PAIRS / TARGET REGIONS
# ============================================================

TRAITS <- c(
  "AF",
  "HFpEF",
  "BMI",
  "OSA"
)

EXPECTED_PAIRS <- c(
  "AF__BMI",
  "AF__HFpEF",
  "AF__OSA",
  "HFpEF__BMI",
  "HFpEF__OSA",
  "BMI__OSA"
)

EXPECTED_FUMA_FOLDERS <- c(
  "STD_AF_BMI_conjFDR",
  "STD_AF_HFpEF_conjFDR",
  "STD_AF_OSA_conjFDR",
  "STD_BMI_OSA_conjFDR",
  "STD_HFpEF_BMI_conjFDR",
  "STD_HFpEF_OSA_conjFDR"
)

# Frozen GRCh37 / hg19 regions from STEP8A.
TARGET_REGIONS <- data.table(
  region = c(
    "chr11_recurrent",
    "chr16_locus2135_FTO",
    "chr3_BMI_centered",
    "chr2_AF_HFpEF"
  ),
  chr = c(
    11L,
    16L,
    3L,
    2L
  ),
  start = c(
    16024798L,
    53393883L,
    185170000L,
    200401427L
  ),
  stop = c(
    17024798L,
    54866095L,
    186170000L,
    201401427L
  )
)

# The eight formal V3 targeted SuSiE tests.
FORMAL_TESTS <- data.table(
  region = c(
    "chr16_locus2135_FTO",
    "chr16_locus2135_FTO",
    "chr16_locus2135_FTO",
    "chr3_BMI_centered",
    "chr3_BMI_centered",
    "chr11_recurrent",
    "chr11_recurrent",
    "chr2_AF_HFpEF"
  ),
  trait1 = c(
    "AF",
    "BMI",
    "AF",
    "HFpEF",
    "AF",
    "AF",
    "AF",
    "AF"
  ),
  trait2 = c(
    "BMI",
    "OSA",
    "OSA",
    "BMI",
    "BMI",
    "OSA",
    "HFpEF",
    "HFpEF"
  )
)

# ============================================================
# 3. HELPERS
# ============================================================

canon_pair <- function(a, b) {

  a <- as.character(a)
  b <- as.character(b)

  ia <- match(a, TRAITS)
  ib <- match(b, TRAITS)

  if (is.na(ia) || is.na(ib)) {
    return(
      paste(
        sort(c(a, b)),
        collapse = "__"
      )
    )
  }

  if (ia < ib) {
    paste(a, b, sep = "__")
  } else {
    paste(b, a, sep = "__")
  }
}

split_pair <- function(pair_key) {

  z <- strsplit(
    pair_key,
    "__",
    fixed = TRUE
  )[[1]]

  if (length(z) != 2L) {
    return(c(NA_character_, NA_character_))
  }

  z
}

clean_gene <- function(x) {

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  y[
    y %in% c(
      "",
      "NA",
      "N/A",
      "NAN",
      "NULL"
    )
  ] <- NA_character_

  y
}

as_num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

as_int <- function(x) {
  suppressWarnings(
    as.integer(
      as.character(x)
    )
  )
}

as_logical_safe <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  y %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "PASS"
  )
}

collapse_unique <- function(x) {

  z <- sort(
    unique(
      as.character(
        x[
          !is.na(x) &
          nzchar(
            trimws(
              as.character(x)
            )
          )
        ]
      )
    )
  )

  if (!length(z)) {
    return("")
  }

  paste(
    z,
    collapse = ";"
  )
}

safe_min <- function(x) {

  z <- x[
    is.finite(x)
  ]

  if (!length(z)) {
    return(NA_real_)
  }

  min(z)
}

pick_one_exact <- function(
  directory,
  preferred_names,
  label,
  recursive = TRUE,
  required = TRUE
) {

  all_files <- list.files(
    directory,
    recursive = recursive,
    full.names = TRUE
  )

  for (nm in preferred_names) {

    hit <- all_files[
      basename(all_files) == nm
    ]

    if (length(hit) == 1L) {
      return(
        normalizePath(
          hit,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }

    if (length(hit) > 1L) {
      stop(
        "More than one exact file matched for ",
        label,
        ": ",
        nm,
        "\n",
        paste(
          hit,
          collapse = "\n"
        )
      )
    }
  }

  if (required) {
    stop(
      "Required FINAL file not found for ",
      label,
      ". Expected one of:\n",
      paste(
        preferred_names,
        collapse = "\n"
      ),
      "\nUnder:\n",
      directory
    )
  }

  NA_character_
}

pick_col <- function(
  nms,
  candidates,
  required = TRUE,
  label = "column"
) {

  low <- tolower(nms)

  for (cand in candidates) {

    idx <- which(
      low == tolower(cand)
    )

    if (length(idx)) {
      return(
        nms[idx[1]]
      )
    }
  }

  if (required) {
    stop(
      "Cannot identify ",
      label,
      ". Available columns:\n",
      paste(
        nms,
        collapse = ", "
      )
    )
  }

  NA_character_
}

interval_overlap <- function(
  a_start,
  a_stop,
  b_start,
  b_stop
) {

  is.finite(a_start) &
  is.finite(a_stop) &
  is.finite(b_start) &
  is.finite(b_stop) &
  a_start <= b_stop &
  a_stop >= b_start
}

parse_test_id <- function(test_id) {

  z <- strsplit(
    as.character(test_id),
    "__",
    fixed = TRUE
  )[[1]]

  if (length(z) < 3L) {
    return(
      list(
        region = NA_character_,
        trait1 = NA_character_,
        trait2 = NA_character_,
        pair_key = NA_character_
      )
    )
  }

  t1 <- z[length(z) - 1L]
  t2 <- z[length(z)]

  region <- paste(
    z[
      seq_len(
        length(z) - 2L
      )
    ],
    collapse = "__"
  )

  list(
    region = region,
    trait1 = t1,
    trait2 = t2,
    pair_key = canon_pair(t1, t2)
  )
}

extract_folder_pair <- function(folder_name) {

  x <- folder_name

  x <- sub(
    "^STD_",
    "",
    x
  )

  x <- sub(
    "_conjFDR$",
    "",
    x
  )

  parts <- strsplit(
    x,
    "_",
    fixed = TRUE
  )[[1]]

  # HFpEF is one token because folder separator is underscore
  # only between traits.
  possible <- list(
    c("AF", "BMI"),
    c("AF", "HFpEF"),
    c("AF", "OSA"),
    c("BMI", "OSA"),
    c("HFpEF", "BMI"),
    c("HFpEF", "OSA")
  )

  for (p in possible) {

    if (
      identical(
        x,
        paste(
          p,
          collapse = "_"
        )
      )
    ) {
      return(
        canon_pair(
          p[1],
          p[2]
        )
      )
    }
  }

  stop(
    "Cannot parse FUMA pair from folder: ",
    folder_name
  )
}

pair_to_display <- function(pair_key) {

  p <- split_pair(pair_key)

  paste(
    p[1],
    p[2],
    sep = "–"
  )
}

# ============================================================
# 4. STRICT FINAL INPUT DISCOVERY
# ============================================================

for (d in c(
  CONJ_DIR,
  FUMA_DIR,
  LAVA_DIR,
  COLOC_DIR,
  SUSIE_DIR,
  EQTL_DIR,
  SQTL_DIR
)) {

  if (!dir.exists(d)) {
    stop(
      "Required final-input directory missing:\n",
      d
    )
  }
}

# Exact final LAVA.
LAVA_FILE <- pick_one_exact(
  LAVA_DIR,
  c(
    "STEP7D_local_rg_all.csv"
  ),
  "FINAL LAVA local rg"
)

# Prefer corrected coloc if present.
COLOC_FILE <- pick_one_exact(
  COLOC_DIR,
  c(
    "STEP8A_coloc_key_results_CORRECTED.csv",
    "STEP8A_coloc_key_results.csv"
  ),
  "FINAL coloc ABF"
)

# CRITICAL:
# V3 = formal 8-test multi-signal SuSiE analysis.
SUSIE_V3_FILE <- pick_one_exact(
  SUSIE_DIR,
  c(
    "STEP8C_V3_coloc_susie_test_summary.csv"
  ),
  "FINAL formal SuSiE V3 summary"
)

# V4 = reliability/SER audit ONLY.
SUSIE_V4_FILE <- pick_one_exact(
  SUSIE_DIR,
  c(
    "STEP8C_V4_SER_audit_summary.csv"
  ),
  "FINAL SuSiE V4 SER audit"
)

EQTL_FILE <- pick_one_exact(
  EQTL_DIR,
  c(
    "STEP9B_all_32_SMR_HEIDI_results.csv.gz"
  ),
  "FINAL GTEx eQTL-SMR"
)

SQTL_FILE <- pick_one_exact(
  SQTL_DIR,
  c(
    "STEP9C2_V3_all_32_FULLGENOME_sQTL_SMR_HEIDI_results.csv.gz"
  ),
  "FINAL full-genome GTEx sQTL-SMR"
)

# ConjFDR files: exact six final non-zscore locus tables.
CONJ_EXPECTED <- c(
  AF_BMI = "AF_BMI_conjfdr_0.05_loci.csv",
  AF_HFpEF = "AF_HFpEF_conjfdr_0.05_loci.csv",
  AF_OSA = "AF_OSA_conjfdr_0.05_loci.csv",
  BMI_OSA = "BMI_OSA_conjfdr_0.05_loci.csv",
  HFpEF_BMI = "HFpEF_BMI_conjfdr_0.05_loci.csv",
  HFpEF_OSA = "HFpEF_OSA_conjfdr_0.05_loci.csv"
)

CONJ_FILES <- character()

all_conj_files <- list.files(
  CONJ_DIR,
  recursive = TRUE,
  full.names = TRUE
)

for (nm in CONJ_EXPECTED) {

  hit <- all_conj_files[
    basename(all_conj_files) == nm
  ]

  if (length(hit) != 1L) {
    stop(
      "Expected exactly one FINAL conjFDR file:\n",
      nm,
      "\nFound: ",
      length(hit)
    )
  }

  CONJ_FILES <- c(
    CONJ_FILES,
    normalizePath(
      hit,
      winslash = "/",
      mustWork = TRUE
    )
  )
}

names(CONJ_FILES) <- names(CONJ_EXPECTED)

# FUMA folders.
actual_fuma_dirs <- list.dirs(
  FUMA_DIR,
  recursive = FALSE,
  full.names = FALSE
)

missing_fuma <- setdiff(
  EXPECTED_FUMA_FOLDERS,
  actual_fuma_dirs
)

if (length(missing_fuma)) {
  stop(
    "Missing required extracted FINAL FUMA folder(s):\n",
    paste(
      missing_fuma,
      collapse = "\n"
    )
  )
}

# ============================================================
# 5. INPUT MANIFEST
# ============================================================

INPUT_MANIFEST <- rbindlist(
  list(
    data.table(
      module = "LAVA",
      role = "final",
      file = LAVA_FILE
    ),
    data.table(
      module = "coloc_ABF",
      role = "final",
      file = COLOC_FILE
    ),
    data.table(
      module = "SuSiE_V3",
      role = "formal_final",
      file = SUSIE_V3_FILE
    ),
    data.table(
      module = "SuSiE_V4",
      role = "reliability_audit",
      file = SUSIE_V4_FILE
    ),
    data.table(
      module = "eQTL_SMR",
      role = "final",
      file = EQTL_FILE
    ),
    data.table(
      module = "sQTL_SMR",
      role = "final_fullgenome",
      file = SQTL_FILE
    ),
    data.table(
      module = "conjFDR",
      role = names(CONJ_FILES),
      file = unname(CONJ_FILES)
    )
  ),
  fill = TRUE
)

fwrite(
  INPUT_MANIFEST,
  file.path(
    OUT_ROOT,
    "00_QC",
    "STEP10B_V3_final_input_manifest.csv"
  )
)

# ============================================================
# 6. READ FINAL conjFDR LOCI
# ============================================================

CONJ_LIST <- list()

for (nm in names(CONJ_FILES)) {

  f <- CONJ_FILES[[nm]]

  x <- fread(
    f,
    showProgress = FALSE
  )

  pair_parts <- strsplit(
    nm,
    "_",
    fixed = TRUE
  )[[1]]

  if (length(pair_parts) != 2L) {
    stop(
      "Unexpected conjFDR pair key: ",
      nm
    )
  }

  pair_key <- canon_pair(
    pair_parts[1],
    pair_parts[2]
  )

  conj_col <- grep(
    "^conjfdr_",
    names(x),
    value = TRUE,
    ignore.case = TRUE
  )

  if (!length(conj_col)) {

    if (
      "min_conjfdr" %in%
      names(x)
    ) {
      conj_col <- "min_conjfdr"
    } else {
      stop(
        "Cannot identify conjFDR value column in ",
        f
      )
    }
  }

  conj_col <- conj_col[1]

  z <- data.table(
    pair_key = pair_key,
    pair_display = pair_to_display(
      pair_key
    ),
    locusnum_conj = as_int(
      x$locusnum
    ),
    lead_snp_conj = as.character(
      x$snpid
    ),
    chr = as_int(
      x$chrnum
    ),
    pos = as_num(
      x$chrpos
    ),
    conjFDR = as_num(
      x[[conj_col]]
    ),
    source_file = normalizePath(
      f,
      winslash = "/",
      mustWork = TRUE
    )
  )

  CONJ_LIST[[
    length(CONJ_LIST) + 1L
  ]] <- z
}

CONJ <- rbindlist(
  CONJ_LIST,
  fill = TRUE
)

stopifnot(
  setequal(
    unique(CONJ$pair_key),
    EXPECTED_PAIRS
  )
)

fwrite(
  CONJ,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_conjFDR_final_loci.csv"
  )
)

# ============================================================
# 7. READ FINAL FUMA LOCI + GENES
# ============================================================

FUMA_LOCUS_LIST <- list()
FUMA_GENE_LIST <- list()
FUMA_LEAD_LIST <- list()
FUMA_MULTI_LOCUS_QC_LIST <- list()

for (folder in EXPECTED_FUMA_FOLDERS) {

  fd <- file.path(
    FUMA_DIR,
    folder
  )

  pair_key <- extract_folder_pair(
    folder
  )

  pair_traits <- split_pair(
    pair_key
  )

  loci_file <- file.path(
    fd,
    "GenomicRiskLoci.txt"
  )

  genes_file <- file.path(
    fd,
    "genes.txt"
  )

  lead_file <- file.path(
    fd,
    "leadSNPs.txt"
  )

  required_fuma <- c(
    loci_file,
    genes_file,
    lead_file
  )

  if (any(!file.exists(required_fuma))) {
    stop(
      "Required FUMA file missing in ",
      fd,
      ":\n",
      paste(
        required_fuma[
          !file.exists(required_fuma)
        ],
        collapse = "\n"
      )
    )
  }

  loci <- fread(
    loci_file,
    showProgress = FALSE
  )

  genes <- fread(
    genes_file,
    showProgress = FALSE
  )

  leads <- fread(
    lead_file,
    showProgress = FALSE
  )

  # ----------------------------------------------------------
  # V3.1 FIX — STANDARDIZE FUMA JOIN KEYS
  # ----------------------------------------------------------
  # FUMA can write GenomicLocus as an integer in
  # GenomicRiskLoci.txt but as character in genes.txt.
  # data.table::merge() requires identical join-key types.
  #
  # We therefore force the locus identifier to CHARACTER in
  # every FUMA table before any merge. This changes only the
  # storage type of the identifier, not its value.
  # ----------------------------------------------------------

  if (!"GenomicLocus" %in% names(loci)) {
    stop(
      "GenomicRiskLoci.txt lacks GenomicLocus in: ",
      fd
    )
  }

  if (!"GenomicLocus" %in% names(genes)) {
    stop(
      "genes.txt lacks GenomicLocus in: ",
      fd
    )
  }

  loci[
    ,
    GenomicLocus :=
      trimws(
        as.character(
          GenomicLocus
        )
      )
  ]

  # ----------------------------------------------------------
  # V3.2 FIX — EXPAND LEGAL MULTI-LOCUS FUMA GENE MAPPINGS
  # ----------------------------------------------------------
  # According to FUMA's SNP2GENE output specification,
  # genes.txt -> GenomicLocus can contain MULTIPLE locus indices
  # separated by ":" when mapped SNPs come from distinct genomic
  # risk loci. Example: "1:2".
  #
  # GenomicRiskLoci.txt, however, contains one locus index per row.
  # Therefore exact matching BEFORE expansion is scientifically
  # incorrect and falsely labels valid multi-locus genes as
  # unmatched.
  #
  # We preserve the original field and expand:
  #   "1:2" -> one copy linked to locus 1
  #          -> one copy linked to locus 2
  #
  # This is a one-to-many normalization only. No gene, SNP,
  # association statistic, or locus is invented or discarded.
  # ----------------------------------------------------------

  genes[
    ,
    GenomicLocus_raw :=
      trimws(
        as.character(
          GenomicLocus
        )
      )
  ]

  genes[
    ,
    FUMA_multi_locus_gene :=
      !is.na(GenomicLocus_raw) &
      grepl(
        ":",
        GenomicLocus_raw,
        fixed = TRUE
      )
  ]

  genes[
    ,
    FUMA_gene_row_id__ :=
      seq_len(.N)
  ]

  locus_split__ <- strsplit(
    genes$GenomicLocus_raw,
    ":",
    fixed = TRUE
  )

  locus_split__ <- lapply(
    locus_split__,
    function(v) {

      if (!length(v)) {
        return(NA_character_)
      }

      v <- trimws(
        as.character(v)
      )

      v[
        v %in% c(
          "",
          "NA",
          "N/A",
          "NaN",
          "NULL"
        )
      ] <- NA_character_

      v <- v[
        !is.na(v)
      ]

      if (!length(v)) {
        return(NA_character_)
      }

      unique(v)
    }
  )

  expansion_n__ <- lengths(
    locus_split__
  )

  expansion_n__[
    expansion_n__ < 1L
  ] <- 1L

  genes_expanded__ <- genes[
    rep(
      seq_len(.N),
      times = expansion_n__
    )
  ]

  genes_expanded__[
    ,
    GenomicLocus :=
      unlist(
        locus_split__,
        use.names = FALSE
      )
  ]

  # Defensive normalization after expansion.
  genes_expanded__[
    ,
    GenomicLocus :=
      trimws(
        as.character(
          GenomicLocus
        )
      )
  ]

  genes_expanded__[
    GenomicLocus %in% c(
      "",
      "NA",
      "N/A",
      "NaN",
      "NULL"
    ),
    GenomicLocus := NA_character_
  ]

  FUMA_MULTI_LOCUS_QC_LIST[[
    length(FUMA_MULTI_LOCUS_QC_LIST) + 1L
  ]] <- data.table(
    source_folder = folder,
    original_gene_rows = nrow(genes),
    multi_locus_gene_rows =
      sum(
        genes$FUMA_multi_locus_gene,
        na.rm = TRUE
      ),
    expanded_gene_locus_rows =
      nrow(
        genes_expanded__
      ),
    extra_rows_created_by_expansion =
      nrow(
        genes_expanded__
      ) -
      nrow(
        genes
      ),
    max_loci_per_gene_row =
      max(
        expansion_n__,
        na.rm = TRUE
      )
  )

  genes <- genes_expanded__

  rm(
    genes_expanded__,
    locus_split__,
    expansion_n__
  )

  if ("GenomicLocus" %in% names(leads)) {
    leads[
      ,
      GenomicLocus :=
        trimws(
          as.character(
            GenomicLocus
          )
        )
    ]
  }

  loci[
    ,
    `:=`(
      source_folder = folder,
      pair_key = pair_key,
      pair_display = pair_to_display(
        pair_key
      ),
      trait1 = pair_traits[1],
      trait2 = pair_traits[2]
    )
  ]

  loci[
    ,
    locus_uid := paste(
      pair_key,
      paste0(
        "chr",
        chr,
        ":",
        start,
        "-",
        end
      ),
      paste0(
        "FUMA",
        GenomicLocus
      ),
      sep = "|"
    )
  ]

  # Standardized genes.
  genes[
    ,
    `:=`(
      source_folder = folder,
      pair_key = pair_key,
      pair_display = pair_to_display(
        pair_key
      ),
      trait1 = pair_traits[1],
      trait2 = pair_traits[2],
      gene = clean_gene(
        symbol
      )
    )
  ]

  genes[
    ,
    locus_uid := paste(
      pair_key,
      paste0(
        "FUMA",
        GenomicLocus
      ),
      sep = "|"
    )
  ]

  # Temporary key; replaced with exact locus_uid after merge.
  locus_key <- loci[
    ,
    .(
      pair_key,
      GenomicLocus,
      locus_uid_full = locus_uid,
      locus_chr = as_int(chr),
      locus_start = as_num(start),
      locus_stop = as_num(end),
      locus_lead = as.character(rsID)
    )
  ]

  # Hard type QC AFTER multi-locus expansion and before the locus-gene join.
  if (
    !identical(
      typeof(genes$GenomicLocus),
      typeof(locus_key$GenomicLocus)
    )
  ) {
    stop(
      "Internal FUMA join-key type mismatch remains in ",
      folder,
      ": genes=",
      typeof(genes$GenomicLocus),
      ", locus_key=",
      typeof(locus_key$GenomicLocus)
    )
  }

  genes <- merge(
    genes,
    locus_key,
    by = c(
      "pair_key",
      "GenomicLocus"
    ),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

  # After correct ":" expansion, every NON-MISSING locus token
  # should map to GenomicRiskLoci.txt. Only a true residual mismatch
  # is considered an error.
  bad_gene_locus <- genes[
    is.na(locus_uid_full) &
    !is.na(GenomicLocus) &
    nzchar(GenomicLocus)
  ]

  if (nrow(bad_gene_locus)) {

    bad_file <- file.path(
      OUT_ROOT,
      "00_QC",
      paste0(
        "FUMA_gene_locus_TRUE_unmatched_after_colon_expansion__",
        folder,
        ".csv"
      )
    )

    fwrite(
      bad_gene_locus,
      bad_file
    )

    stop(
      "After correct FUMA ':' multi-locus expansion, ",
      nrow(bad_gene_locus),
      " locus token(s) still cannot be matched to ",
      "GenomicRiskLoci.txt in ",
      folder,
      ". This is now a TRUE residual mismatch. Inspect: ",
      bad_file
    )
  }

  genes[
    ,
    locus_uid := locus_uid_full
  ]

  genes[
    ,
    locus_uid_full := NULL
  ]

  genes[
    ,
    FUMA_gene_row_id__ := NULL
  ]

  leads[
    ,
    `:=`(
      source_folder = folder,
      pair_key = pair_key,
      pair_display = pair_to_display(
        pair_key
      )
    )
  ]

  FUMA_LOCUS_LIST[[
    length(FUMA_LOCUS_LIST) + 1L
  ]] <- loci

  FUMA_GENE_LIST[[
    length(FUMA_GENE_LIST) + 1L
  ]] <- genes

  FUMA_LEAD_LIST[[
    length(FUMA_LEAD_LIST) + 1L
  ]] <- leads
}

FUMA_LOCI <- rbindlist(
  FUMA_LOCUS_LIST,
  fill = TRUE
)

FUMA_GENES <- rbindlist(
  FUMA_GENE_LIST,
  fill = TRUE
)

FUMA_LEADS <- rbindlist(
  FUMA_LEAD_LIST,
  fill = TRUE
)

FUMA_MULTI_LOCUS_QC <- rbindlist(
  FUMA_MULTI_LOCUS_QC_LIST,
  fill = TRUE
)

FUMA_MULTI_LOCUS_QC_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V3_2_FUMA_multi_locus_gene_expansion_QC.csv"
)

fwrite(
  FUMA_MULTI_LOCUS_QC,
  FUMA_MULTI_LOCUS_QC_FILE
)

FUMA_GENES <- FUMA_GENES[
  !is.na(gene)
]

# ============================================================
# 8. VALIDATE FUMA LOCI AGAINST FINAL conjFDR LOCI
# ============================================================

CONJ_KEY <- CONJ[
  ,
  .(
    pair_key,
    lead_snp_conj,
    chr_conj = chr,
    pos_conj = pos,
    conjFDR_final = conjFDR
  )
]

FUMA_LOCI <- merge(
  FUMA_LOCI,
  CONJ_KEY,
  by.x = c(
    "pair_key",
    "rsID"
  ),
  by.y = c(
    "pair_key",
    "lead_snp_conj"
  ),
  all.x = TRUE
)

FUMA_LOCI[
  ,
  conjFDR_match :=
    !is.na(
      conjFDR_final
    )
]

FUMA_LOCI[
  ,
  conjFDR_value :=
    fifelse(
      !is.na(
        conjFDR_final
      ),
      conjFDR_final,
      as_num(p)
    )
]

FUMA_LOCI[
  ,
  conjFDR_value_source :=
    fifelse(
      conjFDR_match,
      "final_conjFDR_loci_file",
      "FUMA_GenomicRiskLoci_p"
    )
]

# Hard QC: FUMA p is itself based on the submitted conjFDR.
# Exact lead-SNP mismatch is allowed only as a flagged QC item
# because FUMA can merge nearby submitted loci.
CONJ_VALIDATION <- FUMA_LOCI[
  ,
  .(
    n_FUMA_loci = .N,
    n_exact_lead_match =
      sum(
        conjFDR_match
      ),
    fraction_exact_lead_match =
      mean(
        conjFDR_match
      )
  ),
  by = pair_key
]

fwrite(
  CONJ_VALIDATION,
  file.path(
    OUT_ROOT,
    "00_QC",
    "STEP10B_V3_FUMA_vs_conjFDR_validation.csv"
  )
)

# ============================================================
# 9. READ FINAL LAVA
# ============================================================

LAVA <- fread(
  LAVA_FILE,
  showProgress = FALSE
)

required_lava <- c(
  "phen1",
  "phen2",
  "rho",
  "p",
  "locus",
  "chr",
  "start",
  "stop"
)

if (!all(
  required_lava %in%
  names(LAVA)
)) {
  stop(
    "FINAL LAVA table lacks required columns."
  )
}

LAVA[
  ,
  pair_key_std :=
    mapply(
      canon_pair,
      phen1,
      phen2
    )
]

if (
  !"pass_pair_bonferroni" %in%
  names(LAVA)
) {
  stop(
    "FINAL LAVA table lacks pass_pair_bonferroni."
  )
}

LAVA[
  ,
  pass_pair_bonferroni_std :=
    as_logical_safe(
      pass_pair_bonferroni
    )
]

if (
  "pass_BH_all" %in%
  names(LAVA)
) {
  LAVA[
    ,
    pass_BH_all_std :=
      as_logical_safe(
        pass_BH_all
      )
  ]
} else {
  LAVA[
    ,
    pass_BH_all_std := FALSE
  ]
}

fwrite(
  LAVA,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_LAVA_final.csv"
  )
)

# ============================================================
# 10. MAP LAVA TO EACH FUMA LOCUS BY PAIR + INTERVAL
# ============================================================

LAVA_MAP_LIST <- list()

for (i in seq_len(nrow(FUMA_LOCI))) {

  rr <- FUMA_LOCI[i]

  z <- LAVA[
    pair_key_std ==
      rr$pair_key &
    as_int(chr) ==
      as_int(rr$chr) &
    interval_overlap(
      as_num(start),
      as_num(stop),
      as_num(rr$start),
      as_num(rr$end)
    )
  ]

  if (!nrow(z)) {

    LAVA_MAP_LIST[[
      length(LAVA_MAP_LIST) + 1L
    ]] <- data.table(
      locus_uid = rr$locus_uid,
      LAVA_tested = FALSE,
      LAVA_n_overlap_blocks = 0L,
      LAVA_pair_bonf_sig = FALSE,
      LAVA_BH_sig = FALSE,
      LAVA_min_p = NA_real_,
      LAVA_rho_at_min_p = NA_real_,
      LAVA_direction_at_min_p = NA_character_,
      LAVA_overlap_loci = ""
    )

    next
  }

  best <- which.min(
    z$p
  )

  LAVA_MAP_LIST[[
    length(LAVA_MAP_LIST) + 1L
  ]] <- data.table(
    locus_uid = rr$locus_uid,
    LAVA_tested = TRUE,
    LAVA_n_overlap_blocks = nrow(z),
    LAVA_pair_bonf_sig =
      any(
        z$pass_pair_bonferroni_std,
        na.rm = TRUE
      ),
    LAVA_BH_sig =
      any(
        z$pass_BH_all_std,
        na.rm = TRUE
      ),
    LAVA_min_p =
      as_num(
        z$p[best]
      ),
    LAVA_rho_at_min_p =
      as_num(
        z$rho[best]
      ),
    LAVA_direction_at_min_p =
      if (
        "direction" %in%
        names(z)
      ) {
        as.character(
          z$direction[best]
        )
      } else {
        ifelse(
          z$rho[best] >= 0,
          "positive",
          "negative"
        )
      },
    LAVA_overlap_loci =
      collapse_unique(
        z$locus
      )
  )
}

LAVA_MAP <- rbindlist(
  LAVA_MAP_LIST,
  fill = TRUE
)

# ============================================================
# 11. READ FINAL coloc ABF
# ============================================================

COLOC <- fread(
  COLOC_FILE,
  showProgress = FALSE
)

if (!"test_id" %in% names(COLOC)) {
  stop(
    "FINAL coloc file lacks test_id."
  )
}

COLOC_PARSE <- lapply(
  COLOC$test_id,
  parse_test_id
)

COLOC[
  ,
  region :=
    vapply(
      COLOC_PARSE,
      function(x) x$region,
      character(1)
    )
]

COLOC[
  ,
  trait1 :=
    vapply(
      COLOC_PARSE,
      function(x) x$trait1,
      character(1)
    )
]

COLOC[
  ,
  trait2 :=
    vapply(
      COLOC_PARSE,
      function(x) x$trait2,
      character(1)
    )
]

COLOC[
  ,
  pair_key_std :=
    vapply(
      COLOC_PARSE,
      function(x) x$pair_key,
      character(1)
    )
]

COLOC[
  ,
  ABF_class :=
    fifelse(
      as_num(`PP.H4`) >= 0.80 &
      as_num(`PP.H4`) >
        as_num(`PP.H3`),
      "H4_strong",

      fifelse(
        as_num(`PP.H3`) >= 0.80 &
        as_num(`PP.H3`) >
          as_num(`PP.H4`),
        "H3_strong",
        "non_dominant"
      )
    )
]

fwrite(
  COLOC,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_coloc_ABF_final.csv"
  )
)

# ============================================================
# 12. READ FINAL FORMAL SuSiE V3
# ============================================================

SUSIE_V3 <- fread(
  SUSIE_V3_FILE,
  showProgress = FALSE
)

required_v3 <- c(
  "test_id",
  "region",
  "trait1",
  "trait2",
  "n_signal_pairs_default",
  "max_PP_H4_default",
  "max_PP_H3_default",
  "n_signal_pairs_conservative",
  "max_PP_H4_conservative",
  "max_PP_H3_conservative"
)

if (!all(
  required_v3 %in%
  names(SUSIE_V3)
)) {
  stop(
    "Formal STEP8C_V3 summary lacks required columns."
  )
}

SUSIE_V3[
  ,
  pair_key_std :=
    mapply(
      canon_pair,
      trait1,
      trait2
    )
]

SUSIE_V3[
  ,
  SuSiE_class :=
    fifelse(
      as_int(
        n_signal_pairs_default
      ) == 0L,
      "No_signal_pair",

      fifelse(
        as_num(
          max_PP_H3_default
        ) >= 0.80 &
        as_num(
          max_PP_H4_default
        ) < 0.20,
        "H3_dominant",

        fifelse(
          as_num(
            max_PP_H4_default
          ) >= 0.80 &
          as_num(
            max_PP_H4_conservative
          ) >= 0.80,
          "H4_robust_both_priors",

          fifelse(
            as_num(
              max_PP_H4_default
            ) >= 0.80 &
            (
              is.na(
                as_num(
                  max_PP_H4_conservative
                )
              ) |
              as_num(
                max_PP_H4_conservative
              ) < 0.80
            ),
            "H4_prior_sensitive",
            "No_dominant_H3_H4"
          )
        )
      )
    )
]

fwrite(
  SUSIE_V3,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_SuSiE_V3_formal_final.csv"
  )
)

# ============================================================
# 13. READ V4 SER AUDIT AS ANNOTATION ONLY
# ============================================================

SUSIE_V4 <- fread(
  SUSIE_V4_FILE,
  showProgress = FALSE
)

if (!"test_id" %in% names(SUSIE_V4)) {
  stop(
    "STEP8C_V4 SER audit lacks test_id."
  )
}

SUSIE_V4_AUDIT <- SUSIE_V4[
  ,
  .(
    test_id,
    V4_audited = TRUE,
    V4_trait = trait,
    V4_partner_trait = partner_trait,
    V4_R_reliability_flag = V3_R_reliability_flag,
    V4_SER_top_snp = SER_top_snp,
    V4_SER_top_PIP = as_num(SER_top_PIP),
    V4_SER_CS_size = as_num(SER_CS_size),
    V4_SER_partner_max_Jaccard =
      as_num(
        SER_partner_max_Jaccard
      ),
    V4_SER_top_in_partner_CS =
      as_logical_safe(
        SER_top_in_partner_any_EB_CS
      ),
    V4_audit_label =
      as.character(
        audit_label
      ),
    V4_partner_cluster_label =
      as.character(
        partner_cluster_label
      )
  )
]

SUSIE_V3 <- merge(
  SUSIE_V3,
  SUSIE_V4_AUDIT,
  by = "test_id",
  all.x = TRUE
)

SUSIE_V3[
  is.na(V4_audited),
  V4_audited := FALSE
]

fwrite(
  SUSIE_V4,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_SuSiE_V4_SER_audit.csv"
  )
)

# ============================================================
# 14. MAP FUMA LOCI TO THE FOUR FROZEN TARGET REGIONS
# ============================================================

TARGET_MAP_LIST <- list()

for (i in seq_len(nrow(FUMA_LOCI))) {

  rr <- FUMA_LOCI[i]

  hits <- TARGET_REGIONS[
    as_int(chr) ==
      as_int(rr$chr) &
    interval_overlap(
      as_num(start),
      as_num(stop),
      as_num(rr$start),
      as_num(rr$end)
    )
  ]

  TARGET_MAP_LIST[[
    length(TARGET_MAP_LIST) + 1L
  ]] <- data.table(
    locus_uid = rr$locus_uid,
    target_region_count = nrow(hits),
    target_regions =
      if (
        nrow(hits)
      ) {
        collapse_unique(
          hits$region
        )
      } else {
        ""
      }
  )
}

TARGET_MAP <- rbindlist(
  TARGET_MAP_LIST,
  fill = TRUE
)

# ============================================================
# 15. LOCUS-LEVEL coloc / SuSiE MAPPING
# ============================================================

COLOC_LOCUS_LIST <- list()
SUSIE_LOCUS_LIST <- list()

for (i in seq_len(nrow(FUMA_LOCI))) {

  rr <- FUMA_LOCI[i]

  # Regions overlapping this FUMA locus.
  reg_hits <- TARGET_REGIONS[
    as_int(chr) ==
      as_int(rr$chr) &
    interval_overlap(
      as_num(start),
      as_num(stop),
      as_num(rr$start),
      as_num(rr$end)
    ),
    region
  ]

  # --------------------
  # coloc ABF
  # --------------------
  cz <- COLOC[
    pair_key_std ==
      rr$pair_key &
    region %in%
      reg_hits
  ]

  if (!nrow(cz)) {

    COLOC_LOCUS_LIST[[
      length(COLOC_LOCUS_LIST) + 1L
    ]] <- data.table(
      locus_uid = rr$locus_uid,
      coloc_tested = FALSE,
      coloc_test_ids = "",
      coloc_max_PP_H4 = NA_real_,
      coloc_max_PP_H3 = NA_real_,
      coloc_ABF_class = "Not_tested",
      coloc_interpretation = ""
    )

  } else {

    # choose row with the larger of PP.H3/H4 for summary
    strength <- pmax(
      as_num(cz$`PP.H3`),
      as_num(cz$`PP.H4`),
      na.rm = TRUE
    )

    best <- which.max(
      strength
    )

    COLOC_LOCUS_LIST[[
      length(COLOC_LOCUS_LIST) + 1L
    ]] <- data.table(
      locus_uid = rr$locus_uid,
      coloc_tested = TRUE,
      coloc_test_ids =
        collapse_unique(
          cz$test_id
        ),
      coloc_max_PP_H4 =
        max(
          as_num(cz$`PP.H4`),
          na.rm = TRUE
        ),
      coloc_max_PP_H3 =
        max(
          as_num(cz$`PP.H3`),
          na.rm = TRUE
        ),
      coloc_ABF_class =
        as.character(
          cz$ABF_class[best]
        ),
      coloc_interpretation =
        if (
          "interpretation_revised" %in%
          names(cz)
        ) {
          as.character(
            cz$interpretation_revised[best]
          )
        } else {
          ""
        }
    )
  }

  # --------------------
  # formal SuSiE V3
  # --------------------
  sz <- SUSIE_V3[
    pair_key_std ==
      rr$pair_key &
    region %in%
      reg_hits
  ]

  if (!nrow(sz)) {

    SUSIE_LOCUS_LIST[[
      length(SUSIE_LOCUS_LIST) + 1L
    ]] <- data.table(
      locus_uid = rr$locus_uid,
      SuSiE_tested = FALSE,
      SuSiE_test_ids = "",
      SuSiE_class = "Not_tested",
      SuSiE_max_H4_default = NA_real_,
      SuSiE_max_H4_conservative = NA_real_,
      SuSiE_max_H3_default = NA_real_,
      V4_audited = FALSE,
      V4_audit_labels = ""
    )

  } else {

    # At most one V3 test per region/pair in the frozen design.
    ss <- sz[1]

    SUSIE_LOCUS_LIST[[
      length(SUSIE_LOCUS_LIST) + 1L
    ]] <- data.table(
      locus_uid = rr$locus_uid,
      SuSiE_tested = TRUE,
      SuSiE_test_ids =
        collapse_unique(
          sz$test_id
        ),
      SuSiE_class =
        as.character(
          ss$SuSiE_class
        ),
      SuSiE_max_H4_default =
        as_num(
          ss$max_PP_H4_default
        ),
      SuSiE_max_H4_conservative =
        as_num(
          ss$max_PP_H4_conservative
        ),
      SuSiE_max_H3_default =
        as_num(
          ss$max_PP_H3_default
        ),
      V4_audited =
        any(
          sz$V4_audited,
          na.rm = TRUE
        ),
      V4_audit_labels =
        collapse_unique(
          sz$V4_audit_label
        )
    )
  }
}

COLOC_LOCUS_MAP <- rbindlist(
  COLOC_LOCUS_LIST,
  fill = TRUE
)

SUSIE_LOCUS_MAP <- rbindlist(
  SUSIE_LOCUS_LIST,
  fill = TRUE
)

# ============================================================
# 16. BUILD CORE LOCUS MASTER
# ============================================================

LOCUS <- merge(
  FUMA_LOCI,
  LAVA_MAP,
  by = "locus_uid",
  all.x = TRUE
)

LOCUS <- merge(
  LOCUS,
  TARGET_MAP,
  by = "locus_uid",
  all.x = TRUE
)

LOCUS <- merge(
  LOCUS,
  COLOC_LOCUS_MAP,
  by = "locus_uid",
  all.x = TRUE
)

LOCUS <- merge(
  LOCUS,
  SUSIE_LOCUS_MAP,
  by = "locus_uid",
  all.x = TRUE
)

# FUMA mapped gene counts.
GENE_COUNTS <- FUMA_GENES[
  ,
  .(
    FUMA_mapped_gene_n =
      uniqueN(
        gene
      ),
    FUMA_positional_gene_n =
      uniqueN(
        gene[
          as_num(posMapSNPs) >
            0
        ]
      ),
    FUMA_eqtl_mapped_gene_n =
      uniqueN(
        gene[
          as_num(eqtlMapSNPs) >
            0
        ]
      )
  ),
  by = locus_uid
]

LOCUS <- merge(
  LOCUS,
  GENE_COUNTS,
  by = "locus_uid",
  all.x = TRUE
)

for (cc in c(
  "FUMA_mapped_gene_n",
  "FUMA_positional_gene_n",
  "FUMA_eqtl_mapped_gene_n"
)) {
  LOCUS[
    is.na(
      get(cc)
    ),
    (cc) := 0L
  ]
}

# ============================================================
# 17. READ FINAL eQTL / sQTL SMR
# ============================================================

EQTL <- fread(
  EQTL_FILE,
  showProgress = FALSE
)

SQTL <- fread(
  SQTL_FILE,
  showProgress = FALSE
)

required_smr <- c(
  "trait",
  "Gene",
  "Probe_Chr",
  "Probe_bp",
  "p_SMR",
  "p_HEIDI",
  "Primary_pass"
)

if (!all(
  required_smr %in%
  names(EQTL)
)) {
  stop(
    "FINAL eQTL SMR master lacks required columns."
  )
}

if (!all(
  required_smr %in%
  names(SQTL)
)) {
  stop(
    "FINAL sQTL SMR master lacks required columns."
  )
}

standardize_smr <- function(x, layer) {

  z <- data.table(
    layer = layer,
    trait = as.character(
      x$trait
    ),
    tissue =
      if (
        "tissue" %in%
        names(x)
      ) {
        as.character(
          x$tissue
        )
      } else {
        NA_character_
      },
    gene = clean_gene(
      x$Gene
    ),
    ProbeID =
      if (
        "ProbeID" %in%
        names(x)
      ) {
        as.character(
          x$ProbeID
        )
      } else {
        NA_character_
      },
    chr = as_int(
      gsub(
        "^chr",
        "",
        as.character(
          x$Probe_Chr
        ),
        ignore.case = TRUE
      )
    ),
    bp = as_num(
      x$Probe_bp
    ),
    p_SMR = as_num(
      x$p_SMR
    ),
    p_HEIDI = as_num(
      x$p_HEIDI
    ),
    Primary_pass =
      as_logical_safe(
        x$Primary_pass
      )
  )

  z[
    Primary_pass == TRUE &
    !is.na(gene)
  ]
}

EQTL_P <- standardize_smr(
  EQTL,
  "eQTL"
)

SQTL_P <- standardize_smr(
  SQTL,
  "sQTL"
)

fwrite(
  EQTL_P,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_eQTL_primary_final.csv"
  )
)

fwrite(
  SQTL_P,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SourceData_sQTL_primary_final.csv"
  )
)

# ============================================================
# 18. MAP PRIMARY REGULATORY EVIDENCE TO FUMA GENE-LOCUS PAIRS
# ============================================================

# Keep only columns required for gene-locus mapping.
GENE_LOCUS <- unique(
  FUMA_GENES[
    ,
    .(
      locus_uid,
      pair_key,
      pair_display,
      trait1,
      trait2,
      GenomicLocus,
      gene,
      ensg =
        if (
          "ensg" %in%
          names(FUMA_GENES)
        ) {
          as.character(ensg)
        } else {
          NA_character_
        },
      gene_chr =
        as_int(chr),
      gene_start =
        as_num(start),
      gene_end =
        as_num(end),
      locus_chr =
        as_int(locus_chr),
      locus_start =
        as_num(locus_start),
      locus_stop =
        as_num(locus_stop),
      locus_lead =
        as.character(locus_lead),
      FUMA_positional =
        as.integer(
          as_num(posMapSNPs) >
            0
        ),
      FUMA_eqtl_mapped =
        as.integer(
          as_num(eqtlMapSNPs) >
            0
        ),
      FUMA_posMapSNPs =
        as_num(posMapSNPs),
      FUMA_eqtlMapSNPs =
        as_num(eqtlMapSNPs),
      FUMA_eqtlMapminP =
        as_num(eqtlMapminP),
      FUMA_minGwasP =
        as_num(minGwasP)
    )
  ],
  by = c(
    "locus_uid",
    "gene"
  )
)

REG_MAP_LIST <- list()

for (i in seq_len(nrow(GENE_LOCUS))) {

  rr <- GENE_LOCUS[i]

  traits_pair <- c(
    rr$trait1,
    rr$trait2
  )

  # --------------------
  # eQTL
  # --------------------
  e <- EQTL_P[
    gene ==
      rr$gene &
    trait %in%
      traits_pair &
    chr ==
      rr$locus_chr &
    bp >=
      rr$locus_start &
    bp <=
      rr$locus_stop
  ]

  # --------------------
  # sQTL
  # --------------------
  s <- SQTL_P[
    gene ==
      rr$gene &
    trait %in%
      traits_pair &
    chr ==
      rr$locus_chr &
    bp >=
      rr$locus_start &
    bp <=
      rr$locus_stop
  ]

  e1 <- e[
    trait ==
      rr$trait1
  ]

  e2 <- e[
    trait ==
      rr$trait2
  ]

  s1 <- s[
    trait ==
      rr$trait1
  ]

  s2 <- s[
    trait ==
      rr$trait2
  ]

  # strict shared splicing event:
  # same tissue + same ProbeID across both traits.
  shared_s_event_n <- 0L

  if (
    nrow(s1) &&
    nrow(s2)
  ) {

    k1 <- unique(
      paste(
        s1$tissue,
        s1$ProbeID,
        sep = "|"
      )
    )

    k2 <- unique(
      paste(
        s2$tissue,
        s2$ProbeID,
        sep = "|"
      )
    )

    shared_s_event_n <- length(
      intersect(
        k1,
        k2
      )
    )
  }

  REG_MAP_LIST[[
    length(REG_MAP_LIST) + 1L
  ]] <- data.table(
    locus_uid = rr$locus_uid,
    gene = rr$gene,

    eQTL_trait1 =
      as.integer(
        nrow(e1) > 0
      ),

    eQTL_trait2 =
      as.integer(
        nrow(e2) > 0
      ),

    eQTL_shared_gene =
      as.integer(
        nrow(e1) > 0 &
        nrow(e2) > 0
      ),

    eQTL_min_p_trait1 =
      safe_min(
        e1$p_SMR
      ),

    eQTL_min_p_trait2 =
      safe_min(
        e2$p_SMR
      ),

    eQTL_tissues_trait1 =
      collapse_unique(
        e1$tissue
      ),

    eQTL_tissues_trait2 =
      collapse_unique(
        e2$tissue
      ),

    sQTL_trait1 =
      as.integer(
        nrow(s1) > 0
      ),

    sQTL_trait2 =
      as.integer(
        nrow(s2) > 0
      ),

    sQTL_shared_gene =
      as.integer(
        nrow(s1) > 0 &
        nrow(s2) > 0
      ),

    sQTL_shared_event_n =
      shared_s_event_n,

    sQTL_min_p_trait1 =
      safe_min(
        s1$p_SMR
      ),

    sQTL_min_p_trait2 =
      safe_min(
        s2$p_SMR
      ),

    sQTL_tissues_trait1 =
      collapse_unique(
        s1$tissue
      ),

    sQTL_tissues_trait2 =
      collapse_unique(
        s2$tissue
      )
  )
}

REG_MAP <- rbindlist(
  REG_MAP_LIST,
  fill = TRUE
)

# ============================================================
# 19. GENE-LOCUS MASTER
# ============================================================

GENE_LOCUS <- merge(
  GENE_LOCUS,
  REG_MAP,
  by = c(
    "locus_uid",
    "gene"
  ),
  all.x = TRUE
)

LOCUS_ATTACH <- LOCUS[
  ,
  .(
    locus_uid,
    pair_key,
    pair_display,
    locus_chr = as_int(chr),
    locus_start = as_num(start),
    locus_stop = as_num(end),
    lead_snp = as.character(rsID),
    conjFDR_value,
    LAVA_tested,
    LAVA_pair_bonf_sig,
    LAVA_min_p,
    LAVA_rho_at_min_p,
    LAVA_direction_at_min_p,
    target_regions,
    coloc_tested,
    coloc_max_PP_H4,
    coloc_max_PP_H3,
    coloc_ABF_class,
    coloc_interpretation,
    SuSiE_tested,
    SuSiE_class,
    SuSiE_max_H4_default,
    SuSiE_max_H4_conservative,
    SuSiE_max_H3_default,
    V4_audited,
    V4_audit_labels
  )
]

GENE_LOCUS <- merge(
  GENE_LOCUS,
  LOCUS_ATTACH,
  by = c(
    "locus_uid",
    "pair_key",
    "pair_display"
  ),
  all.x = TRUE,
  suffixes = c(
    "",
    "_locus"
  )
)

# Evidence-structure descriptors.
GENE_LOCUS[
  ,
  regulatory_layer_n :=
    as.integer(
      eQTL_trait1 +
      eQTL_trait2 >
        0
    ) +
    as.integer(
      sQTL_trait1 +
      sQTL_trait2 >
        0
    )
]

GENE_LOCUS[
  ,
  regulatory_trait_n :=
    as.integer(
      eQTL_trait1 +
      sQTL_trait1 >
        0
    ) +
    as.integer(
      eQTL_trait2 +
      sQTL_trait2 >
        0
    )
]

GENE_LOCUS[
  ,
  shared_regulatory_gene :=
    as.integer(
      eQTL_shared_gene == 1 |
      sQTL_shared_gene == 1
    )
]

GENE_LOCUS[
  ,
  FUMA_mapping_mode_n :=
    FUMA_positional +
    FUMA_eqtl_mapped
]

GENE_LOCUS[
  ,
  locus_genetic_layer_n :=
    1L + # conjFDR / FUMA locus is the starting layer
    as.integer(
      LAVA_pair_bonf_sig
    ) +
    as.integer(
      coloc_tested
    ) +
    as.integer(
      SuSiE_tested
    )
]

# Descriptive candidate class; NOT causal proof.
GENE_LOCUS[
  ,
  candidate_pattern :=
    fifelse(
      shared_regulatory_gene == 1 &
      regulatory_layer_n == 2,
      "Shared_gene_with_eQTL_and_sQTL",

      fifelse(
        shared_regulatory_gene == 1,
        "Shared_gene_with_one_regulatory_layer",

        fifelse(
          regulatory_trait_n == 2,
          "Both_traits_regulatory_nonshared_layer",

          fifelse(
            regulatory_layer_n >= 1,
            "One_trait_regulatory_support",
            "FUMA_mapping_only"
          )
        )
      )
    )
]

# Best regulatory p for lexicographic ordering only.
GENE_LOCUS[
  ,
  min_primary_regulatory_p :=
    apply(
      cbind(
        eQTL_min_p_trait1,
        eQTL_min_p_trait2,
        sQTL_min_p_trait1,
        sQTL_min_p_trait2
      ),
      1,
      function(z) {
        z <- as.numeric(z)
        z <- z[
          is.finite(z)
        ]
        if (!length(z)) {
          return(NA_real_)
        }
        min(z)
      }
    )
]

# ============================================================
# 20. LOCUS-LEVEL REGULATORY SUMMARY
# ============================================================

LOCUS_REG <- GENE_LOCUS[
  ,
  .(
    eQTL_gene_n_trait1 =
      sum(
        eQTL_trait1,
        na.rm = TRUE
      ),

    eQTL_gene_n_trait2 =
      sum(
        eQTL_trait2,
        na.rm = TRUE
      ),

    eQTL_shared_gene_n =
      sum(
        eQTL_shared_gene,
        na.rm = TRUE
      ),

    sQTL_gene_n_trait1 =
      sum(
        sQTL_trait1,
        na.rm = TRUE
      ),

    sQTL_gene_n_trait2 =
      sum(
        sQTL_trait2,
        na.rm = TRUE
      ),

    sQTL_shared_gene_n =
      sum(
        sQTL_shared_gene,
        na.rm = TRUE
      ),

    sQTL_shared_event_n =
      sum(
        sQTL_shared_event_n,
        na.rm = TRUE
      ),

    shared_regulatory_genes =
      collapse_unique(
        gene[
          shared_regulatory_gene ==
            1
        ]
      )
  ),
  by = locus_uid
]

LOCUS <- merge(
  LOCUS,
  LOCUS_REG,
  by = "locus_uid",
  all.x = TRUE
)

for (cc in c(
  "eQTL_gene_n_trait1",
  "eQTL_gene_n_trait2",
  "eQTL_shared_gene_n",
  "sQTL_gene_n_trait1",
  "sQTL_gene_n_trait2",
  "sQTL_shared_gene_n",
  "sQTL_shared_event_n"
)) {

  LOCUS[
    is.na(
      get(cc)
    ),
    (cc) := 0L
  ]
}

LOCUS[
  is.na(
    shared_regulatory_genes
  ),
  shared_regulatory_genes := ""
]

# ============================================================
# 21. FINAL LOCUS INTERPRETATION
# ============================================================
#
# This field is descriptive and conservative.
# It prioritizes formal SuSiE when available.
# ============================================================

LOCUS[
  ,
  locus_interpretation :=
    fifelse(
      SuSiE_tested &
      SuSiE_class ==
        "H4_robust_both_priors",
      "Shared causal-signal evidence (robust SuSiE)",

      fifelse(
        SuSiE_tested &
        SuSiE_class ==
          "H4_prior_sensitive",
        "Shared-signal evidence under default prior; prior-sensitive",

        fifelse(
          SuSiE_tested &
          SuSiE_class ==
            "H3_dominant",
          "Distinct causal signals (H3-dominant)",

          fifelse(
            SuSiE_tested &
            SuSiE_class ==
              "No_signal_pair",
            "Fine-mapping unresolved / signal-pair unavailable",

            fifelse(
              SuSiE_tested,
              "No dominant H3/H4 SuSiE conclusion",

              fifelse(
                LAVA_pair_bonf_sig,
                "Significant local genetic correlation; no targeted SuSiE test",
                "conjFDR/FUMA locus without additional strict local evidence"
              )
            )
          )
        )
      )
    )
]

# ============================================================
# 22. BUILD THE EIGHT FORMAL TARGET-TEST TABLE
# ============================================================

FORMAL_TESTS[
  ,
  pair_key :=
    mapply(
      canon_pair,
      trait1,
      trait2
    )
]

FORMAL_TESTS[
  ,
  test_id :=
    paste(
      region,
      trait1,
      trait2,
      sep = "__"
    )
]

FOCAL_LIST <- list()

for (i in seq_len(nrow(FORMAL_TESTS))) {

  ft <- FORMAL_TESTS[i]

  reg <- TARGET_REGIONS[
    region ==
      ft$region
  ]

  # FUMA loci overlapping this target region for the same pair.
  fl <- LOCUS[
    pair_key ==
      ft$pair_key &
    as_int(chr) ==
      as_int(reg$chr) &
    interval_overlap(
      as_num(start),
      as_num(end),
      as_num(reg$start),
      as_num(reg$stop)
    )
  ]

  cz <- COLOC[
    test_id ==
      ft$test_id
  ]

  sz <- SUSIE_V3[
    test_id ==
      ft$test_id
  ]

  # LAVA rows directly overlapping target region.
  lz <- LAVA[
    pair_key_std ==
      ft$pair_key &
    as_int(chr) ==
      as_int(reg$chr) &
    interval_overlap(
      as_num(start),
      as_num(stop),
      as_num(reg$start),
      as_num(reg$stop)
    )
  ]

  # Regulatory hits across target window, not only the exact FUMA interval.
  e1 <- EQTL_P[
    trait ==
      ft$trait1 &
    chr ==
      reg$chr &
    bp >=
      reg$start &
    bp <=
      reg$stop
  ]

  e2 <- EQTL_P[
    trait ==
      ft$trait2 &
    chr ==
      reg$chr &
    bp >=
      reg$start &
    bp <=
      reg$stop
  ]

  s1 <- SQTL_P[
    trait ==
      ft$trait1 &
    chr ==
      reg$chr &
    bp >=
      reg$start &
    bp <=
      reg$stop
  ]

  s2 <- SQTL_P[
    trait ==
      ft$trait2 &
    chr ==
      reg$chr &
    bp >=
      reg$start &
    bp <=
      reg$stop
  ]

  shared_e_gene <- intersect(
    unique(e1$gene),
    unique(e2$gene)
  )

  shared_s_gene <- intersect(
    unique(s1$gene),
    unique(s2$gene)
  )

  shared_s_event <- intersect(
    unique(
      paste(
        s1$tissue,
        s1$ProbeID,
        sep = "|"
      )
    ),
    unique(
      paste(
        s2$tissue,
        s2$ProbeID,
        sep = "|"
      )
    )
  )

  # LAVA summary.
  lava_sig <- if (nrow(lz)) {
    any(
      lz$pass_pair_bonferroni_std,
      na.rm = TRUE
    )
  } else {
    FALSE
  }

  lava_min_p <- if (nrow(lz)) {
    safe_min(
      lz$p
    )
  } else {
    NA_real_
  }

  lava_rho <- NA_real_

  if (
    nrow(lz) &&
    any(
      is.finite(
        lz$p
      )
    )
  ) {
    ib <- which.min(
      lz$p
    )
    lava_rho <- as_num(
      lz$rho[ib]
    )
  }

  FOCAL_LIST[[
    length(FOCAL_LIST) + 1L
  ]] <- data.table(
    test_id = ft$test_id,
    region = ft$region,
    pair_key = ft$pair_key,
    pair_display = pair_to_display(
      ft$pair_key
    ),
    trait1 = ft$trait1,
    trait2 = ft$trait2,
    chr = reg$chr,
    start = reg$start,
    stop = reg$stop,

    conjFDR_FUMA_locus_present =
      nrow(fl) > 0,

    n_FUMA_loci_overlap =
      nrow(fl),

    FUMA_lead_SNPs =
      collapse_unique(
        fl$rsID
      ),

    FUMA_mapped_gene_n =
      if (
        nrow(fl)
      ) {
        sum(
          fl$FUMA_mapped_gene_n,
          na.rm = TRUE
        )
      } else {
        0L
      },

    LAVA_tested =
      nrow(lz) > 0,

    LAVA_pair_bonf_sig =
      lava_sig,

    LAVA_min_p =
      lava_min_p,

    LAVA_rho_at_min_p =
      lava_rho,

    coloc_tested =
      nrow(cz) > 0,

    coloc_PP_H3 =
      if (
        nrow(cz)
      ) {
        as_num(
          cz$`PP.H3`[1]
        )
      } else {
        NA_real_
      },

    coloc_PP_H4 =
      if (
        nrow(cz)
      ) {
        as_num(
          cz$`PP.H4`[1]
        )
      } else {
        NA_real_
      },

    coloc_class =
      if (
        nrow(cz)
      ) {
        as.character(
          cz$ABF_class[1]
        )
      } else {
        "Not_tested"
      },

    SuSiE_tested =
      nrow(sz) > 0,

    SuSiE_class =
      if (
        nrow(sz)
      ) {
        as.character(
          sz$SuSiE_class[1]
        )
      } else {
        "Not_tested"
      },

    SuSiE_H4_default =
      if (
        nrow(sz)
      ) {
        as_num(
          sz$max_PP_H4_default[1]
        )
      } else {
        NA_real_
      },

    SuSiE_H4_conservative =
      if (
        nrow(sz)
      ) {
        as_num(
          sz$max_PP_H4_conservative[1]
        )
      } else {
        NA_real_
      },

    SuSiE_H3_default =
      if (
        nrow(sz)
      ) {
        as_num(
          sz$max_PP_H3_default[1]
        )
      } else {
        NA_real_
      },

    V4_SER_audited =
      if (
        nrow(sz)
      ) {
        isTRUE(
          sz$V4_audited[1]
        )
      } else {
        FALSE
      },

    V4_audit_label =
      if (
        nrow(sz)
      ) {
        as.character(
          sz$V4_audit_label[1]
        )
      } else {
        ""
      },

    eQTL_trait1_gene_n =
      uniqueN(
        e1$gene
      ),

    eQTL_trait2_gene_n =
      uniqueN(
        e2$gene
      ),

    eQTL_shared_gene_n =
      length(
        shared_e_gene
      ),

    eQTL_shared_genes =
      collapse_unique(
        shared_e_gene
      ),

    sQTL_trait1_gene_n =
      uniqueN(
        s1$gene
      ),

    sQTL_trait2_gene_n =
      uniqueN(
        s2$gene
      ),

    sQTL_shared_gene_n =
      length(
        shared_s_gene
      ),

    sQTL_shared_genes =
      collapse_unique(
        shared_s_gene
      ),

    sQTL_shared_event_n =
      length(
        shared_s_event
      )
  )
}

FOCAL <- rbindlist(
  FOCAL_LIST,
  fill = TRUE
)

# ============================================================
# 23. FORMAL TARGET-TEST INTERPRETATION
# ============================================================

FOCAL[
  ,
  final_test_interpretation :=
    fifelse(
      SuSiE_class ==
        "H4_robust_both_priors",
      "Robust shared causal-signal evidence",

      fifelse(
        SuSiE_class ==
          "H4_prior_sensitive",
        "Shared-signal evidence under default prior; prior-sensitive",

        fifelse(
          SuSiE_class ==
            "H3_dominant",
          "Distinct causal signals (H3-dominant)",

          fifelse(
            SuSiE_class ==
              "No_signal_pair",
            "Fine-mapping unresolved / signal-pair unavailable",

            fifelse(
              SuSiE_class ==
                "No_dominant_H3_H4",
              "No dominant H3/H4 multi-signal conclusion",
              "Not formally tested by SuSiE"
            )
          )
        )
      )
    )
]

# ============================================================
# 24. LEXICOGRAPHIC CANDIDATE ORDER — NO ARBITRARY SCORE
# ============================================================

# Use evidence structure, not a summed causal score.
GENE_LOCUS[
  ,
  rank_shared_reg :=
    -shared_regulatory_gene
]

GENE_LOCUS[
  ,
  rank_reg_layers :=
    -regulatory_layer_n
]

GENE_LOCUS[
  ,
  rank_reg_traits :=
    -regulatory_trait_n
]

GENE_LOCUS[
  ,
  rank_shared_event :=
    -as.integer(
      sQTL_shared_event_n >
        0
    )
]

GENE_LOCUS[
  ,
  rank_fuma_modes :=
    -FUMA_mapping_mode_n
]

GENE_LOCUS[
  ,
  rank_lava :=
    -as.integer(
      LAVA_pair_bonf_sig
    )
]

GENE_LOCUS[
  ,
  rank_susie :=
    fifelse(
      SuSiE_class ==
        "H4_robust_both_priors",
      0L,
      fifelse(
        SuSiE_class ==
          "H4_prior_sensitive",
        1L,
        fifelse(
          SuSiE_class ==
            "H3_dominant",
          3L,
          fifelse(
            SuSiE_tested,
            2L,
            4L
          )
        )
      )
    )
]

setorder(
  GENE_LOCUS,
  rank_shared_reg,
  rank_reg_layers,
  rank_reg_traits,
  rank_shared_event,
  rank_susie,
  rank_lava,
  rank_fuma_modes,
  min_primary_regulatory_p,
  gene
)

GENE_LOCUS[
  ,
  Evidence_rank :=
    seq_len(.N)
]

# Remove internal ordering helpers from publication table.
RANK_HELPERS <- c(
  "rank_shared_reg",
  "rank_reg_layers",
  "rank_reg_traits",
  "rank_shared_event",
  "rank_fuma_modes",
  "rank_lava",
  "rank_susie"
)

# Candidate publication table:
# focus on genes with some strict regulatory evidence OR
# multi-mode FUMA mapping in genetically informative loci.
PRIORITY_CANDIDATES <- GENE_LOCUS[
  regulatory_layer_n > 0 |
  shared_regulatory_gene == 1 |
  FUMA_mapping_mode_n >= 2
]

PRIORITY_CANDIDATES[
  ,
  (RANK_HELPERS) := NULL
]

# Top 50 by evidence-structure order.
TOP50 <- PRIORITY_CANDIDATES[
  1:min(
    50L,
    .N
  )
]

# ============================================================
# 25. MAIN / SUPPLEMENTARY TABLE EXPORTS
# ============================================================

TABLE1_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table1_focal_8test_multilayer_evidence.csv"
)

fwrite(
  FOCAL,
  TABLE1_FILE
)

TABLE2_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table2_priority_locus_gene_candidates_TOP50.csv"
)

fwrite(
  TOP50,
  TABLE2_FILE
)

LOCUS_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS1_all_FUMA_loci_multilayer_evidence.csv"
)

fwrite(
  LOCUS,
  LOCUS_FILE
)

GENE_LOCUS_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS2_all_locus_gene_evidence.csv"
)

GENE_LOCUS_EXPORT <- copy(
  GENE_LOCUS
)

GENE_LOCUS_EXPORT[
  ,
  (RANK_HELPERS) := NULL
]

fwrite(
  GENE_LOCUS_EXPORT,
  GENE_LOCUS_FILE
)

FUMA_GENE_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS3_FUMA_gene_mapping_final.csv"
)

fwrite(
  FUMA_GENES,
  FUMA_GENE_FILE
)

FORMAL_SUSIE_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS4_formal_SuSiE_V3_with_V4_audit.csv"
)

fwrite(
  SUSIE_V3,
  FORMAL_SUSIE_FILE
)

# ============================================================
# 26. MAIN FIGURE 1 SOURCE DATA
#     EIGHT FORMAL TESTS × MULTI-LAYER EVIDENCE
# ============================================================

FIG1_LIST <- list()

for (i in seq_len(nrow(FOCAL))) {

  r <- FOCAL[i]

  row_label <- paste0(
    r$region,
    " | ",
    r$pair_display
  )

  # --------------------
  # conjFDR / FUMA locus
  # --------------------
  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "conjFDR locus",
    status =
      if (
        r$conjFDR_FUMA_locus_present
      ) {
        "Shared_locus"
      } else {
        "No_overlap"
      },
    plot_code =
      if (
        r$conjFDR_FUMA_locus_present
      ) {
        2L
      } else {
        0L
      },
    label =
      if (
        r$conjFDR_FUMA_locus_present
      ) {
        "Locus"
      } else {
        "—"
      }
  )

  # --------------------
  # LAVA
  # --------------------
  lava_status <-
    if (!r$LAVA_tested) {
      "Not_tested"
    } else if (r$LAVA_pair_bonf_sig) {
      ifelse(
        is.finite(r$LAVA_rho_at_min_p) &&
        r$LAVA_rho_at_min_p < 0,
        "Significant_negative",
        "Significant_positive"
      )
    } else {
      "Tested_not_significant"
    }

  lava_code <-
    if (lava_status == "Not_tested") {
      NA_integer_
    } else if (
      lava_status %in%
      c(
        "Significant_positive",
        "Significant_negative"
      )
    ) {
      2L
    } else {
      0L
    }

  lava_label <-
    if (
      lava_status ==
        "Significant_positive"
    ) {
      "rg+"
    } else if (
      lava_status ==
        "Significant_negative"
    ) {
      "rg−"
    } else if (
      lava_status ==
        "Tested_not_significant"
    ) {
      "NS"
    } else {
      "NA"
    }

  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "LAVA",
    status = lava_status,
    plot_code = lava_code,
    label = lava_label
  )

  # --------------------
  # coloc ABF
  # --------------------
  abf_status <-
    if (!r$coloc_tested) {
      "Not_tested"
    } else if (
      r$coloc_class ==
        "H4_strong"
    ) {
      "H4_support"
    } else if (
      r$coloc_class ==
        "H3_strong"
    ) {
      "H3_support"
    } else {
      "Non_dominant"
    }

  abf_code <-
    if (
      abf_status ==
        "Not_tested"
    ) {
      NA_integer_
    } else if (
      abf_status ==
        "H4_support"
    ) {
      2L
    } else if (
      abf_status ==
        "H3_support"
    ) {
      -1L
    } else {
      1L
    }

  abf_label <-
    if (
      abf_status ==
        "H4_support"
    ) {
      "H4"
    } else if (
      abf_status ==
        "H3_support"
    ) {
      "H3"
    } else if (
      abf_status ==
        "Non_dominant"
    ) {
      "Amb"
    } else {
      "NA"
    }

  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "coloc ABF",
    status = abf_status,
    plot_code = abf_code,
    label = abf_label
  )

  # --------------------
  # formal SuSiE
  # --------------------
  sus_status <- r$SuSiE_class

  sus_code <-
    if (
      sus_status ==
        "H4_robust_both_priors"
    ) {
      2L
    } else if (
      sus_status ==
        "H4_prior_sensitive"
    ) {
      1L
    } else if (
      sus_status ==
        "H3_dominant"
    ) {
      -1L
    } else if (
      sus_status %in%
        c(
          "No_signal_pair",
          "No_dominant_H3_H4"
        )
    ) {
      0L
    } else {
      NA_integer_
    }

  sus_label <-
    if (
      sus_status ==
        "H4_robust_both_priors"
    ) {
      "H4"
    } else if (
      sus_status ==
        "H4_prior_sensitive"
    ) {
      "H4*"
    } else if (
      sus_status ==
        "H3_dominant"
    ) {
      "H3"
    } else if (
      sus_status ==
        "No_signal_pair"
    ) {
      "No CS"
    } else if (
      sus_status ==
        "No_dominant_H3_H4"
    ) {
      "ND"
    } else {
      "NA"
    }

  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "SuSiE",
    status = sus_status,
    plot_code = sus_code,
    label = sus_label
  )

  # --------------------
  # eQTL-SMR
  # --------------------
  eqtl_status <-
    if (
      r$eQTL_shared_gene_n > 0
    ) {
      "Shared_gene"
    } else if (
      r$eQTL_trait1_gene_n > 0 |
      r$eQTL_trait2_gene_n > 0
    ) {
      "One_trait_only"
    } else {
      "No_strict_signal"
    }

  eqtl_code <-
    if (
      eqtl_status ==
        "Shared_gene"
    ) {
      2L
    } else if (
      eqtl_status ==
        "One_trait_only"
    ) {
      1L
    } else {
      0L
    }

  eqtl_label <-
    if (
      eqtl_status ==
        "Shared_gene"
    ) {
      paste0(
        "Shared(",
        r$eQTL_shared_gene_n,
        ")"
      )
    } else if (
      eqtl_status ==
        "One_trait_only"
    ) {
      "1-trait"
    } else {
      "—"
    }

  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "eQTL-SMR",
    status = eqtl_status,
    plot_code = eqtl_code,
    label = eqtl_label
  )

  # --------------------
  # sQTL-SMR
  # --------------------
  sqtl_status <-
    if (
      r$sQTL_shared_gene_n > 0
    ) {
      "Shared_gene"
    } else if (
      r$sQTL_trait1_gene_n > 0 |
      r$sQTL_trait2_gene_n > 0
    ) {
      "One_trait_only"
    } else {
      "No_strict_signal"
    }

  sqtl_code <-
    if (
      sqtl_status ==
        "Shared_gene"
    ) {
      2L
    } else if (
      sqtl_status ==
        "One_trait_only"
    ) {
      1L
    } else {
      0L
    }

  sqtl_label <-
    if (
      sqtl_status ==
        "Shared_gene"
    ) {
      paste0(
        "Shared(",
        r$sQTL_shared_gene_n,
        ")"
      )
    } else if (
      sqtl_status ==
        "One_trait_only"
    ) {
      "1-trait"
    } else {
      "—"
    }

  FIG1_LIST[[
    length(FIG1_LIST) + 1L
  ]] <- data.table(
    row_label = row_label,
    test_id = r$test_id,
    evidence_layer = "sQTL-SMR",
    status = sqtl_status,
    plot_code = sqtl_code,
    label = sqtl_label
  )
}

FIG1 <- rbindlist(
  FIG1_LIST,
  fill = TRUE
)

FIG1[
  ,
  evidence_layer :=
    factor(
      evidence_layer,
      levels = c(
        "conjFDR locus",
        "LAVA",
        "coloc ABF",
        "SuSiE",
        "eQTL-SMR",
        "sQTL-SMR"
      )
    )
]

row_order <- paste0(
  FORMAL_TESTS$region,
  " | ",
  vapply(
    FORMAL_TESTS$pair_key,
    pair_to_display,
    character(1)
  )
)

FIG1[
  ,
  row_label :=
    factor(
      row_label,
      levels = rev(
        row_order
      )
    )
]

FIG1_SOURCE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_Figure1_focal_multilayer_evidence_matrix.csv"
)

fwrite(
  FIG1,
  FIG1_SOURCE
)

# Figure legend definition stored separately, not drawn.
FIG1_LEGEND <- data.table(
  plot_code = c(
    2,
    1,
    0,
    -1,
    NA
  ),
  meaning = c(
    "Convergent/shared or significant supporting evidence",
    "Partial/prior-sensitive/one-trait evidence",
    "Tested without strict supporting evidence",
    "Evidence favoring distinct causal signals (H3)",
    "Not tested / not evaluable"
  )
)

fwrite(
  FIG1_LEGEND,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "Legend_Figure1_evidence_code.csv"
  )
)

# ============================================================
# 27. MAIN FIGURE 1 — FOCAL MULTI-LAYER EVIDENCE MATRIX
# ============================================================

FIG1[
  ,
  code_factor :=
    factor(
      plot_code,
      levels = c(
        -1,
        0,
        1,
        2
      )
    )
]

p1 <- ggplot(
  FIG1,
  aes(
    x = evidence_layer,
    y = row_label
  )
) +
  geom_tile(
    aes(
      fill = code_factor
    ),
    colour = "#D2D5DA",
    linewidth = 0.45,
    width = 0.94,
    height = 0.88
  ) +
  geom_text(
    aes(
      label = label
    ),
    family = "Arial",
    size = 3.25,
    colour = "#111111"
  ) +
  scale_fill_manual(
    values = c(
      "-1" = "#D8A0A0",
      "0" = "#EFEFEF",
      "1" = "#D9C48C",
      "2" = "#7CA7A0"
    ),
    na.value = "white",
    drop = FALSE
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0.4,
      colour = "black",
      size = 10
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 9.2
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      14,
      18,
      14,
      12
    )
  )

FIG1_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure1_focal_multilayer_evidence_matrix.tiff"
)

ggsave(
  FIG1_TIFF,
  p1,
  width = 7.4,
  height = 5.9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

FIG1_PDF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure1_focal_multilayer_evidence_matrix.pdf"
)

ggsave(
  FIG1_PDF,
  p1,
  width = 7.4,
  height = 5.9,
  units = "in",
  device = cairo_pdf
)

# ============================================================
# 28. MAIN FIGURE 2 SOURCE DATA
#     TOP GENE-LOCUS CANDIDATES
# ============================================================

# Select at most 30 distinct gene-locus candidates by the
# lexicographic evidence ordering already frozen above.
TOP30_GENE_LOCUS <- TOP50[
  1:min(
    30L,
    .N
  )
]

TOP30_GENE_LOCUS[
  ,
  row_label := paste0(
    gene,
    " | ",
    pair_display,
    " | chr",
    locus_chr,
    ":",
    format(
      locus_start,
      scientific = FALSE,
      trim = TRUE
    )
  )
]

# Matrix columns are gene/locus-linked evidence.
FIG2_LIST <- list()

for (i in seq_len(nrow(TOP30_GENE_LOCUS))) {

  r <- TOP30_GENE_LOCUS[i]

  add_row <- function(
    layer,
    present,
    label
  ) {

    FIG2_LIST[[
      length(FIG2_LIST) + 1L
    ]] <<- data.table(
      row_label = r$row_label,
      gene = r$gene,
      locus_uid = r$locus_uid,
      pair_display = r$pair_display,
      evidence_layer = layer,
      present = as.integer(
        present
      ),
      label = label
    )
  }

  add_row(
    "FUMA-pos",
    r$FUMA_positional == 1,
    ifelse(
      r$FUMA_positional == 1,
      "●",
      ""
    )
  )

  add_row(
    "FUMA-eQTL",
    r$FUMA_eqtl_mapped == 1,
    ifelse(
      r$FUMA_eqtl_mapped == 1,
      "●",
      ""
    )
  )

  add_row(
    paste0(
      "eQTL-",
      r$trait1
    ),
    r$eQTL_trait1 == 1,
    ifelse(
      r$eQTL_trait1 == 1,
      "●",
      ""
    )
  )

  add_row(
    paste0(
      "eQTL-",
      r$trait2
    ),
    r$eQTL_trait2 == 1,
    ifelse(
      r$eQTL_trait2 == 1,
      "●",
      ""
    )
  )

  add_row(
    paste0(
      "sQTL-",
      r$trait1
    ),
    r$sQTL_trait1 == 1,
    ifelse(
      r$sQTL_trait1 == 1,
      "●",
      ""
    )
  )

  add_row(
    paste0(
      "sQTL-",
      r$trait2
    ),
    r$sQTL_trait2 == 1,
    ifelse(
      r$sQTL_trait2 == 1,
      "●",
      ""
    )
  )

  add_row(
    "Shared-reg",
    r$shared_regulatory_gene == 1,
    ifelse(
      r$shared_regulatory_gene == 1,
      "Shared",
      ""
    )
  )
}

FIG2 <- rbindlist(
  FIG2_LIST,
  fill = TRUE
)

# Because trait-specific columns differ by pair, use a universal
# ordered set of possible layers.
fig2_levels <- c(
  "FUMA-pos",
  "FUMA-eQTL",
  "eQTL-AF",
  "eQTL-HFpEF",
  "eQTL-BMI",
  "eQTL-OSA",
  "sQTL-AF",
  "sQTL-HFpEF",
  "sQTL-BMI",
  "sQTL-OSA",
  "Shared-reg"
)

FIG2[
  ,
  evidence_layer :=
    factor(
      evidence_layer,
      levels = fig2_levels
    )
]

FIG2[
  ,
  row_label :=
    factor(
      row_label,
      levels = rev(
        TOP30_GENE_LOCUS$row_label
      )
    )
]

FIG2_SOURCE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_Figure2_top_gene_locus_regulatory_matrix.csv"
)

fwrite(
  FIG2,
  FIG2_SOURCE
)

# ============================================================
# 29. MAIN FIGURE 2 — GENE-LOCUS REGULATORY MATRIX
# ============================================================

p2 <- ggplot(
  FIG2,
  aes(
    x = evidence_layer,
    y = row_label
  )
) +
  geom_tile(
    fill = "#F4F4F4",
    colour = "#D9D9D9",
    linewidth = 0.35,
    width = 0.92,
    height = 0.88
  ) +
  geom_point(
    data = FIG2[
      present == 1
    ],
    aes(
      x = evidence_layer,
      y = row_label
    ),
    shape = 21,
    size = 3.7,
    stroke = 0.45,
    fill = "#6D8FA3",
    colour = "#263943"
  ) +
  scale_x_discrete(
    position = "top",
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10.5,
    base_family = "Arial"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0.4,
      colour = "black",
      size = 9.3
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 7.8
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      14,
      14,
      12,
      10
    )
  )

FIG2_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure2_top_gene_locus_regulatory_matrix.tiff"
)

ggsave(
  FIG2_TIFF,
  p2,
  width = 8.1,
  height = 8.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# 30. SUPPLEMENTARY FIGURE
#     ALL INFORMATIVE FUMA LOCI
# ============================================================

LOCUS_PLOT <- LOCUS[
  LAVA_pair_bonf_sig |
  coloc_tested |
  SuSiE_tested |
  eQTL_shared_gene_n > 0 |
  sQTL_shared_gene_n > 0
]

LOCUS_PLOT[
  ,
  row_label := paste0(
    pair_display,
    " | chr",
    chr,
    ":",
    start,
    "-",
    end
  )
]

LOCUS_PLOT[
  ,
  genetic_layers :=
    1L +
    as.integer(
      LAVA_pair_bonf_sig
    ) +
    as.integer(
      coloc_tested
    ) +
    as.integer(
      SuSiE_tested
    )
]

LOCUS_PLOT[
  ,
  regulatory_layers :=
    as.integer(
      eQTL_gene_n_trait1 +
      eQTL_gene_n_trait2 >
        0
    ) +
    as.integer(
      sQTL_gene_n_trait1 +
      sQTL_gene_n_trait2 >
        0
    )
]

setorder(
  LOCUS_PLOT,
  -genetic_layers,
  -regulatory_layers,
  conjFDR_value
)

if (nrow(LOCUS_PLOT) > 50L) {
  LOCUS_PLOT <- LOCUS_PLOT[
    1:50
  ]
}

LOCUS_PLOT[
  ,
  row_label :=
    factor(
      row_label,
      levels = rev(
        row_label
      )
    )
]

FIGS_SOURCE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_FigureS_all_informative_loci.csv"
)

fwrite(
  LOCUS_PLOT,
  FIGS_SOURCE
)

pS <- ggplot(
  LOCUS_PLOT,
  aes(
    x = -log10(
      pmax(
        conjFDR_value,
        .Machine$double.xmin
      )
    ),
    y = row_label
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = -log10(
        pmax(
          conjFDR_value,
          .Machine$double.xmin
        )
      ),
      y = row_label,
      yend = row_label
    ),
    linewidth = 0.5,
    colour = "#C6CBD0"
  ) +
  geom_point(
    aes(
      size = genetic_layers
    ),
    shape = 21,
    stroke = 0.4,
    fill = "#718FA2",
    colour = "#263943"
  ) +
  labs(
    x = expression(-log[10]("conjFDR")),
    y = NULL
  ) +
  theme_classic(
    base_size = 10.5,
    base_family = "Arial"
  ) +
  theme(
    axis.text.y = element_text(
      size = 7.3,
      colour = "black"
    ),
    axis.text.x = element_text(
      colour = "black"
    ),
    legend.position = "none",
    panel.border = element_blank()
  )

ggsave(
  file.path(
    OUT_ROOT,
    "05_SUPPLEMENTARY_FIGURES",
    "FigureS_informative_loci_landscape.tiff"
  ),
  pS,
  width = 7.4,
  height = 9.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# 31. QC SUMMARY
# ============================================================

QC <- data.table(
  item = c(
    "conjFDR_final_locus_rows",
    "FUMA_loci",
    "FUMA_mapped_gene_rows",
    "FUMA_unique_genes",
    "LAVA_rows",
    "coloc_ABF_tests",
    "formal_SuSiE_V3_tests",
    "V4_SER_audit_rows",
    "eQTL_primary_rows",
    "sQTL_primary_rows",
    "gene_locus_rows",
    "shared_eQTL_gene_locus_rows",
    "shared_sQTL_gene_locus_rows",
    "strict_shared_sQTL_event_locus_gene_rows",
    "formal_focal_tests",
    "priority_candidate_rows",
    "top50_rows"
  ),

  value = c(
    nrow(CONJ),
    nrow(FUMA_LOCI),
    nrow(FUMA_GENES),
    uniqueN(FUMA_GENES$gene),
    nrow(LAVA),
    nrow(COLOC),
    nrow(SUSIE_V3),
    nrow(SUSIE_V4),
    nrow(EQTL_P),
    nrow(SQTL_P),
    nrow(GENE_LOCUS),
    sum(
      GENE_LOCUS$eQTL_shared_gene,
      na.rm = TRUE
    ),
    sum(
      GENE_LOCUS$sQTL_shared_gene,
      na.rm = TRUE
    ),
    sum(
      GENE_LOCUS$sQTL_shared_event_n > 0,
      na.rm = TRUE
    ),
    nrow(FOCAL),
    nrow(PRIORITY_CANDIDATES),
    nrow(TOP50)
  )
)

QC_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V3_QC_summary.csv"
)

fwrite(
  QC,
  QC_FILE
)

# ============================================================
# 32. METHOD PROVENANCE / INTERPRETATION NOTES
# ============================================================

PROVENANCE <- data.table(
  item = c(
    "locus_definition",
    "gene_mapping",
    "LAVA_mapping",
    "coloc_mapping",
    "SuSiE_formal_result",
    "SuSiE_reliability_audit",
    "eQTL_primary_definition",
    "sQTL_primary_definition",
    "gene_candidate_order",
    "causal_claim_rule"
  ),
  value = c(
    "FINAL extracted FUMA GenomicRiskLoci from six standardized conjFDR pairs",
    "FUMA genes.txt joined by source trait-pair and GenomicLocus",
    "Same trait-pair plus genomic interval overlap; strict evidence uses pair-specific Bonferroni",
    "Frozen STEP8A target region plus same trait-pair",
    "STEP8C_V3 formal eight-test multi-signal SuSiE summary",
    "STEP8C_V4 SER audit appended only to the four flagged/reliability-audited tests",
    "STEP9B Primary_pass only",
    "STEP9C2_V3 full-genome Primary_pass only",
    "Lexicographic evidence structure; no arbitrary summed causal score",
    "No gene or locus is labelled causal solely from this integration"
  )
)

PROVENANCE_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V3_method_provenance.csv"
)

fwrite(
  PROVENANCE,
  PROVENANCE_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_ROOT,
    "06_LOGS",
    "STEP10B_V3_sessionInfo.txt"
  )
)

# ============================================================
# 33. FINAL CONSOLE REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP10B_V3.2 FINAL LOCUS-CENTRIC INTEGRATION COMPLETED\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nQC:\n"
)

print(
  QC
)

cat(
  "\nEight formal targeted tests:\n"
)

print(
  FOCAL[
    ,
    .(
      region,
      pair_display,
      conjFDR_FUMA_locus_present,
      LAVA_pair_bonf_sig,
      LAVA_rho_at_min_p,
      coloc_PP_H3,
      coloc_PP_H4,
      SuSiE_class,
      eQTL_shared_gene_n,
      eQTL_shared_genes,
      sQTL_shared_gene_n,
      sQTL_shared_genes,
      sQTL_shared_event_n,
      final_test_interpretation
    )
  ]
)

cat(
  "\nTop locus-gene candidates by evidence structure:\n"
)

print(
  TOP50[
    1:min(
      30L,
      .N
    ),
    .(
      Evidence_rank,
      gene,
      pair_display,
      locus_chr,
      locus_start,
      lead_snp,
      candidate_pattern,
      FUMA_positional,
      FUMA_eqtl_mapped,
      eQTL_trait1,
      eQTL_trait2,
      eQTL_shared_gene,
      sQTL_trait1,
      sQTL_trait2,
      sQTL_shared_gene,
      sQTL_shared_event_n,
      LAVA_pair_bonf_sig,
      SuSiE_class,
      min_primary_regulatory_p
    )
  ]
)

cat(
  "\nIMPORTANT INTERPRETATION RULES:\n",
  "1. Locus-level evidence is never assigned to a gene unless that gene is FUMA-mapped to that locus.\n",
  "2. SuSiE V3 is the formal 8-test result; V4 is reliability/SER annotation only.\n",
  "3. H3-dominant evidence is shown as distinct causal signals, not negative evidence.\n",
  "4. No-signal-pair is power/fine-mapping unresolved, not proof of no sharing.\n",
  "5. Shared regulatory gene requires strict Primary_pass evidence in BOTH traits inside the same locus.\n",
  "6. No arbitrary causal score is used.\n",
  sep = ""
)

cat(
  "\nUPLOAD THESE CORE OUTPUTS:\n",
  "1) ", QC_FILE, "\n",
  "1b) ", FUMA_MULTI_LOCUS_QC_FILE, "\n",
  "2) ", TABLE1_FILE, "\n",
  "3) ", TABLE2_FILE, "\n",
  "4) ", LOCUS_FILE, "\n",
  "5) ", GENE_LOCUS_FILE, "\n",
  "6) ", FIG1_SOURCE, "\n",
  "7) ", FIG1_TIFF, "\n",
  "8) ", FIG2_SOURCE, "\n",
  "9) ", FIG2_TIFF, "\n",
  sep = ""
)

cat(
  "====================================================\n"
)
