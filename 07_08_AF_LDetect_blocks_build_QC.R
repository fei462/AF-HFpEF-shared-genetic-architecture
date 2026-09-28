# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G2A_V4_AUTO_GUNZIP_CHAIN_AND_BUILD_QC(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP11G2A V4 — AF-priority SCAVENGE: robust LDetect + auto-gunzip chain + build QC
# Project: AF–HFpEF–BMI–OSA
#
# PURPOSE
#   Convert the STEP11G1 genome-wide-significant AF signals into the exact
#   fine-mapping partition used by the AF/scATAC literature:
#
#       Berisa-Pickrell EUR LDetect blocks (GRCh37/hg19)
#       -> retain blocks with >=1 AF SNP at P < 5e-8
#       -> map selected block intervals back to GRCh38
#       -> verify against original GCST90624412 GRCh38 coordinates
#
# IMPORTANT
#   - This step DOES NOT run SuSiE yet.
#   - The STEP11G1 +/-1 Mb loci are retained for audit only.
#   - STEP11G2B will fine-map the selected LDetect blocks with:
#         susie_rss()
#         L = 1
#         uniform prior
#         1000G EUR LD
#   - This mirrors the AF single-cell fine-mapping / SCAVENGE strategy.
# ==============================================================================

rm(list = ls())
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = max(3600, getOption("timeout"))
)

# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

DATA_ROOT <- "D:/A/data"

GW_SIG_P <- 5e-8

# hg19 extended MHC — keep excluded for this project.
MHC_CHR <- 6L
MHC_START <- 25000000L
MHC_END <- 34000000L

# Published Berisa-Pickrell EUR LDetect block file.
LDETECT_URL <- paste0(
  "https://bitbucket.org/nygcresearch/ldetect-data/",
  "raw/master/EUR/fourier_ls-all.bed"
)

# Reverse chain is needed because:
#   GWAS source = GRCh38
#   LDetect / 1000G fine-mapping reference = GRCh37/hg19
HG19_TO_HG38_URL <- paste0(
  "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/",
  "hg19ToHg38.over.chain.gz"
)

# ==============================================================================
# 1. INPUT / OUTPUT PATHS
# ==============================================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

G1_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "01_AF_FINEMAP_PREP"
)

G1_QC_DIR <- file.path(
  G1_ROOT,
  "00_QC"
)

G1_LOCI_DIR <- file.path(
  G1_ROOT,
  "02_LOCI"
)

G1_READINESS <- file.path(
  G1_QC_DIR,
  "STEP11G1_readiness.csv"
)

G1_BUILD_AUDIT <- file.path(
  G1_QC_DIR,
  "STEP11G1_AF_coordinate_build_audit.csv"
)

G1_GWS <- file.path(
  G1_LOCI_DIR,
  "STEP11G1_AF_GWS_variants_mapped_hg19.csv"
)

G1_OLD_MANIFEST <- file.path(
  G1_LOCI_DIR,
  "STEP11G1_AF_finemap_locus_manifest.csv"
)

OUT_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE",
  "02_AF_LDETECT_BLOCKS"
)

REF_DIR <- file.path(
  OUT_ROOT,
  "00_REFERENCE"
)

QC_DIR <- file.path(
  OUT_ROOT,
  "01_QC"
)

BLOCK_DIR <- file.path(
  OUT_ROOT,
  "02_SELECTED_BLOCKS"
)

for (d in c(
  OUT_ROOT,
  REF_DIR,
  QC_DIR,
  BLOCK_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. PACKAGES
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)

bioc_pkgs <- c(
  "GenomicRanges",
  "IRanges",
  "rtracklayer"
)

if (any(
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }

  missing_bioc <- bioc_pkgs[
    !vapply(
      bioc_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(IRanges)
)
suppressPackageStartupMessages(
  library(rtracklayer)
)

# ==============================================================================
# 3. HELPERS
# ==============================================================================

norm_path <- function(
  x,
  mustWork = FALSE
) {
  normalizePath(
    x,
    winslash = "/",
    mustWork = mustWork
  )
}

safe_download <- function(
  url,
  dest
) {

  dir.create(
    dirname(dest),
    recursive = TRUE,
    showWarnings = FALSE
  )

  cat(
    "\nDownloading:\n",
    url,
    "\n-> ",
    dest,
    "\n",
    sep = ""
  )

  status <- tryCatch(
    {
      download.file(
        url,
        destfile = dest,
        mode = "wb",
        quiet = FALSE
      )
      0L
    },
    error = function(e) {
      message(
        "Download error: ",
        conditionMessage(e)
      )
      1L
    }
  )

  isTRUE(status == 0L) &&
    file.exists(dest) &&
    file.info(dest)$size > 0
}

chr_num <- function(x) {
  suppressWarnings(
    as.integer(
      gsub(
        "^chr",
        "",
        as.character(x),
        ignore.case = TRUE
      )
    )
  )
}

# ==============================================================================
# 4. VERIFY STEP11G1
# ==============================================================================

required_g1 <- c(
  G1_READINESS,
  G1_BUILD_AUDIT,
  G1_GWS,
  G1_OLD_MANIFEST
)

if (!all(
  file.exists(required_g1)
)) {
  stop(
    paste0(
      "Missing STEP11G1 file(s):\n",
      paste(
        required_g1[
          !file.exists(required_g1)
        ],
        collapse = "\n"
      )
    )
  )
}

G1_READY <- fread(
  G1_READINESS
)

if (
  !all(
    G1_READY$pass %in% c(
      TRUE,
      "TRUE",
      1
    )
  )
) {
  print(G1_READY)

  stop(
    "STEP11G1 readiness is not fully PASS."
  )
}

BUILD_AUDIT <- fread(
  G1_BUILD_AUDIT
)

direct_row <- BUILD_AUDIT[
  mapping ==
    "direct_source_coordinates_as_hg19"
]

lift_row <- BUILD_AUDIT[
  mapping ==
    "GRCh38_to_hg19_liftOver"
]

if (
  nrow(direct_row) != 1L ||
  nrow(lift_row) != 1L
) {
  stop(
    "Unexpected STEP11G1 build-audit format."
  )
}

if (
  !is.finite(
    lift_row$allele_pair_match_rate
  ) ||
  lift_row$allele_pair_match_rate < 0.90
) {
  stop(
    "STEP11G1 GRCh38->hg19 allele concordance is <90%; do not fine-map."
  )
}

if (
  lift_row$allele_pair_match_rate <=
    direct_row$allele_pair_match_rate
) {
  stop(
    "STEP11G1 did not clearly support GRCh38 source coordinates."
  )
}

cat(
  "\nSTEP11G1 build decision confirmed:\n",
  "Direct-as-hg19 allele match = ",
  signif(
    direct_row$allele_pair_match_rate,
    5
  ),
  "\n",
  "GRCh38->hg19 allele match   = ",
  signif(
    lift_row$allele_pair_match_rate,
    5
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 5. READ STEP11G1 GWS VARIANTS
# ==============================================================================

GWS <- fread(
  G1_GWS,
  showProgress = FALSE
)

required_gws_cols <- c(
  "CHR_SOURCE",
  "BP_SOURCE",
  "A1",
  "A2",
  "SE",
  "P",
  "BETA",
  "EAF",
  "N",
  "Z",
  "CHR_HG19",
  "BP_HG19"
)

missing_cols <- setdiff(
  required_gws_cols,
  names(GWS)
)

if (length(missing_cols)) {
  stop(
    "STEP11G1 GWS file is missing columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

# Ensure project-wide primary threshold remains frozen.
GWS <- GWS[
  is.finite(P) &
    P < GW_SIG_P &
    CHR_HG19 %between% c(
      1L,
      22L
    ) &
    is.finite(BP_HG19) &
    BP_HG19 > 0
]

# Defensive MHC exclusion in target hg19 space.
GWS <- GWS[
  !(
    CHR_HG19 == MHC_CHR &
      BP_HG19 >= MHC_START &
      BP_HG19 <= MHC_END
  )
]

GWS[
  ,
  gws_id :=
    .I
]

cat(
  "\nNon-MHC AF GWS variants entering LDetect assignment: ",
  format(
    nrow(GWS),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 6. REPAIR THE REDUNDANT STEP11G1 lead_chr_hg19 COLUMN FOR AUDIT ONLY
# ==============================================================================

OLD_MANIFEST <- fread(
  G1_OLD_MANIFEST
)

if (
  all(
    c(
      "CHR",
      "lead_chr_hg19"
    ) %in%
      names(OLD_MANIFEST)
  )
) {

  OLD_MANIFEST[
    is.na(lead_chr_hg19),
    lead_chr_hg19 :=
      CHR
  ]

  if (
    "lead_p" %in%
      names(OLD_MANIFEST)
  ) {
    OLD_MANIFEST[
      ,
      lead_p_numeric_underflow_zero :=
        lead_p == 0
    ]
  }

  fwrite(
    OLD_MANIFEST,
    file.path(
      QC_DIR,
      "STEP11G2A_STEP11G1_plusminus1Mb_manifest_repaired_AUDIT_ONLY.csv"
    )
  )
}

# ==============================================================================
# 7. FIND / DOWNLOAD ORIGINAL EUR LDETECT BLOCKS
# ==============================================================================

find_existing_ldetect <- function() {

  candidates <- list.files(
    DATA_ROOT,
    pattern = "^(fourier_ls-all[.]bed|EUR_LD_blocks[.]bed)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(candidates)) {
    return(
      NA_character_
    )
  }

  # Prefer original published filename.
  original <- candidates[
    tolower(
      basename(candidates)
    ) ==
      "fourier_ls-all.bed"
  ]

  if (length(original)) {
    return(
      norm_path(
        original[
          order(
            nchar(original)
          )
        ][1],
        TRUE
      )
    )
  }

  norm_path(
    candidates[
      order(
        nchar(candidates)
      )
    ][1],
    TRUE
  )
}

LDETECT_FILE <- find_existing_ldetect()

if (
  is.na(LDETECT_FILE) ||
  !file.exists(LDETECT_FILE)
) {

  LDETECT_FILE <- file.path(
    REF_DIR,
    "fourier_ls-all.bed"
  )

  ok <- safe_download(
    LDETECT_URL,
    LDETECT_FILE
  )

  if (!ok) {
    stop(
      paste0(
        "Could not download the original Berisa-Pickrell EUR LDetect BED.\n",
        "Please manually place EUR/fourier_ls-all.bed at:\n",
        LDETECT_FILE,
        "\nThen rerun this script."
      )
    )
  }
}

LDETECT_FILE <- norm_path(
  LDETECT_FILE,
  TRUE
)

cat(
  "\nLDetect file:\n",
  LDETECT_FILE,
  "\n",
  sep = ""
)

# ==============================================================================
# 8. PARSE LDETECT BED (hg19 / GRCh37) — ROBUST V2
# ==============================================================================
#
# IMPORTANT:
# The original Berisa-Pickrell file has lines such as:
#   chr     start     stop
#   chr1    10583     1892607
#
# In the downloaded file, spaces may surround TAB characters.
# fread() automatic separator detection can therefore choose "space"
# and split "\t" into its own field, causing numeric columns to become NA.
#
# V2 avoids auto separator detection completely:
#   - read each line
#   - trim
#   - split by one-or-more whitespace characters
#   - keep exactly the first 3 fields: chr/start/stop
# ==============================================================================

ld_lines <- readLines(
  LDETECT_FILE,
  warn = FALSE
)

ld_lines <- trimws(
  ld_lines
)

ld_lines <- ld_lines[
  nzchar(ld_lines) &
    !grepl(
      "^#",
      ld_lines
    )
]

if (!length(ld_lines)) {
  stop(
    "LDetect BED contains no non-empty data lines."
  )
}

# Save the first lines exactly as read for audit.
writeLines(
  head(
    ld_lines,
    10
  ),
  file.path(
    QC_DIR,
    "STEP11G2A_V2_LDetect_raw_first10.txt"
  )
)

# Split on ANY run of whitespace: spaces, tabs, or mixed whitespace.
ld_tokens <- strsplit(
  ld_lines,
  "[[:space:]]+",
  perl = TRUE
)

token_n <- lengths(
  ld_tokens
)

TOKEN_QC <- data.table(
  line_number = seq_along(ld_lines),
  n_fields = token_n,
  raw_line = ld_lines
)

fwrite(
  TOKEN_QC[
    seq_len(
      min(
        20L,
        .N
      )
    )
  ],
  file.path(
    QC_DIR,
    "STEP11G2A_V2_LDetect_tokenization_preview.csv"
  )
)

if (any(token_n < 3L)) {

  bad <- which(
    token_n < 3L
  )

  fwrite(
    TOKEN_QC[
      bad
    ],
    file.path(
      QC_DIR,
      "STEP11G2A_V2_LDetect_bad_lines.csv"
    )
  )

  stop(
    paste0(
      "LDetect parsing found ",
      length(bad),
      " line(s) with fewer than 3 whitespace-delimited fields. ",
      "Inspect STEP11G2A_V2_LDetect_bad_lines.csv."
    )
  )
}

LD_PARSE <- rbindlist(
  lapply(
    seq_along(ld_tokens),
    function(i) {

      x <- ld_tokens[[i]]

      data.table(
        raw_chr = x[1],
        raw_start = x[2],
        raw_stop = x[3]
      )
    }
  )
)

# Remove the header by CONTENT, not by assuming it is line 1.
LD_PARSE[
  ,
  CHR := chr_num(
    trimws(
      raw_chr
    )
  )
]

LD_PARSE[
  ,
  bed_start0 := suppressWarnings(
    as.integer(
      trimws(
        raw_start
      )
    )
  )
]

LD_PARSE[
  ,
  bed_stop := suppressWarnings(
    as.integer(
      trimws(
        raw_stop
      )
    )
  )
]

# Audit parsing BEFORE filtering.
PARSE_QC <- data.table(
  metric = c(
    "raw_nonempty_lines",
    "lines_with_at_least_3_fields",
    "numeric_autosomal_rows_before_interval_filter",
    "rows_with_valid_numeric_start",
    "rows_with_valid_numeric_stop"
  ),
  value = c(
    length(ld_lines),
    sum(token_n >= 3L),
    sum(
      LD_PARSE$CHR %between% c(1L, 22L),
      na.rm = TRUE
    ),
    sum(
      is.finite(
        LD_PARSE$bed_start0
      )
    ),
    sum(
      is.finite(
        LD_PARSE$bed_stop
      )
    )
  )
)

fwrite(
  PARSE_QC,
  file.path(
    QC_DIR,
    "STEP11G2A_V2_LDetect_parse_QC.csv"
  )
)

# Header has CHR/start/stop strings -> NA after numeric conversion.
# Keep only valid autosomal numeric intervals.
LD <- LD_PARSE[
  CHR %between% c(
    1L,
    22L
  ) &
    is.finite(
      bed_start0
    ) &
    is.finite(
      bed_stop
    ) &
    bed_start0 >= 0L &
    bed_stop > bed_start0,
  .(
    CHR,
    bed_start0,
    bed_stop
  )
]

if (!nrow(LD)) {

  cat(
    "\nLDetect parse QC:\n"
  )
  print(
    PARSE_QC
  )

  cat(
    "\nFirst parsed rows:\n"
  )
  print(
    head(
      LD_PARSE,
      10
    )
  )

  stop(
    paste0(
      "LDetect V2 parser produced zero valid intervals. ",
      "Inspect STEP11G2A_V2_LDetect_parse_QC.csv and ",
      "STEP11G2A_V2_LDetect_tokenization_preview.csv."
    )
  )
}

# Remove any accidental duplicate intervals.
LD <- unique(
  LD,
  by = c(
    "CHR",
    "bed_start0",
    "bed_stop"
  )
)

# BED convention:
#   start = 0-based
#   stop  = end coordinate / half-open boundary
#
# Convert to the 1-based closed interval used for variant overlap.
LD[
  ,
  `:=`(
    start_hg19 = bed_start0 + 1L,
    end_hg19 = bed_stop
  )
]

setorder(
  LD,
  CHR,
  start_hg19,
  end_hg19
)

LD[
  ,
  block_id :=
    sprintf(
      "EUR_LD_B%04d",
      seq_len(.N)
    )
]

LD[
  ,
  width_bp :=
    end_hg19 -
    start_hg19 +
    1L
]

# Additional structural QC: chromosome order and overlap.
LD[
  ,
  previous_end :=
    shift(
      end_hg19
    ),
  by = CHR
]

LD[
  ,
  overlaps_previous :=
    !is.na(
      previous_end
    ) &
    start_hg19 <
      previous_end
]

N_OVERLAP <- sum(
  LD$overlaps_previous,
  na.rm = TRUE
)

LD[
  ,
  previous_end := NULL
]

LDETECT_QC <- data.table(
  metric = c(
    "n_blocks",
    "n_chromosomes",
    "min_width_bp",
    "median_width_bp",
    "max_width_bp",
    "n_adjacent_overlaps"
  ),
  value = c(
    nrow(LD),
    uniqueN(
      LD$CHR
    ),
    min(
      LD$width_bp
    ),
    median(
      LD$width_bp
    ),
    max(
      LD$width_bp
    ),
    N_OVERLAP
  )
)

fwrite(
  LD,
  file.path(
    REF_DIR,
    "STEP11G2A_EUR_LDetect_blocks_hg19.csv"
  )
)

fwrite(
  LDETECT_QC,
  file.path(
    QC_DIR,
    "STEP11G2A_LDetect_reference_QC.csv"
  )
)

cat(
  "\n============================================================\n",
  "LDETECT V2 PARSE CHECK\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nRaw/parse QC:\n"
)
print(
  PARSE_QC
)

cat(
  "\nLDetect structure QC:\n"
)
print(
  LDETECT_QC
)

cat(
  "\nFirst 5 parsed EUR blocks:\n"
)
print(
  LD[
    1:min(
      5L,
      .N
    ),
    .(
      block_id,
      CHR,
      bed_start0,
      bed_stop,
      start_hg19,
      end_hg19,
      width_bp
    )
  ]
)

# The published EUR file contains ~1700 blocks (normally 1703 autosomal rows).
# Keep a modest range rather than hard-coding exactly 1703, but fail on
# clearly wrong parsing.
if (
  nrow(LD) < 1650L ||
  nrow(LD) > 1750L ||
  uniqueN(
    LD$CHR
  ) != 22L ||
  !all(
    1:22 %in%
      LD$CHR
  )
) {

  stop(
    paste0(
      "Unexpected EUR LDetect reference structure after robust parsing. ",
      "Observed ",
      nrow(LD),
      " blocks across ",
      uniqueN(LD$CHR),
      " chromosomes; expected approximately 1700 autosomal blocks."
    )
  )
}

# ==============================================================================

# 9. ASSIGN AF GWS VARIANTS TO hg19 LDETECT BLOCKS
# ==============================================================================

g_gr <- GRanges(
  seqnames = paste0(
    "chr",
    GWS$CHR_HG19
  ),
  ranges = IRanges(
    start = GWS$BP_HG19,
    width = 1L
  )
)

b_gr <- GRanges(
  seqnames = paste0(
    "chr",
    LD$CHR
  ),
  ranges = IRanges(
    start = LD$start_hg19,
    end = LD$end_hg19
  )
)

ov <- findOverlaps(
  g_gr,
  b_gr,
  type = "within",
  ignore.strand = TRUE
)

ASSIGN <- data.table(
  gws_id = queryHits(
    ov
  ),
  block_index = subjectHits(
    ov
  )
)

# Each point should map to at most one approximately independent block.
dup_gws <- ASSIGN[
  ,
  .N,
  by = gws_id
][
  N > 1L
]

if (nrow(dup_gws)) {
  stop(
    "Some GWS variants overlap >1 LDetect block; inspect reference BED."
  )
}

ASSIGN[
  ,
  block_id :=
    LD$block_id[
      block_index
    ]
]

GWS_ASSIGNED <- merge(
  GWS,
  ASSIGN[
    ,
    .(
      gws_id,
      block_id,
      block_index
    )
  ],
  by = "gws_id",
  all.x = TRUE
)

ASSIGNED_RATE <- mean(
  !is.na(
    GWS_ASSIGNED$block_id
  )
)

cat(
  "\nGWS -> LDetect assignment rate: ",
  sprintf(
    "%.4f%%",
    100 * ASSIGNED_RATE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(
    ASSIGNED_RATE
  ) ||
  ASSIGNED_RATE < 0.95
) {
  stop(
    paste0(
      "Only ",
      round(
        100 * ASSIGNED_RATE,
        2
      ),
      "% of GWS variants were assigned to EUR LDetect blocks; expected >=95%."
    )
  )
}

fwrite(
  GWS_ASSIGNED,
  file.path(
    BLOCK_DIR,
    "STEP11G2A_AF_GWS_to_LDetect_assignment.csv"
  )
)

# ==============================================================================
# 10. FREEZE SELECTED AF BLOCKS
# ==============================================================================

BLOCK_GWS_SUMMARY <- GWS_ASSIGNED[
  !is.na(block_id),
  {

    o <- order(
      P,
      -abs(Z),
      BP_HG19
    )

    lead <- o[1]

    .(
      n_gws_variants = .N,
      n_p_underflow_zero = sum(
        P == 0,
        na.rm = TRUE
      ),
      lead_chr_hg19 = CHR_HG19[lead],
      lead_bp_hg19 = BP_HG19[lead],
      lead_chr_source = CHR_SOURCE[lead],
      lead_bp_source = BP_SOURCE[lead],
      lead_p = P[lead],
      lead_z = Z[lead],
      lead_beta = BETA[lead],
      lead_se = SE[lead],
      lead_a1 = A1[lead],
      lead_a2 = A2[lead]
    )
  },
  by = block_id
]

SELECTED <- merge(
  LD,
  BLOCK_GWS_SUMMARY,
  by = "block_id",
  all = FALSE,
  sort = FALSE
)

setorder(
  SELECTED,
  CHR,
  start_hg19
)

SELECTED[
  ,
  selected_order :=
    seq_len(.N)
]

SELECTED[
  ,
  selected_block_label :=
    sprintf(
      "AF_LDB%03d",
      selected_order
    )
]

setcolorder(
  SELECTED,
  c(
    "selected_block_label",
    "block_id",
    "selected_order",
    "CHR",
    "start_hg19",
    "end_hg19",
    "width_bp",
    "n_gws_variants",
    "n_p_underflow_zero",
    "lead_chr_hg19",
    "lead_bp_hg19",
    "lead_p",
    "lead_z",
    "lead_beta",
    "lead_se",
    "lead_a1",
    "lead_a2",
    "lead_chr_source",
    "lead_bp_source",
    "bed_start0",
    "bed_stop"
  )
)

SELECTED_FILE <- file.path(
  BLOCK_DIR,
  "STEP11G2A_AF_selected_LDetect_blocks_hg19.csv"
)

fwrite(
  SELECTED,
  SELECTED_FILE
)

cat(
  "\nSelected LDetect blocks containing >=1 AF GWS SNP: ",
  nrow(SELECTED),
  "\n",
  sep = ""
)

# ==============================================================================
# 11. FIND / DOWNLOAD hg19 -> hg38 CHAIN — SSL-SAFE V3
# ==============================================================================
#
# The previous run already proved that LDetect parsing and AF block selection
# reached this point. The only failure was R/libcurl SSL certificate handling
# against the UCSC server.
#
# V3 strategy:
#   A) Prefer an existing local hg19ToHg38 chain anywhere under D:/A/data
#   B) Try official UCSC URL with R download.file()
#   C) Try Windows curl.exe with --ssl-no-revoke
#   D) Try Windows PowerShell Invoke-WebRequest
#   E) If all automated methods fail, stop cleanly and tell the user the exact
#      official URL + exact destination. No upstream work is invalidated.
#
# No insecure "-k/--insecure" fallback is used.
# ==============================================================================

# Official UCSC endpoint confirmed to contain hg19ToHg38.over.chain.gz.
HG19_TO_HG38_URLS <- c(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)

CHAIN_DEST <- file.path(
  REF_DIR,
  "hg19ToHg38.over.chain.gz"
)

find_local_hg19_to_hg38 <- function() {

  # Search both compressed and uncompressed chain names.
  hit <- list.files(
    DATA_ROOT,
    pattern = "^hg19ToHg38[.]over[.]chain([.]gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  hit <- hit[
    file.exists(hit) &
      file.info(hit)$size > 10000
  ]

  if (!length(hit)) {
    return(
      NA_character_
    )
  }

  # Prefer shortest path; typically the user's canonical data copy.
  hit <- hit[
    order(
      nchar(
        norm_path(
          hit
        )
      )
    )
  ]

  norm_path(
    hit[1],
    TRUE
  )
}

prepare_chain_for_rtracklayer <- function(path) {

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size < 10000
  ) {
    return(
      list(
        ok = FALSE,
        import_path = NA_character_,
        source_path = path,
        message = "missing_or_too_small"
      )
    )
  }

  path <- norm_path(
    path,
    TRUE
  )

  # rtracklayer::import.chain() in the current environment is not
  # transparently decompressing the .gz file. Keep the original .gz
  # and create a plain .chain beside it.
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {

    plain_path <- sub(
      "\\.gz$",
      "",
      path,
      ignore.case = TRUE
    )

    need_unpack <- (
      !file.exists(plain_path) ||
      is.na(file.info(plain_path)$size) ||
      file.info(plain_path)$size < 10000
    )

    if (need_unpack) {

      cat(
        "\nDecompressing chain for rtracklayer:\n",
        path,
        "\n-> ",
        plain_path,
        "\n",
        sep = ""
      )

      con_in <- gzfile(
        path,
        open = "rb"
      )

      con_out <- file(
        plain_path,
        open = "wb"
      )

      ok_unpack <- tryCatch(
        {
          repeat {

            buf <- readBin(
              con_in,
              what = "raw",
              n = 1024 * 1024
            )

            if (!length(buf)) {
              break
            }

            writeBin(
              buf,
              con_out
            )
          }

          TRUE
        },
        error = function(e) {
          message(
            "Chain decompression failed: ",
            conditionMessage(e)
          )
          FALSE
        },
        finally = {
          try(
            close(con_in),
            silent = TRUE
          )
          try(
            close(con_out),
            silent = TRUE
          )
        }
      )

      if (
        !isTRUE(ok_unpack) ||
        !file.exists(plain_path) ||
        file.info(plain_path)$size < 10000
      ) {
        return(
          list(
            ok = FALSE,
            import_path = plain_path,
            source_path = path,
            message = "gunzip_failed"
          )
        )
      }
    }

    import_path <- norm_path(
      plain_path,
      TRUE
    )

  } else {

    import_path <- path
  }

  ok_import <- tryCatch(
    {
      tmp <- rtracklayer::import.chain(
        import_path
      )

      length(tmp) > 0
    },
    error = function(e) {
      message(
        "Chain import validation failed for ",
        import_path,
        ": ",
        conditionMessage(e)
      )
      FALSE
    }
  )

  list(
    ok = isTRUE(ok_import),
    import_path = import_path,
    source_path = path,
    message = if (isTRUE(ok_import)) "PASS" else "import_failed"
  )
}

validate_chain_candidate <- function(path) {

  prep <- prepare_chain_for_rtracklayer(
    path
  )

  isTRUE(
    prep$ok
  )
}

download_chain_R <- function(url, dest) {

  unlink(
    dest[
      file.exists(dest)
    ]
  )

  ok <- tryCatch(
    {
      suppressWarnings(
        download.file(
          url,
          destfile = dest,
          mode = "wb",
          quiet = FALSE,
          method = "libcurl"
        )
      )
      validate_chain_candidate(dest)
    },
    error = function(e) {
      message(
        "R/libcurl download failed: ",
        conditionMessage(e)
      )
      FALSE
    }
  )

  isTRUE(ok)
}

download_chain_curl_windows <- function(url, dest) {

  curl_exe <- Sys.which(
    "curl.exe"
  )

  if (!nzchar(curl_exe)) {
    message(
      "curl.exe not found; skipping Windows curl fallback."
    )
    return(FALSE)
  }

  unlink(
    dest[
      file.exists(dest)
    ]
  )

  args <- c(
    "-L",
    "--fail",
    "--show-error",
    "--silent",
    "--ssl-no-revoke",
    "--retry",
    "3",
    "--retry-delay",
    "2",
    "-o",
    dest,
    url
  )

  status <- tryCatch(
    system2(
      curl_exe,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      message(
        "Windows curl error: ",
        conditionMessage(e)
      )
      1L
    }
  )

  # system2(stdout=TRUE) may return character output with a status attribute.
  exit_status <- attr(
    status,
    "status"
  )

  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  if (!identical(
    as.integer(exit_status),
    0L
  )) {
    message(
      "Windows curl returned status ",
      exit_status
    )
    return(FALSE)
  }

  validate_chain_candidate(dest)
}

download_chain_powershell <- function(url, dest) {

  ps <- Sys.which(
    "powershell.exe"
  )

  if (!nzchar(ps)) {
    message(
      "powershell.exe not found; skipping PowerShell fallback."
    )
    return(FALSE)
  }

  unlink(
    dest[
      file.exists(dest)
    ]
  )

  # Single quotes are safe for these controlled paths/URLs.
  cmd <- paste0(
    "$ProgressPreference='SilentlyContinue'; ",
    "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; ",
    "Invoke-WebRequest -UseBasicParsing -Uri '",
    url,
    "' -OutFile '",
    gsub(
      "'",
      "''",
      norm_path(
        dest,
        FALSE
      )
    ),
    "'"
  )

  out <- tryCatch(
    system2(
      ps,
      args = c(
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-Command",
        shQuote(
          cmd
        )
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      message(
        "PowerShell download error: ",
        conditionMessage(e)
      )
      structure(
        character(0),
        status = 1L
      )
    }
  )

  exit_status <- attr(
    out,
    "status"
  )

  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  if (!identical(
    as.integer(exit_status),
    0L
  )) {
    message(
      "PowerShell returned status ",
      exit_status
    )
    return(FALSE)
  }

  validate_chain_candidate(dest)
}

# ------------------------------------------------------------------
# 11A. First preference: an already-existing local chain.
# ------------------------------------------------------------------

CHAIN_19_38 <- find_local_hg19_to_hg38()

if (
  !is.na(CHAIN_19_38) &&
  validate_chain_candidate(
    CHAIN_19_38
  )
) {

  cat(
    "\nUsing existing local hg19->hg38 chain:\n",
    CHAIN_19_38,
    "\n",
    sep = ""
  )

  CHAIN_DOWNLOAD_METHOD <- "existing_local"

} else {

  CHAIN_19_38 <- CHAIN_DEST
  CHAIN_DOWNLOAD_METHOD <- NA_character_

  cat(
    "\nNo valid local hg19ToHg38 chain found.\n",
    "Trying automated official UCSC download fallbacks...\n",
    sep = ""
  )

  # ----------------------------------------------------------------
  # 11B. Try all official URLs with R/libcurl.
  # ----------------------------------------------------------------

  for (url in HG19_TO_HG38_URLS) {

    cat(
      "\nAttempt 1 — R/libcurl:\n",
      url,
      "\n",
      sep = ""
    )

    if (
      download_chain_R(
        url,
        CHAIN_DEST
      )
    ) {
      CHAIN_DOWNLOAD_METHOD <- "R_libcurl"
      break
    }
  }

  # ----------------------------------------------------------------
  # 11C. Windows curl.exe fallback.
  # ----------------------------------------------------------------

  if (is.na(CHAIN_DOWNLOAD_METHOD)) {

    for (url in HG19_TO_HG38_URLS) {

      cat(
        "\nAttempt 2 — Windows curl.exe --ssl-no-revoke:\n",
        url,
        "\n",
        sep = ""
      )

      if (
        download_chain_curl_windows(
          url,
          CHAIN_DEST
        )
      ) {
        CHAIN_DOWNLOAD_METHOD <- "Windows_curl_ssl_no_revoke"
        break
      }
    }
  }

  # ----------------------------------------------------------------
  # 11D. PowerShell TLS1.2 fallback.
  # ----------------------------------------------------------------

  if (is.na(CHAIN_DOWNLOAD_METHOD)) {

    for (url in HG19_TO_HG38_URLS) {

      cat(
        "\nAttempt 3 — PowerShell Invoke-WebRequest:\n",
        url,
        "\n",
        sep = ""
      )

      if (
        download_chain_powershell(
          url,
          CHAIN_DEST
        )
      ) {
        CHAIN_DOWNLOAD_METHOD <- "PowerShell_InvokeWebRequest"
        break
      }
    }
  }

  # ----------------------------------------------------------------
  # 11E. Clean manual-download stop if network/certificate blocks all.
  # ----------------------------------------------------------------

  if (
    is.na(CHAIN_DOWNLOAD_METHOD) ||
    !validate_chain_candidate(
      CHAIN_DEST
    )
  ) {

    MANUAL_NOTE <- file.path(
      QC_DIR,
      "STEP11G2A_V3_MANUAL_CHAIN_DOWNLOAD_REQUIRED.txt"
    )

    writeLines(
      c(
        "STEP11G2A V3 stopped ONLY because automated UCSC chain download was blocked.",
        "",
        "The LDetect parser and upstream AF/LDetect analysis are not invalidated.",
        "",
        "Download this official UCSC file in your browser:",
        HG19_TO_HG38_URLS[1],
        "",
        "Save it EXACTLY as:",
        norm_path(
          CHAIN_DEST,
          FALSE
        ),
        "",
        "Then rerun STEP11G2A V3 from the beginning.",
        "The script will find and validate the local chain automatically."
      ),
      MANUAL_NOTE
    )

    stop(
      paste0(
        "Automated hg19ToHg38 chain download is blocked by SSL/certificate settings.\n",
        "Nothing upstream needs to be rerun manually.\n\n",
        "Please download the official UCSC file in your browser:\n",
        HG19_TO_HG38_URLS[1],
        "\n\nSave it EXACTLY as:\n",
        norm_path(
          CHAIN_DEST,
          FALSE
        ),
        "\n\nThen rerun this V3 script.\n",
        "Instruction file saved at:\n",
        MANUAL_NOTE
      )
    )
  }

  CHAIN_19_38 <- norm_path(
    CHAIN_DEST,
    TRUE
  )
}

# Final definitive validation + transparent decompression.
CHAIN_PREP <- prepare_chain_for_rtracklayer(
  CHAIN_19_38
)

if (!isTRUE(CHAIN_PREP$ok)) {
  stop(
    paste0(
      "hg19ToHg38 chain exists, but preparation/import failed.\n",
      "Source file: ",
      CHAIN_19_38,
      "\nPrepared import path: ",
      CHAIN_PREP$import_path,
      "\nStatus: ",
      CHAIN_PREP$message
    )
  )
}

CHAIN_IMPORT_PATH <- CHAIN_PREP$import_path

CHAIN_QC <- data.table(
  field = c(
    "chain_source_file",
    "chain_import_file",
    "source_was_gz",
    "download_method",
    "source_file_size_bytes",
    "import_file_size_bytes",
    "rtracklayer_import_pass"
  ),
  value = c(
    CHAIN_19_38,
    CHAIN_IMPORT_PATH,
    as.character(
      grepl(
        "\\.gz$",
        CHAIN_19_38,
        ignore.case = TRUE
      )
    ),
    CHAIN_DOWNLOAD_METHOD,
    as.character(
      file.info(
        CHAIN_19_38
      )$size
    ),
    as.character(
      file.info(
        CHAIN_IMPORT_PATH
      )$size
    ),
    "TRUE"
  )
)

fwrite(
  CHAIN_QC,
  file.path(
    QC_DIR,
    "STEP11G2A_V4_chain_QC.csv"
  )
)

cat(
  "\n============================================================\n",
  "hg19 -> hg38 CHAIN READY\n",
  "============================================================\n",
  "Downloaded/source file: ",
  CHAIN_19_38,
  "\n",
  "rtracklayer import file: ",
  CHAIN_IMPORT_PATH,
  "\n",
  "Method: ",
  CHAIN_DOWNLOAD_METHOD,
  "\n",
  "Compressed size: ",
  format(
    file.info(
      CHAIN_19_38
    )$size,
    big.mark = ","
  ),
  " bytes\n",
  "Uncompressed size: ",
  format(
    file.info(
      CHAIN_IMPORT_PATH
    )$size,
    big.mark = ","
  ),
  " bytes\n",
  sep = ""
)

CHAIN <- rtracklayer::import.chain(
  CHAIN_IMPORT_PATH
)

# ==============================================================================

# 12. MAP SELECTED LDETECT BLOCKS hg19 -> GRCh38
# ==============================================================================

sel_gr19 <- GRanges(
  seqnames = paste0(
    "chr",
    SELECTED$CHR
  ),
  ranges = IRanges(
    start = SELECTED$start_hg19,
    end = SELECTED$end_hg19
  )
)

names(
  sel_gr19
) <- SELECTED$block_id

lifted_blocks <- rtracklayer::liftOver(
  sel_gr19,
  CHAIN
)

BLOCK_SOURCE_SEGMENTS_LIST <- vector(
  "list",
  length(
    lifted_blocks
  )
)

for (
  i in seq_along(
    lifted_blocks
  )
) {

  gr_i <- lifted_blocks[[i]]

  if (!length(gr_i)) {
    BLOCK_SOURCE_SEGMENTS_LIST[[i]] <- NULL
    next
  }

  chr_i <- chr_num(
    seqnames(
      gr_i
    )
  )

  keep_i <- which(
    chr_i %between%
      c(
        1L,
        22L
      ) &
      chr_i ==
        SELECTED$CHR[i]
  )

  if (!length(keep_i)) {
    BLOCK_SOURCE_SEGMENTS_LIST[[i]] <- NULL
    next
  }

  gr_i <- gr_i[
    keep_i
  ]

  BLOCK_SOURCE_SEGMENTS_LIST[[i]] <- data.table(
    block_id = SELECTED$block_id[i],
    selected_block_label =
      SELECTED$selected_block_label[i],
    CHR_SOURCE = chr_num(
      seqnames(
        gr_i
      )
    ),
    source_start = start(
      gr_i
    ),
    source_end = end(
      gr_i
    ),
    source_width = width(
      gr_i
    ),
    source_segment_index =
      seq_along(
        gr_i
      )
  )
}

BLOCK_SOURCE_SEGMENTS <- rbindlist(
  BLOCK_SOURCE_SEGMENTS_LIST,
  fill = TRUE
)

if (!nrow(BLOCK_SOURCE_SEGMENTS)) {
  stop(
    "No selected LDetect block mapped from hg19 to GRCh38."
  )
}

fwrite(
  BLOCK_SOURCE_SEGMENTS,
  file.path(
    BLOCK_DIR,
    "STEP11G2A_selected_LDetect_blocks_GRCh38_source_segments.csv"
  )
)

# ==============================================================================
# 13. PER-BLOCK REVERSE-LIFTOVER QC
# ==============================================================================

BLOCK_MAP_QC <- BLOCK_SOURCE_SEGMENTS[
  ,
  .(
    n_source_segments = .N,
    total_source_width =
      sum(
        source_width
      ),
    source_chr_n =
      uniqueN(
        CHR_SOURCE
      ),
    source_chr =
      paste(
        sort(
          unique(
            CHR_SOURCE
          )
        ),
        collapse = ","
      )
  ),
  by = .(
    block_id,
    selected_block_label
  )
]

BLOCK_MAP_QC <- merge(
  SELECTED[
    ,
    .(
      block_id,
      selected_block_label,
      CHR,
      start_hg19,
      end_hg19,
      width_bp,
      n_gws_variants
    )
  ],
  BLOCK_MAP_QC,
  by = c(
    "block_id",
    "selected_block_label"
  ),
  all.x = TRUE
)

BLOCK_MAP_QC[
  ,
  mapped_to_source :=
    !is.na(
      n_source_segments
    ) &
    n_source_segments >= 1L &
    source_chr_n == 1L
]

fwrite(
  BLOCK_MAP_QC,
  file.path(
    QC_DIR,
    "STEP11G2A_selected_block_reverse_liftover_QC.csv"
  )
)

BLOCK_MAP_RATE <- mean(
  BLOCK_MAP_QC$mapped_to_source
)

if (
  !is.finite(
    BLOCK_MAP_RATE
  ) ||
  BLOCK_MAP_RATE < 0.95
) {
  stop(
    "Fewer than 95% of selected LDetect blocks mapped cleanly to GRCh38."
  )
}

# ==============================================================================
# 14. ROUND-TRIP VALIDATION USING ACTUAL GWS SOURCE COORDINATES
# ==============================================================================

# We already know for each GWS SNP:
#   hg19 coordinate -> selected LDetect block
#   original GRCh38 source coordinate
#
# The original source coordinate should fall inside at least one
# reverse-lifted GRCh38 segment of that same LDetect block.

PTS <- GWS_ASSIGNED[
  !is.na(block_id),
  .(
    gws_id,
    block_id,
    CHR_SOURCE,
    source_start = BP_SOURCE,
    source_end = BP_SOURCE
  )
]

SEGS <- BLOCK_SOURCE_SEGMENTS[
  ,
  .(
    block_id,
    CHR_SOURCE,
    source_start,
    source_end
  )
]

setkey(
  SEGS,
  block_id,
  CHR_SOURCE,
  source_start,
  source_end
)

ROUNDTRIP_HIT <- foverlaps(
  x = PTS,
  y = SEGS,
  by.x = c(
    "block_id",
    "CHR_SOURCE",
    "source_start",
    "source_end"
  ),
  by.y = c(
    "block_id",
    "CHR_SOURCE",
    "source_start",
    "source_end"
  ),
  type = "within",
  nomatch = 0L
)

ROUNDTRIP_MATCHED_IDS <- unique(
  ROUNDTRIP_HIT$gws_id
)

ROUNDTRIP_RATE <- length(
  ROUNDTRIP_MATCHED_IDS
) /
  nrow(
    PTS
  )

ROUNDTRIP_QC <- data.table(
  metric = c(
    "n_GWS_assigned_to_selected_blocks",
    "n_GWS_inside_reverse_lifted_source_segments",
    "roundtrip_source_coordinate_rate",
    "n_selected_blocks",
    "selected_block_reverse_liftover_rate"
  ),
  value = c(
    nrow(PTS),
    length(
      ROUNDTRIP_MATCHED_IDS
    ),
    ROUNDTRIP_RATE,
    nrow(
      SELECTED
    ),
    BLOCK_MAP_RATE
  )
)

fwrite(
  ROUNDTRIP_QC,
  file.path(
    QC_DIR,
    "STEP11G2A_roundtrip_coordinate_QC.csv"
  )
)

cat(
  "\nRound-trip GWS coordinate validation: ",
  sprintf(
    "%.4f%%",
    100 * ROUNDTRIP_RATE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(
    ROUNDTRIP_RATE
  ) ||
  ROUNDTRIP_RATE < 0.90
) {
  stop(
    paste0(
      "Round-trip source-coordinate QC is <90% (",
      round(
        100 * ROUNDTRIP_RATE,
        2
      ),
      "%). Do not start SuSiE."
    )
  )
}

# ==============================================================================
# 15. BLOCK-LEVEL SUMMARY
# ==============================================================================

BLOCK_SUMMARY <- SELECTED[
  ,
  .(
    n_selected_blocks = .N,
    median_block_width_mb =
      median(
        width_bp
      ) /
      1e6,
    p90_block_width_mb =
      as.numeric(
        quantile(
          width_bp,
          0.90
        )
      ) /
      1e6,
    max_block_width_mb =
      max(
        width_bp
      ) /
      1e6,
    median_GWS_per_block =
      median(
        n_gws_variants
      ),
    max_GWS_per_block =
      max(
        n_gws_variants
      ),
    n_blocks_with_P_underflow =
      sum(
        n_p_underflow_zero > 0
      )
  )
]

fwrite(
  BLOCK_SUMMARY,
  file.path(
    QC_DIR,
    "STEP11G2A_selected_block_summary.csv"
  )
)

# ==============================================================================
# 16. READINESS
# ==============================================================================

READINESS <- data.table(
  check = c(
    "STEP11G1_all_pass",
    "STEP11G1_GRCh38_source_build_confirmed",
    "LDetect_reference_loaded",
    "LDetect_reference_approximately_1700_blocks",
    "LDetect_chr1_22_present",
    "GWS_assignment_rate_ge_95pct",
    "selected_AF_blocks_created",
    "hg19ToHg38_chain_available",
    "selected_block_reverse_liftover_rate_ge_95pct",
    "GWS_roundtrip_source_coordinate_rate_ge_90pct",
    "STEP11G2B_ready"
  ),
  pass = c(
    all(
      G1_READY$pass %in%
        c(
          TRUE,
          "TRUE",
          1
        )
    ),
    lift_row$allele_pair_match_rate >
      direct_row$allele_pair_match_rate &&
      lift_row$allele_pair_match_rate >=
      0.90,
    file.exists(
      LDETECT_FILE
    ),
    nrow(LD) >= 1600L &&
      nrow(LD) <= 1800L,
    all(
      1:22 %in%
        LD$CHR
    ),
    ASSIGNED_RATE >= 0.95,
    nrow(
      SELECTED
    ) > 0L,
    file.exists(
      CHAIN_19_38
    ),
    BLOCK_MAP_RATE >= 0.95,
    ROUNDTRIP_RATE >= 0.90,
    TRUE
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G2A_readiness.csv"
  )
)

# ==============================================================================
# 17. METHOD PROVENANCE
# ==============================================================================

METHOD <- data.table(
  field = c(
    "primary_AF_GWAS",
    "source_GWAS_build",
    "fine_mapping_target_build",
    "fine_mapping_partition",
    "LDetect_population",
    "block_selection_rule",
    "MHC_rule",
    "next_model",
    "next_model_L",
    "next_prior",
    "next_LD_reference",
    "SCAVENGE_input_plan"
  ),
  value = c(
    "GCST90624412 EUR",
    "GRCh38",
    "GRCh37/hg19",
    "Berisa-Pickrell approximately independent LD blocks",
    "EUR",
    "retain LDetect block if it contains >=1 AF SNP with P < 5e-8",
    "exclude hg19 chr6:25-34 Mb",
    "susieR::susie_rss",
    "L=1",
    "uniform prior",
    "1000 Genomes Phase 3 EUR",
    "uniform-prior fine-mapped PIP -> gchromVAR -> SCAVENGE across all GSE238242 ATAC cells"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G2A_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G2A_sessionInfo.txt"
  )
)

# ==============================================================================
# 18. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP11G2A COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "EUR LDetect blocks total: ",
  nrow(LD),
  "\n",
  "AF-selected blocks:       ",
  nrow(SELECTED),
  "\n",
  "GWS assignment rate:      ",
  sprintf(
    "%.2f%%",
    100 * ASSIGNED_RATE
  ),
  "\n",
  "Block hg19->hg38 rate:    ",
  sprintf(
    "%.2f%%",
    100 * BLOCK_MAP_RATE
  ),
  "\n",
  "GWS round-trip rate:      ",
  sprintf(
    "%.2f%%",
    100 * ROUNDTRIP_RATE
  ),
  "\n\n",
  sep = ""
)

cat(
  "Selected block summary:\n"
)

print(
  BLOCK_SUMMARY
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

cat(
  "\nUPLOAD THESE FILES:\n",
  file.path(
    QC_DIR,
    "STEP11G2A_readiness.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G2A_LDetect_reference_QC.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G2A_roundtrip_coordinate_QC.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G2A_selected_block_reverse_liftover_QC.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G2A_selected_block_summary.csv"
  ),
  "\n",
  SELECTED_FILE,
  "\n",
  file.path(
    BLOCK_DIR,
    "STEP11G2A_selected_LDetect_blocks_GRCh38_source_segments.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nIf all readiness checks are TRUE, proceed to STEP11G2B SuSiE RSS fine-mapping.\n"
)
