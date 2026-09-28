# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11G1_AF_SCAVENGE_FINEMAP_PREP(1).R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11G1 — AF-priority SCAVENGE fine-mapping preparation
# Project: AF–HFpEF–BMI–OSA shared genetics
#
# GOAL
#   Prepare the PRIMARY AF GWAS (GCST90624412, EUR) for
#   genome-wide locus-level fine-mapping before SCAVENGE.
#
# WHY THIS STEP
#   SCAVENGE is designed to integrate scATAC-seq with fine-mapped
#   posterior probabilities (PIP/PP). Raw GWAS p-values are not the
#   preferred input because LD can obscure causal-cell inference.
#
# WHAT THIS SCRIPT DOES
#   1) Finds the frozen primary AF GWAS: GCST90624412.tsv.gz
#   2) Detects columns robustly
#   3) Cleans autosomal biallelic variants
#   4) Finds the existing 1000G EUR hg19 PLINK reference
#   5) Determines whether AF coordinates already match hg19 or are
#      better explained by GRCh38->hg19 liftOver
#      using genome-wide-significant variants only
#   6) Defines independent AF fine-mapping loci by merging ±1 Mb
#      windows around P < 5e-8 variants in hg19 space
#   7) Excludes the extended MHC (chr6:25–34 Mb) for LD-based analysis
#   8) Writes a frozen locus manifest for STEP11G2 SuSiE fine-mapping
#
# IMPORTANT
#   - This step DOES NOT run SuSiE yet.
#   - This step DOES NOT use the S-LDSC CM result to select loci.
#   - SCAVENGE will later be run across ALL 11,986 cells, not CM only.
#   - Primary AF GWAS remains GCST90624412 (EUR).
# ============================================================

rm(list = ls())
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = max(3600, getOption("timeout"))
)

# ============================================================
# 0. USER SETTINGS
# ============================================================

DATA_ROOT <- "D:/A/data"

AF_BASENAME <- "GCST90624412.tsv.gz"

# Frozen primary threshold
GW_SIG_P <- 5e-8

# Locus definition:
# +/- 1 Mb around each genome-wide significant AF variant,
# then merge overlapping windows.
LOCUS_FLANK_BP <- 1000000L

# Extended MHC exclusion, hg19 / GRCh37
MHC_CHR <- 6L
MHC_START <- 25000000L
MHC_END <- 34000000L

# Fallback total sample size used in prior project scripts if per-SNP N absent.
AF_N_FALLBACK <- 1840341

# ============================================================
# 1. OUTPUTS
# ============================================================

ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

SCAVENGE_ROOT <- file.path(
  ROOT,
  "04_SCAVENGE"
)

OUT_DIR <- file.path(
  SCAVENGE_ROOT,
  "01_AF_FINEMAP_PREP"
)

QC_DIR <- file.path(
  OUT_DIR,
  "00_QC"
)

INPUT_DIR <- file.path(
  OUT_DIR,
  "01_STANDARDIZED_GWAS"
)

LOCI_DIR <- file.path(
  OUT_DIR,
  "02_LOCI"
)

for (d in c(
  SCAVENGE_ROOT,
  OUT_DIR,
  QC_DIR,
  INPUT_DIR,
  LOCI_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 2. PACKAGES
# ============================================================

cran_pkgs <- c(
  "data.table"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(data.table)

# rtracklayer is only required if direct hg19 matching is poor.
HAS_RTRACKLAYER <- requireNamespace(
  "rtracklayer",
  quietly = TRUE
)

HAS_GRANGES <- requireNamespace(
  "GenomicRanges",
  quietly = TRUE
) &&
  requireNamespace(
    "IRanges",
    quietly = TRUE
)

# ============================================================
# 3. HELPERS
# ============================================================

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

find_named_file <- function(
  basename_target
) {

  direct <- file.path(
    DATA_ROOT,
    basename_target
  )

  if (file.exists(direct)) {
    return(
      norm_path(
        direct,
        TRUE
      )
    )
  }

  all_files <- list.files(
    DATA_ROOT,
    recursive = TRUE,
    full.names = TRUE
  )

  bn <- basename(
    all_files
  )

  hit <- all_files[
    tolower(bn) ==
      tolower(basename_target)
  ]

  if (!length(hit)) {
    stop(
      "Could not locate ",
      basename_target,
      " under ",
      DATA_ROOT
    )
  }

  # Prefer shortest path / least-derived copy.
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

pick_col <- function(
  cols,
  candidates,
  required = TRUE,
  label = "column"
) {

  low <- tolower(
    cols
  )

  for (cand in candidates) {

    j <- which(
      low == tolower(cand)
    )

    if (length(j)) {
      return(
        cols[j[1]]
      )
    }
  }

  if (required) {
    stop(
      paste0(
        "Could not detect ",
        label,
        ". Tried: ",
        paste(
          candidates,
          collapse = ", "
        ),
        "\nAvailable columns:\n",
        paste(
          cols,
          collapse = " | "
        )
      )
    )
  }

  NA_character_
}

is_avail <- function(x) {
  length(x) == 1L &&
    !is.na(x) &&
    nzchar(x)
}

allele_pair_key <- function(a1, a2) {

  a1 <- toupper(
    as.character(a1)
  )
  a2 <- toupper(
    as.character(a2)
  )

  lo <- ifelse(
    a1 <= a2,
    a1,
    a2
  )

  hi <- ifelse(
    a1 <= a2,
    a2,
    a1
  )

  paste0(
    lo,
    "/",
    hi
  )
}

valid_base_allele <- function(x) {
  !is.na(x) &
    toupper(x) %chin% c(
      "A",
      "C",
      "G",
      "T"
    )
}

# ============================================================
# 4. LOCATE PRIMARY AF GWAS
# ============================================================

AF_FILE <- find_named_file(
  AF_BASENAME
)

cat(
  "\nPrimary AF GWAS:\n",
  AF_FILE,
  "\n",
  sep = ""
)

# ============================================================
# 5. DETECT AF COLUMNS
# ============================================================

hdr <- fread(
  AF_FILE,
  nrows = 0,
  check.names = FALSE,
  showProgress = FALSE
)

COLS <- names(
  hdr
)

COL_CHR <- pick_col(
  COLS,
  c(
    "chromosome",
    "chr",
    "CHR",
    "#chrom"
  ),
  TRUE,
  "AF chromosome"
)

COL_BP <- pick_col(
  COLS,
  c(
    "base_pair_location",
    "position",
    "pos",
    "BP",
    "bp"
  ),
  TRUE,
  "AF base-pair position"
)

COL_A1 <- pick_col(
  COLS,
  c(
    "effect_allele",
    "A1",
    "EA",
    "alt"
  ),
  TRUE,
  "AF effect allele"
)

COL_A2 <- pick_col(
  COLS,
  c(
    "other_allele",
    "A2",
    "NEA",
    "non_effect_allele",
    "ref"
  ),
  TRUE,
  "AF non-effect allele"
)

COL_BETA <- pick_col(
  COLS,
  c(
    "beta",
    "BETA",
    "effect",
    "estimate"
  ),
  FALSE,
  "AF beta"
)

COL_OR <- pick_col(
  COLS,
  c(
    "odds_ratio",
    "OR",
    "or"
  ),
  FALSE,
  "AF odds ratio"
)

if (
  !is_avail(COL_BETA) &&
  !is_avail(COL_OR)
) {
  stop(
    "AF file contains neither beta nor odds_ratio."
  )
}

COL_SE <- pick_col(
  COLS,
  c(
    "standard_error",
    "SE",
    "se",
    "stderr",
    "sebeta"
  ),
  TRUE,
  "AF standard error"
)

COL_P <- pick_col(
  COLS,
  c(
    "p_value",
    "P",
    "p",
    "pval",
    "pvalue"
  ),
  TRUE,
  "AF p value"
)

COL_EAF <- pick_col(
  COLS,
  c(
    "effect_allele_frequency",
    "EAF",
    "eaf",
    "freq",
    "af"
  ),
  FALSE,
  "AF EAF"
)

COL_N <- pick_col(
  COLS,
  c(
    "n",
    "N",
    "n_total",
    "N_total",
    "sample_size"
  ),
  FALSE,
  "AF sample size"
)

COLUMN_MAP <- data.table(
  field = c(
    "CHR",
    "BP",
    "A1",
    "A2",
    "BETA",
    "OR",
    "SE",
    "P",
    "EAF",
    "N"
  ),
  source_column = c(
    COL_CHR,
    COL_BP,
    COL_A1,
    COL_A2,
    COL_BETA,
    COL_OR,
    COL_SE,
    COL_P,
    COL_EAF,
    COL_N
  )
)

fwrite(
  COLUMN_MAP,
  file.path(
    QC_DIR,
    "STEP11G1_AF_column_map.csv"
  )
)

# ============================================================
# 6. READ REQUIRED AF FIELDS ONLY
# ============================================================

select_cols <- unique(
  na.omit(
    c(
      COL_CHR,
      COL_BP,
      COL_A1,
      COL_A2,
      COL_BETA,
      COL_OR,
      COL_SE,
      COL_P,
      COL_EAF,
      COL_N
    )
  )
)

cat(
  "\nReading AF GWAS columns...\n"
)

AF_RAW <- fread(
  AF_FILE,
  select = select_cols,
  showProgress = TRUE
)

# ============================================================
# 7. STANDARDIZE
# ============================================================

AF <- data.table(
  CHR_SOURCE = suppressWarnings(
    as.integer(
      gsub(
        "^chr",
        "",
        as.character(
          AF_RAW[[COL_CHR]]
        ),
        ignore.case = TRUE
      )
    )
  ),

  BP_SOURCE = suppressWarnings(
    as.integer(
      AF_RAW[[COL_BP]]
    )
  ),

  A1 = toupper(
    as.character(
      AF_RAW[[COL_A1]]
    )
  ),

  A2 = toupper(
    as.character(
      AF_RAW[[COL_A2]]
    )
  ),

  SE = suppressWarnings(
    as.numeric(
      AF_RAW[[COL_SE]]
    )
  ),

  P = suppressWarnings(
    as.numeric(
      AF_RAW[[COL_P]]
    )
  )
)

if (is_avail(COL_BETA)) {

  AF[
    ,
    BETA :=
      suppressWarnings(
        as.numeric(
          AF_RAW[[COL_BETA]]
        )
      )
  ]

} else {

  orv <- suppressWarnings(
    as.numeric(
      AF_RAW[[COL_OR]]
    )
  )

  AF[
    ,
    BETA :=
      log(
        orv
      )
  ]
}

if (is_avail(COL_EAF)) {

  AF[
    ,
    EAF :=
      suppressWarnings(
        as.numeric(
          AF_RAW[[COL_EAF]]
        )
      )
  ]

} else {

  AF[
    ,
    EAF :=
      NA_real_
  ]
}

if (is_avail(COL_N)) {

  AF[
    ,
    N :=
      suppressWarnings(
        as.numeric(
          AF_RAW[[COL_N]]
        )
      )
  ]

  AF[
    !is.finite(N) |
      N <= 0,
    N := AF_N_FALLBACK
  ]

} else {

  AF[
    ,
    N :=
      AF_N_FALLBACK
  ]
}

AF[
  ,
  Z :=
    BETA /
    SE
]

AF[
  ,
  allele_key :=
    allele_pair_key(
      A1,
      A2
    )
]

# ============================================================
# 8. QC FILTERS — NO RESULT-DRIVEN RELAXATION
# ============================================================

AF_QC_COUNTS <- list()

add_count <- function(stage, n) {
  AF_QC_COUNTS[[
    length(AF_QC_COUNTS) + 1L
  ]] <<- data.table(
    stage = stage,
    n = n
  )
}

add_count(
  "raw_rows",
  nrow(AF)
)

AF <- AF[
  CHR_SOURCE %between% c(1L, 22L)
]

add_count(
  "autosomes_chr1_22",
  nrow(AF)
)

AF <- AF[
  !is.na(BP_SOURCE) &
    BP_SOURCE > 0L
]

add_count(
  "valid_position",
  nrow(AF)
)

AF <- AF[
  valid_base_allele(A1) &
    valid_base_allele(A2) &
    A1 != A2
]

add_count(
  "biallelic_ACGT",
  nrow(AF)
)

AF <- AF[
  is.finite(BETA) &
    is.finite(SE) &
    SE > 0 &
    is.finite(P) &
    P > 0 &
    P <= 1 &
    is.finite(Z)
]

add_count(
  "valid_beta_se_p_z",
  nrow(AF)
)

# Deduplicate exact source-coordinate/allele-pair records:
# keep the most significant row if duplicates exist.
setorder(
  AF,
  CHR_SOURCE,
  BP_SOURCE,
  allele_key,
  P
)

AF <- AF[
  ,
  .SD[1],
  by = .(
    CHR_SOURCE,
    BP_SOURCE,
    allele_key
  )
]

add_count(
  "deduplicated",
  nrow(AF)
)

AF_QC <- rbindlist(
  AF_QC_COUNTS
)

fwrite(
  AF_QC,
  file.path(
    QC_DIR,
    "STEP11G1_AF_GWAS_QC_counts.csv"
  )
)

# ============================================================
# 9. AUTO-DETECT EXISTING 1000G EUR PLINK REFERENCE
# ============================================================

find_ref_chr1_bim <- function() {

  preferred_roots <- c(
    file.path(
      DATA_ROOT,
      "CELLULAR",
      "04_SLDSC",
      "STEP10D1",
      "02_LDSC_reference",
      "plink"
    ),
    file.path(
      DATA_ROOT,
      "STEP10_CELLULAR",
      "04_SLDSC",
      "STEP10D1",
      "02_LDSC_reference",
      "plink"
    ),
    DATA_ROOT
  )

  hits <- character(0)

  for (
    rr in preferred_roots[
      dir.exists(
        preferred_roots
      )
    ]
  ) {

    z <- list.files(
      rr,
      pattern = "^1000G[.]EUR[.]QC[.]1[.]bim$",
      recursive = TRUE,
      full.names = TRUE
    )

    hits <- unique(
      c(
        hits,
        z
      )
    )

    if (length(hits)) {
      break
    }
  }

  hits
}

REF_CHR1 <- find_ref_chr1_bim()

if (!length(REF_CHR1)) {
  stop(
    "Could not locate 1000G.EUR.QC.1.bim."
  )
}

REF_CHR1 <- REF_CHR1[
  order(
    nchar(
      norm_path(
        REF_CHR1
      )
    )
  )
][1]

REF_PREFIX <- sub(
  "1[.]bim$",
  "",
  norm_path(
    REF_CHR1,
    TRUE
  )
)

ref_complete <- all(
  vapply(
    1:22,
    function(chr) {
      all(
        file.exists(
          paste0(
            REF_PREFIX,
            chr,
            c(
              ".bed",
              ".bim",
              ".fam"
            )
          )
        )
      )
    },
    logical(1)
  )
)

if (!ref_complete) {
  stop(
    "Detected 1000G EUR prefix is not complete for chr1-22:\n",
    REF_PREFIX
  )
}

cat(
  "\n1000G EUR prefix:\n",
  REF_PREFIX,
  "\n",
  sep = ""
)

# ============================================================
# 10. READ REFERENCE BIM POSITIONS
# ============================================================

REF_LIST <- vector(
  "list",
  22L
)

for (chr in 1:22) {

  bim <- fread(
    paste0(
      REF_PREFIX,
      chr,
      ".bim"
    ),
    header = FALSE,
    col.names = c(
      "CHR",
      "SNP",
      "CM",
      "BP",
      "REF_A1",
      "REF_A2"
    ),
    showProgress = FALSE
  )

  bim[
    ,
    allele_key :=
      allele_pair_key(
        REF_A1,
        REF_A2
      )
  ]

  REF_LIST[[chr]] <- bim[
    ,
    .(
      CHR,
      SNP,
      BP,
      REF_A1,
      REF_A2,
      allele_key
    )
  ]
}

REF <- rbindlist(
  REF_LIST
)

setkey(
  REF,
  CHR,
  BP
)

# ============================================================
# 11. GWS VARIANTS FOR BUILD AUDIT
# ============================================================

GWS <- AF[
  P < GW_SIG_P
]

# exclude MHC from fine-mapping locus definition
GWS[
  ,
  in_MHC :=
    CHR_SOURCE == MHC_CHR &
    BP_SOURCE >= MHC_START &
    BP_SOURCE <= MHC_END
]

GWS_NONMHC <- GWS[
  in_MHC == FALSE
]

if (!nrow(GWS_NONMHC)) {
  stop(
    "No non-MHC AF variants reached P < ",
    GW_SIG_P
  )
}

cat(
  "\nGenome-wide significant AF variants (non-MHC): ",
  nrow(GWS_NONMHC),
  "\n",
  sep = ""
)

# ============================================================
# 12. DIRECT COORDINATE MATCH AGAINST hg19 REFERENCE
# ============================================================

GWS_DIRECT <- merge(
  GWS_NONMHC,
  REF,
  by.x = c(
    "CHR_SOURCE",
    "BP_SOURCE"
  ),
  by.y = c(
    "CHR",
    "BP"
  ),
  all.x = TRUE,
  allow.cartesian = TRUE
)

GWS_DIRECT[
  ,
  position_matched :=
    !is.na(SNP)
]

GWS_DIRECT[
  ,
  allele_matched :=
    position_matched &
    allele_key.x == allele_key.y
]

DIRECT_POS_RATE <- mean(
  GWS_DIRECT$position_matched
)

DIRECT_ALLELE_RATE <- mean(
  GWS_DIRECT$allele_matched
)

# ============================================================
# 13. OPTIONAL GRCh38 -> hg19 AUDIT
# ============================================================

CHAIN_FILE <- file.path(
  DATA_ROOT,
  "hg38ToHg19.over.chain"
)

LIFT_POS_RATE <- NA_real_
LIFT_ALLELE_RATE <- NA_real_
GWS_LIFT <- NULL

if (
  file.exists(CHAIN_FILE) &&
  HAS_RTRACKLAYER &&
  HAS_GRANGES
) {

  suppressPackageStartupMessages(
    library(GenomicRanges)
  )
  suppressPackageStartupMessages(
    library(IRanges)
  )
  suppressPackageStartupMessages(
    library(rtracklayer)
  )

  CHAIN <- rtracklayer::import.chain(
    CHAIN_FILE
  )

  gr38 <- GRanges(
    seqnames = paste0(
      "chr",
      GWS_NONMHC$CHR_SOURCE
    ),
    ranges = IRanges(
      start = GWS_NONMHC$BP_SOURCE,
      width = 1L
    )
  )

  lifted <- rtracklayer::liftOver(
    gr38,
    CHAIN
  )

  nmap <- lengths(
    lifted
  )

  unique_idx <- which(
    nmap == 1L
  )

  lift_dt <- data.table(
    row_id = unique_idx,
    CHR_HG19 = suppressWarnings(
      as.integer(
        gsub(
          "^chr",
          "",
          as.character(
            seqnames(
              unlist(
                lifted[
                  unique_idx
                ],
                use.names = FALSE
              )
            )
          )
        )
      )
    ),
    BP_HG19 = start(
      unlist(
        lifted[
          unique_idx
        ],
        use.names = FALSE
      )
    )
  )

  base_lift <- copy(
    GWS_NONMHC
  )

  base_lift[
    ,
    row_id :=
      seq_len(.N)
  ]

  GWS_LIFT <- merge(
    base_lift,
    lift_dt,
    by = "row_id",
    all.x = TRUE
  )

  GWS_LIFT <- merge(
    GWS_LIFT,
    REF,
    by.x = c(
      "CHR_HG19",
      "BP_HG19"
    ),
    by.y = c(
      "CHR",
      "BP"
    ),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

  GWS_LIFT[
    ,
    position_matched :=
      !is.na(SNP)
  ]

  GWS_LIFT[
    ,
    allele_matched :=
      position_matched &
      allele_key.x == allele_key.y
  ]

  LIFT_POS_RATE <- mean(
    GWS_LIFT$position_matched,
    na.rm = FALSE
  )

  LIFT_ALLELE_RATE <- mean(
    GWS_LIFT$allele_matched,
    na.rm = FALSE
  )
}

# ============================================================
# 14. SOURCE BUILD DECISION
# ============================================================

BUILD_AUDIT <- data.table(
  mapping = c(
    "direct_source_coordinates_as_hg19",
    "GRCh38_to_hg19_liftOver"
  ),
  position_match_rate = c(
    DIRECT_POS_RATE,
    LIFT_POS_RATE
  ),
  allele_pair_match_rate = c(
    DIRECT_ALLELE_RATE,
    LIFT_ALLELE_RATE
  )
)

fwrite(
  BUILD_AUDIT,
  file.path(
    QC_DIR,
    "STEP11G1_AF_coordinate_build_audit.csv"
  )
)

# Prefer allele-concordant mapping.
# Require a clear winner; do not guess if both are poor.
direct_score <- DIRECT_ALLELE_RATE
lift_score <- LIFT_ALLELE_RATE

if (
  is.finite(direct_score) &&
  direct_score >= 0.70 &&
  (
    !is.finite(lift_score) ||
    direct_score >= lift_score + 0.10
  )
) {

  SOURCE_BUILD_DECISION <- "SOURCE_COORDINATES_ARE_hg19"
  MAPPING_METHOD <- "DIRECT"

} else if (
  is.finite(lift_score) &&
  lift_score >= 0.70 &&
  lift_score >= direct_score + 0.10
) {

  SOURCE_BUILD_DECISION <- "SOURCE_COORDINATES_ARE_GRCh38"
  MAPPING_METHOD <- "LIFTOVER"

} else {

  SOURCE_BUILD_DECISION <- "UNRESOLVED"
  MAPPING_METHOD <- "STOP"
}

cat(
  "\nCoordinate build audit:\n"
)

print(
  BUILD_AUDIT
)

cat(
  "\nDecision: ",
  SOURCE_BUILD_DECISION,
  "\n",
  sep = ""
)

if (MAPPING_METHOD == "STOP") {
  stop(
    paste0(
      "Could not safely resolve the AF source coordinate build.\n",
      "Inspect STEP11G1_AF_coordinate_build_audit.csv before fine-mapping."
    )
  )
}

# ============================================================
# 15. ASSIGN hg19 COORDINATES TO GWS VARIANTS
# ============================================================

if (MAPPING_METHOD == "DIRECT") {

  GWS_MAP <- copy(
    GWS_NONMHC
  )

  GWS_MAP[
    ,
    `:=`(
      CHR_HG19 = CHR_SOURCE,
      BP_HG19 = BP_SOURCE
    )
  ]

} else {

  # Recreate one-row-per-source-variant unique lift map,
  # without the reference-merge duplication.
  gr38 <- GRanges(
    seqnames = paste0(
      "chr",
      GWS_NONMHC$CHR_SOURCE
    ),
    ranges = IRanges(
      start = GWS_NONMHC$BP_SOURCE,
      width = 1L
    )
  )

  lifted <- rtracklayer::liftOver(
    gr38,
    CHAIN
  )

  nmap <- lengths(
    lifted
  )

  GWS_MAP <- copy(
    GWS_NONMHC
  )

  GWS_MAP[
    ,
    `:=`(
      CHR_HG19 = NA_integer_,
      BP_HG19 = NA_integer_,
      liftover_mapping_count = nmap
    )
  ]

  ui <- which(
    nmap == 1L
  )

  ugr <- unlist(
    lifted[
      ui
    ],
    use.names = FALSE
  )

  GWS_MAP[
    ui,
    CHR_HG19 :=
      suppressWarnings(
        as.integer(
          gsub(
            "^chr",
            "",
            as.character(
              seqnames(
                ugr
              )
            )
          )
        )
      )
  ]

  GWS_MAP[
    ui,
    BP_HG19 :=
      start(
        ugr
      )
  ]

  GWS_MAP <- GWS_MAP[
    !is.na(CHR_HG19) &
      !is.na(BP_HG19) &
      CHR_HG19 %between% c(
        1L,
        22L
      )
  ]
}

# ============================================================
# 16. DEFINE MERGED +/-1 Mb AF FINE-MAPPING LOCI
# ============================================================

SEEDS <- GWS_MAP[
  ,
  .(
    CHR = CHR_HG19,
    BP = BP_HG19,
    P,
    BETA,
    SE,
    A1,
    A2,
    CHR_SOURCE,
    BP_SOURCE
  )
]

SEEDS[
  ,
  `:=`(
    start = pmax(
      1L,
      BP - LOCUS_FLANK_BP
    ),
    end = BP + LOCUS_FLANK_BP
  )
]

setorder(
  SEEDS,
  CHR,
  start,
  end
)

# Merge overlapping intervals chromosome by chromosome.
merge_windows_chr <- function(dt_chr) {

  dt_chr <- copy(
    dt_chr
  )

  setorder(
    dt_chr,
    start,
    end
  )

  out <- list()

  cur_start <- dt_chr$start[1]
  cur_end <- dt_chr$end[1]
  members <- 1L

  group_id <- integer(
    nrow(dt_chr)
  )

  gid <- 1L
  group_id[1] <- gid

  if (nrow(dt_chr) > 1L) {

    for (i in 2:nrow(dt_chr)) {

      if (
        dt_chr$start[i] <= cur_end
      ) {

        cur_end <- max(
          cur_end,
          dt_chr$end[i]
        )

      } else {

        gid <- gid + 1L
        cur_start <- dt_chr$start[i]
        cur_end <- dt_chr$end[i]
      }

      group_id[i] <- gid
    }
  }

  dt_chr[
    ,
    group_id :=
      group_id
  ]

  dt_chr
}

SEEDS_GROUPED <- rbindlist(
  lapply(
    split(
      SEEDS,
      by = "CHR",
      keep.by = TRUE
    ),
    merge_windows_chr
  ),
  use.names = TRUE
)

LOCI <- SEEDS_GROUPED[
  ,
  {
    lead_idx <- which.min(P)

    .(
      locus_start_hg19 = min(start),
      locus_end_hg19 = max(end),
      n_gws_variants = .N,
      lead_chr_hg19 = CHR[lead_idx],
      lead_bp_hg19 = BP[lead_idx],
      lead_p = P[lead_idx],
      lead_beta = BETA[lead_idx],
      lead_a1 = A1[lead_idx],
      lead_a2 = A2[lead_idx],
      lead_chr_source = CHR_SOURCE[lead_idx],
      lead_bp_source = BP_SOURCE[lead_idx],
      source_min_bp = min(BP_SOURCE),
      source_max_bp = max(BP_SOURCE)
    )
  },
  by = .(
    CHR,
    group_id
  )
]

setorder(
  LOCI,
  CHR,
  locus_start_hg19
)

LOCI[
  ,
  locus_id :=
    sprintf(
      "AF_L%03d",
      seq_len(.N)
    )
]

# Source-coordinate extraction window for STEP11G2.
# This is intentionally generous around all source GWS members in the merged locus.
LOCI[
  ,
  `:=`(
    source_extract_start = pmax(
      1L,
      source_min_bp - LOCUS_FLANK_BP
    ),
    source_extract_end = source_max_bp + LOCUS_FLANK_BP
  )
]

setcolorder(
  LOCI,
  c(
    "locus_id",
    "CHR",
    "locus_start_hg19",
    "locus_end_hg19",
    "n_gws_variants",
    "lead_chr_hg19",
    "lead_bp_hg19",
    "lead_p",
    "lead_beta",
    "lead_a1",
    "lead_a2",
    "lead_chr_source",
    "lead_bp_source",
    "source_extract_start",
    "source_extract_end",
    "group_id",
    "source_min_bp",
    "source_max_bp"
  )
)

# ============================================================
# 17. SAVE STANDARDIZED GWAS + GWS + LOCUS MANIFEST
# ============================================================

AF_STANDARD_FILE <- file.path(
  INPUT_DIR,
  "AF_GCST90624412_clean_for_finemap.tsv.gz"
)

fwrite(
  AF[
    ,
    .(
      CHR_SOURCE,
      BP_SOURCE,
      A1,
      A2,
      BETA,
      SE,
      Z,
      P,
      EAF,
      N,
      allele_key
    )
  ],
  AF_STANDARD_FILE,
  sep = "\t",
  compress = "gzip"
)

GWS_FILE <- file.path(
  LOCI_DIR,
  "STEP11G1_AF_GWS_variants_mapped_hg19.csv"
)

fwrite(
  GWS_MAP,
  GWS_FILE
)

LOCI_FILE <- file.path(
  LOCI_DIR,
  "STEP11G1_AF_finemap_locus_manifest.csv"
)

fwrite(
  LOCI,
  LOCI_FILE
)

# ============================================================
# 18. LOCUS DIRECTORY SKELETONS
# ============================================================

for (i in seq_len(nrow(LOCI))) {

  ld <- file.path(
    LOCI_DIR,
    LOCI$locus_id[i]
  )

  dir.create(
    ld,
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    c(
      paste0(
        "locus_id=",
        LOCI$locus_id[i]
      ),
      paste0(
        "chr_hg19=",
        LOCI$CHR[i]
      ),
      paste0(
        "start_hg19=",
        LOCI$locus_start_hg19[i]
      ),
      paste0(
        "end_hg19=",
        LOCI$locus_end_hg19[i]
      ),
      paste0(
        "lead_bp_hg19=",
        LOCI$lead_bp_hg19[i]
      ),
      paste0(
        "lead_p=",
        format(
          LOCI$lead_p[i],
          scientific = TRUE
        )
      ),
      paste0(
        "source_build_decision=",
        SOURCE_BUILD_DECISION
      )
    ),
    file.path(
      ld,
      "LOCUS_INFO.txt"
    )
  )
}

# ============================================================
# 19. READINESS
# ============================================================

READINESS <- data.table(
  check = c(
    "primary_AF_GCST90624412_found",
    "required_AF_columns_detected",
    "AF_GWAS_cleaned",
    "1000G_EUR_chr1_22_complete",
    "nonMHC_GWS_variants_present",
    "coordinate_build_resolved",
    "hg19_GWS_mapping_created",
    "fine_mapping_loci_created",
    "standardized_AF_GWAS_saved",
    "locus_manifest_saved"
  ),
  pass = c(
    file.exists(AF_FILE) &&
      basename(AF_FILE) == AF_BASENAME,

    all(
      !is.na(
        c(
          COL_CHR,
          COL_BP,
          COL_A1,
          COL_A2,
          COL_SE,
          COL_P
        )
      )
    ) &&
      (
        is_avail(COL_BETA) ||
        is_avail(COL_OR)
      ),

    nrow(AF) > 100000L,

    ref_complete,

    nrow(GWS_NONMHC) > 0L,

    MAPPING_METHOD %in% c(
      "DIRECT",
      "LIFTOVER"
    ),

    nrow(GWS_MAP) > 0L,

    nrow(LOCI) > 0L,

    file.exists(
      AF_STANDARD_FILE
    ),

    file.exists(
      LOCI_FILE
    )
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11G1_readiness.csv"
  )
)

# ============================================================
# 20. PROVENANCE
# ============================================================

METHOD <- data.table(
  field = c(
    "trait",
    "primary_GWAS",
    "primary_population",
    "GW_significance_threshold",
    "MHC_exclusion",
    "locus_flank",
    "locus_merging",
    "source_coordinate_decision",
    "target_finemap_build",
    "LD_reference",
    "fine_mapping_method_next",
    "SCAVENGE_scope_next"
  ),
  value = c(
    "Atrial fibrillation",
    "GCST90624412",
    "European ancestry",
    "P < 5e-8",
    "chr6:25-34 Mb, hg19",
    "+/- 1 Mb around GWS variants",
    "overlapping windows merged within chromosome",
    SOURCE_BUILD_DECISION,
    "GRCh37/hg19",
    "1000 Genomes Phase 3 EUR",
    "SuSiE RSS per frozen locus in STEP11G2",
    "all 11,986 GSE238242 ATAC cells; S-LDSC CM result used for validation, not filtering"
  )
)

fwrite(
  METHOD,
  file.path(
    QC_DIR,
    "STEP11G1_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11G1_sessionInfo.txt"
  )
)

# ============================================================
# 21. SUMMARY
# ============================================================

cat(
  "\n============================================================\n",
  "STEP11G1 COMPLETE — AF SCAVENGE FINE-MAPPING PREPARATION\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Primary AF GWAS:\n",
  AF_FILE,
  "\n\n",
  "Clean AF variants: ",
  format(
    nrow(AF),
    big.mark = ","
  ),
  "\n",
  "Non-MHC GWS variants: ",
  format(
    nrow(GWS_NONMHC),
    big.mark = ","
  ),
  "\n",
  "Fine-mapping loci: ",
  nrow(LOCI),
  "\n",
  "Coordinate decision: ",
  SOURCE_BUILD_DECISION,
  "\n\n",
  sep = ""
)

cat(
  "Build audit:\n"
)

print(
  BUILD_AUDIT
)

cat(
  "\nFirst fine-mapping loci:\n"
)

print(
  head(
    LOCI,
    20
  )
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
    "STEP11G1_readiness.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G1_AF_coordinate_build_audit.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G1_AF_GWAS_QC_counts.csv"
  ),
  "\n",
  file.path(
    QC_DIR,
    "STEP11G1_method_provenance.csv"
  ),
  "\n",
  LOCI_FILE,
  "\n",
  GWS_FILE,
  "\n",
  sep = ""
)

cat(
  "\nDo NOT start STEP11G2 until the build/locus QC is reviewed.\n"
)
