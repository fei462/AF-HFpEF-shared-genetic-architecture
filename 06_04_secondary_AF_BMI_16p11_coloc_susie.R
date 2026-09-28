# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP10B_V4_SECONDARY_AF_BMI_16P11_COLOC_SUSIE_FIXED_V9.R
# See repository README.md for execution order and external dependencies.


# ============================================================
# STEP10B_V4_SECONDARY_AF_BMI_16P11_COLOC_SUSIE_FIXED_V9
# Secondary/post hoc locus-level follow-up — executable FIXED V9 (SuSiE summary parser repaired; bounded nonconvergence retained)
#
# Locus:
#   chr16:28,538,336-29,008,079 (GRCh37 / hg19)
#   AF-BMI shared locus
#   FUMA lead SNP: rs12448482
#
# RATIONALE
# ---------
# This is NOT one of the eight originally pre-specified STEP8
# targeted tests. It is a SECONDARY, hypothesis-generating
# follow-up motivated by convergence of:
#   - conjFDR
#   - significant LAVA local rg
#   - FUMA mapping
#   - eQTL/sQTL regulatory evidence
#
# ANALYSES
# --------
# 1) Dense regional SNP extraction from final GWAS .ma files.
# 2) Strict allele harmonization to the LD reference.
# 3) PLINK-based signed LD from the same European reference family.
# 4) coloc.abf:
#      p12 = 1e-5 (default formal)
#      p12 = 1e-6 (conservative)
# 5) multi-signal SuSiE / coloc.susie:
#      R_finite = actual LD reference N
#      R_mismatch = "eb"
#      p12 = 1e-5 and 1e-6
# 6) LD mismatch diagnostics.
# 7) SER one-effect sensitivity audit for BOTH traits.
# 8) Standalone publication-quality figures + source data.
#
# IMPORTANT
# ---------
# - This script does NOT alter any previous STEP8 result.
# - This secondary analysis must be reported separately from the
#   eight pre-specified tests.
# - No threshold is relaxed to obtain a positive result.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "coloc",
  "susieR",
  "ggplot2"
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
library(coloc)
library(susieR)
library(ggplot2)

# Hard compatibility gate: this secondary follow-up MUST use the same
# mismatch-aware SuSiE-RSS framework as frozen STEP8C-V3.
susie_args_available <- names(formals(susieR::susie_rss))
if (!all(c("R_finite", "R_mismatch") %in% susie_args_available)) {
  stop(
    paste0(
      "Installed susieR is too old for the frozen mismatch-aware method.\n",
      "Required susie_rss() arguments: R_finite and R_mismatch.\n",
      "Please update susieR, restart R, and rerun this script."
    )
  )
}

cat(
  "\nPackage versions:\n",
  "  coloc  = ", as.character(packageVersion("coloc")), "\n",
  "  susieR = ", as.character(packageVersion("susieR")), "\n",
  sep = ""
)

# ============================================================
# 1. FROZEN PATHS AND PARAMETERS
# ============================================================

DATA_ROOT <- "D:/A/data"

GWAS_MA_DIR <- file.path(
  DATA_ROOT,
  "STEP9_SMR",
  "02_gwas_ma"
)

AF_MA <- file.path(
  GWAS_MA_DIR,
  "AF.ma"
)

BMI_MA <- file.path(
  GWAS_MA_DIR,
  "BMI.ma"
)

# Optional manual override. Leave NA to auto-detect.
# If auto-detection still fails, replace NA_character_ with the exact path
# to your PLINK 1.9 executable, e.g. "D:/A/data/.../plink.exe".
PLINK_EXE_MANUAL <- NA_character_

STEP8_ROOT <- file.path(
  DATA_ROOT,
  "STEP8_COLOC"
)

# Optional manual override for the 1000G EUR PLINK reference PREFIX.
# Leave NA for automatic discovery. If needed, give prefix only:
# e.g. "D:/.../1000G_EUR_autosomes" (without .bed/.bim/.fam).
PLINK_REF_PREFIX_MANUAL <- NA_character_

# Historical fallback from an earlier STEP8 layout.
# IMPORTANT: V4 uses it only when the BED/BIM/FAM triplet truly exists.
REF_FALLBACK <- "D:/A/data/STEP8_COLOC/00_LD_REFERENCE/g1000_eur/g1000_eur"

OUT_ROOT <- file.path(
  DATA_ROOT,
  "STEP10B",
  "STEP10B_V4_SECONDARY_AF_BMI_16P11_COLOC_SUSIE"
)

OUT_DIRS <- c(
  "00_QC",
  "01_SOURCE_DATA",
  "02_MAIN_TABLES",
  "03_SUPPLEMENTARY_TABLES",
  "04_MAIN_FIGURES",
  "05_SUPPLEMENTARY_FIGURES",
  "06_LD_WORK",
  "07_LOGS"
)

for (d in OUT_DIRS) {
  dir.create(
    file.path(
      OUT_ROOT,
      d
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ----------------------------
# Frozen secondary locus
# ----------------------------

CHR <- 16L
REGION_START <- 28538336L
REGION_STOP  <- 29008079L
FUMA_LEAD_SNP <- "rs12448482"

# ----------------------------
# AF European GWAS
# GCST90624412
# ----------------------------

AF_CASES <- 228926
AF_CONTROLS <- 1611415
AF_N_TOTAL <- AF_CASES + AF_CONTROLS
AF_CASE_FRACTION <- AF_CASES / AF_N_TOTAL

# ----------------------------
# coloc priors
# ----------------------------

P1 <- 1e-4
P2 <- 1e-4
P12_DEFAULT <- 1e-5
P12_CONSERVATIVE <- 1e-6

# ----------------------------
# LD / harmonization QC
# ----------------------------

MIN_COMMON_SNPS <- 100L
REF_MAF_MIN <- 0.01
PLINK_GENO_MAX <- 0.05

# SuSiE
SUSIE_L <- 10L
SUSIE_COVERAGE <- 0.95

# ============================================================
# 2. BASIC HELPERS
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

prefix_ok <- function(prefix) {

  all(
    file.exists(
      paste0(
        prefix,
        c(
          ".bed",
          ".bim",
          ".fam"
        )
      )
    )
  )
}

# PLINK --bfile expects a PREFIX, not a literal file.
# Validate the .bed/.bim/.fam triplet, then normalize the bare prefix
# with mustWork=FALSE. Using mustWork=TRUE on the prefix itself is wrong
# because no extensionless file is expected to exist.
plink_bfile_prefix <- function(prefix) {

  if (!prefix_ok(prefix)) {
    stop(
      paste0(
        "Invalid PLINK bfile prefix: ",
        prefix,
        "\nExpected all three files:\n",
        prefix, ".bed\n",
        prefix, ".bim\n",
        prefix, ".fam"
      )
    )
  }

  norm_path(
    prefix,
    FALSE
  )
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

safe_int <- function(x) {
  suppressWarnings(
    as.integer(
      as.character(x)
    )
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

complement_allele <- function(a) {

  a <- toupper(
    as.character(a)
  )

  out <- rep(
    NA_character_,
    length(a)
  )

  out[a == "A"] <- "T"
  out[a == "T"] <- "A"
  out[a == "C"] <- "G"
  out[a == "G"] <- "C"

  out
}

is_palindromic <- function(a1, a2) {

  x <- paste0(
    toupper(a1),
    toupper(a2)
  )

  x %in% c(
    "AT",
    "TA",
    "CG",
    "GC"
  )
}

# ============================================================
# 3. LOCATE PLINK 1.9
# ============================================================

find_plink <- function() {

  # 0) Manual override if user provided one.
  if (
    !is.na(PLINK_EXE_MANUAL) &&
    nzchar(PLINK_EXE_MANUAL) &&
    file.exists(PLINK_EXE_MANUAL)
  ) {
    hits <- PLINK_EXE_MANUAL
  } else {

    # 1) PATH + broad project roots.
    hits <- c(
      unname(Sys.which("plink")),
      unname(Sys.which("plink.exe"))
    )

    roots <- unique(
      c(
        STEP8_ROOT,
        DATA_ROOT,
        "D:/A",
        "E:/app"
      )
    )

    for (rr in roots[dir.exists(roots)]) {
      h <- tryCatch(
        list.files(
          rr,
          pattern = "^plink(1\\.9|_1\\.9)?(\\.exe)?$",
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        ),
        error = function(e) character(0)
      )

      hits <- unique(
        c(
          hits,
          h
        )
      )
    }
  }

  hits <- unique(
    hits[
      nzchar(hits) &
      file.exists(hits)
    ]
  )

  if (!length(hits)) {
    stop(
      paste0(
        "PLINK 1.9 executable was not found automatically.\n",
        "Searched PATH and:\n",
        paste(
          unique(
            c(
              STEP8_ROOT,
              DATA_ROOT,
              "D:/A",
              "E:/app"
            )
          ),
          collapse = "\n"
        ),
        "\n\nIf you know the exact location, set at the top of this script:\n",
        "PLINK_EXE_MANUAL <- 'D:/.../plink.exe'\n",
        "Do not substitute PLINK 2.x."
      )
    )
  }

  # Prefer literal plink.exe / shortest path, but verify version.
  exact <- hits[
    tolower(
      basename(hits)
    ) == "plink.exe"
  ]

  if (length(exact)) {
    hits <- c(
      exact,
      setdiff(
        hits,
        exact
      )
    )
  }

  hits <- hits[
    order(
      nchar(hits)
    )
  ]

  version_audit <- list()

  for (h in hits) {

    z <- tryCatch(
      system2(
        h,
        args = "--version",
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) conditionMessage(e)
    )

    version_audit[[
      length(version_audit) + 1L
    ]] <- data.table(
      candidate = norm_path(
        h,
        FALSE
      ),
      version_text = paste(
        z,
        collapse = " "
      )
    )

    if (
      length(z) &&
      any(
        grepl(
          "PLINK",
          z,
          ignore.case = TRUE
        )
      ) &&
      !any(
        grepl(
          "PLINK v2|PLINK 2",
          z,
          ignore.case = TRUE
        )
      )
    ) {
      audit <- rbindlist(
        version_audit,
        fill = TRUE
      )

      fwrite(
        audit,
        file.path(
          OUT_ROOT,
          "00_QC",
          "STEP10B_V4_PLINK_candidate_audit.csv"
        )
      )

      return(
        norm_path(
          h,
          TRUE
        )
      )
    }
  }

  audit <- rbindlist(
    version_audit,
    fill = TRUE
  )

  fwrite(
    audit,
    file.path(
      OUT_ROOT,
      "00_QC",
      "STEP10B_V4_PLINK_candidate_audit.csv"
    )
  )

  stop(
    paste0(
      "PLINK candidates were found, but no valid PLINK 1.9 executable passed the version check.\n",
      "See: ",
      file.path(
        OUT_ROOT,
        "00_QC",
        "STEP10B_V4_PLINK_candidate_audit.csv"
      )
    )
  )
}

PLINK_EXE <- find_plink()

cat(
  "\nPLINK 1.9 executable:\n",
  PLINK_EXE,
  "\n",
  sep = ""
)

# ============================================================
# 4. LOCATE STEP8-CONSISTENT EUR LD REFERENCE
# ============================================================
#
# V4 discovery priority:
#   1) manual prefix, if supplied;
#   2) the actual STEP8C-V3 reference manifest, if present;
#   3) known project-level merged reference from STEP9A;
#   4) valid BED/BIM/FAM triplets recursively found under
#      STEP8_COLOC / STEP9_SMR / STEP3_MR / D:/A/data.
#
# A candidate is accepted only if:
#   - BED/BIM/FAM all exist;
#   - reference N is 450-550 (covers the frozen 489-person STEP8
#     panel as well as conventional ~503-person 1000G EUR panels);
#   - it contains at least MIN_COMMON_SNPS variants in the exact
#     secondary chr16 region.
# ============================================================

count_lines_fast_local <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }

  length(
    readLines(
      path,
      warn = FALSE
    )
  )
}

reference_candidate_qc <- function(prefix, source_label) {

  if (!prefix_ok(prefix)) {
    return(NULL)
  }

  fam_file <- paste0(
    prefix,
    ".fam"
  )

  bim_file <- paste0(
    prefix,
    ".bim"
  )

  nref <- tryCatch(
    count_lines_fast_local(
      fam_file
    ),
    error = function(e) NA_integer_
  )

  if (
    !is.finite(nref) ||
    nref < 450L ||
    nref > 550L
  ) {
    return(NULL)
  }

  bim <- tryCatch(
    fread(
      bim_file,
      header = FALSE,
      select = c(
        1,
        2,
        4
      ),
      col.names = c(
        "CHR_RAW",
        "SNP",
        "BP"
      ),
      showProgress = FALSE
    ),
    error = function(e) NULL
  )

  if (
    is.null(bim) ||
    !nrow(bim)
  ) {
    return(NULL)
  }

  bim[
    ,
    CHR_NUM :=
      suppressWarnings(
        as.integer(
          CHR_RAW
        )
      )
  ]

  target_chr <- as.integer(
    CHR
  )

  target_n <- bim[
    CHR_NUM == target_chr &
      BP >= REGION_START &
      BP <= REGION_STOP,
    .N
  ]

  if (
    !is.finite(target_n) ||
    target_n < MIN_COMMON_SNPS
  ) {
    return(NULL)
  }

  n_chr <- uniqueN(
    bim[
      CHR_NUM >= 1L &
        CHR_NUM <= 22L,
      CHR_NUM
    ]
  )

  # Prefer historically used STEP8 reference, then STEP9 merged
  # reference, then clearly named 1000G EUR candidates.
  path_score <-
    400L * grepl(
      "STEP8_COLOC",
      prefix,
      ignore.case = TRUE
    ) +
    300L * grepl(
      "STEP9_SMR",
      prefix,
      ignore.case = TRUE
    ) +
    200L * grepl(
      "1000|g1000|1kg|1kgp",
      prefix,
      ignore.case = TRUE
    ) +
    150L * grepl(
      "eur|europe",
      prefix,
      ignore.case = TRUE
    ) +
    100L * grepl(
      "merged|autosome",
      prefix,
      ignore.case = TRUE
    ) +
    75L * as.integer(
      nref == 489L
    ) +
    50L * as.integer(
      n_chr >= 20L
    )

  data.table(
    prefix = norm_path(
      prefix,
      FALSE
    ),
    source = source_label,
    n_reference = nref,
    n_variants = nrow(
      bim
    ),
    n_chr = n_chr,
    n_target_region_variants = target_n,
    path_score = path_score
  )
}

find_reference <- function() {

  candidate_rows <- list()

  add_candidate <- function(prefix, label) {

    q <- reference_candidate_qc(
      prefix,
      label
    )

    if (!is.null(q)) {
      candidate_rows[[
        length(
          candidate_rows
        ) + 1L
      ]] <<- q
    }

    invisible(
      NULL
    )
  }

  # ----------------------------------------------------------
  # 4A. Manual override
  # ----------------------------------------------------------
  if (
    !is.na(
      PLINK_REF_PREFIX_MANUAL
    ) &&
    nzchar(
      PLINK_REF_PREFIX_MANUAL
    )
  ) {
    add_candidate(
      PLINK_REF_PREFIX_MANUAL,
      "manual_override"
    )
  }

  # ----------------------------------------------------------
  # 4B. Reuse the actual STEP8C-V3 reference manifest
  # ----------------------------------------------------------
  manifest_files <- tryCatch(
    list.files(
      STEP8_ROOT,
      pattern = "V3.*region_reference_bfiles.*\\.csv$|region_reference_bfiles.*\\.csv$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    ),
    error = function(e) character(0)
  )

  if (length(manifest_files)) {

    for (mf in manifest_files) {

      tab <- tryCatch(
        fread(
          mf,
          showProgress = FALSE
        ),
        error = function(e) NULL
      )

      if (
        !is.null(tab) &&
        "prefix" %in%
          names(tab)
      ) {

        manifest_prefixes <- unique(
          as.character(
            tab$prefix
          )
        )

        for (pr in manifest_prefixes) {
          add_candidate(
            pr,
            paste0(
              "STEP8C_manifest:",
              basename(
                mf
              )
            )
          )
        }
      }
    }
  }

  # ----------------------------------------------------------
  # 4C. Known project-level reference prefixes
  # ----------------------------------------------------------
  known_prefixes <- c(
    REF_FALLBACK,
    file.path(
      DATA_ROOT,
      "STEP9_SMR",
      "00_reference",
      "1000G_EUR_merged",
      "1000G_EUR_autosomes"
    )
  )

  for (pr in known_prefixes) {
    add_candidate(
      pr,
      "known_project_reference"
    )
  }

  # ----------------------------------------------------------
  # 4D. Recursive discovery of complete PLINK triplets
  # ----------------------------------------------------------
  roots <- unique(
    c(
      STEP8_ROOT,
      file.path(
        DATA_ROOT,
        "STEP9_SMR"
      ),
      file.path(
        DATA_ROOT,
        "STEP3_MR"
      ),
      DATA_ROOT
    )
  )

  bed_files <- character(0)

  for (
    rr in roots[
      dir.exists(
        roots
      )
    ]
  ) {

    z <- tryCatch(
      list.files(
        rr,
        pattern = "\\.bed$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      ),
      error = function(e) character(0)
    )

    bed_files <- unique(
      c(
        bed_files,
        z
      )
    )
  }

  scanned_prefixes <- unique(
    sub(
      "\\.bed$",
      "",
      bed_files,
      ignore.case = TRUE
    )
  )

  scanned_prefixes <- scanned_prefixes[
    vapply(
      scanned_prefixes,
      prefix_ok,
      logical(1)
    )
  ]

  for (pr in scanned_prefixes) {
    add_candidate(
      pr,
      "recursive_scan"
    )
  }

  if (!length(candidate_rows)) {
    stop(
      paste0(
        "No valid STEP8-consistent 1000G EUR PLINK reference was found.\\n",
        "Auto-detection required:\\n",
        "  * complete .bed/.bim/.fam triplet\\n",
        "  * reference N 450-550\\n",
        "  * at least ",
        MIN_COMMON_SNPS,
        " variants in chr16:",
        REGION_START,
        "-",
        REGION_STOP,
        "\\n\\nIf you know the exact prefix, set:\\n",
        "PLINK_REF_PREFIX_MANUAL <- 'D:/.../reference_prefix'\\n",
        "(prefix only; do not append .bed/.bim/.fam)"
      )
    )
  }

  candidates <- unique(
    rbindlist(
      candidate_rows,
      fill = TRUE
    ),
    by = "prefix"
  )

  # Explicitly reward actual manifest provenance above recursive copies.
  candidates[
    grepl(
      "^STEP8C_manifest:",
      source
    ),
    path_score :=
      path_score + 1000L
  ]

  # data.table::setorder() accepts column names only, not expressions.
  # Create a temporary path-length column for the final tie-breaker.
  candidates[
    ,
    path_length :=
      nchar(
        prefix
      )
  ]

  setorder(
    candidates,
    -path_score,
    -n_target_region_variants,
    -n_variants,
    path_length
  )

  candidate_audit_file <- file.path(
    OUT_ROOT,
    "00_QC",
    "STEP10B_V4_reference_candidate_audit.csv"
  )

  fwrite(
    candidates,
    candidate_audit_file
  )

  best <- candidates[
    1
  ]

  # Keep the helper out of downstream provenance objects.
  candidates[
    ,
    path_length := NULL
  ]
  best[
    ,
    path_length := NULL
  ]

  cat(
    "\\n1000G EUR reference selected:\\n",
    best$prefix,
    "\\nReference N: ",
    best$n_reference,
    "\\nReference variants: ",
    best$n_variants,
    "\\nTarget-region variants: ",
    best$n_target_region_variants,
    "\\nSelection source: ",
    best$source,
    "\\n",
    sep = ""
  )

  list(
    prefix = best$prefix,
    n_reference = best$n_reference,
    selection = best$source,
    candidates = candidates
  )
}

REF_INFO <- find_reference()

REF_PREFIX <- REF_INFO$prefix
N_REF <- REF_INFO$n_reference

if (
  length(
    REF_INFO$candidates
  )
) {
  fwrite(
    REF_INFO$candidates,
    file.path(
      OUT_ROOT,
      "00_QC",
      "STEP10B_V4_available_LD_reference_candidates.csv"
    )
  )
}

# ============================================================
# 5. INPUT HARD QC
# ============================================================

if (!file.exists(AF_MA)) {
  stop(
    "Missing final AF.ma:\n",
    AF_MA
  )
}

if (!file.exists(BMI_MA)) {
  stop(
    "Missing final BMI.ma:\n",
    BMI_MA
  )
}

if (!prefix_ok(REF_PREFIX)) {
  stop(
    "Selected LD reference prefix is invalid:\n",
    REF_PREFIX
  )
}

cat(
  "\n====================================================\n",
  "SECONDARY AF-BMI 16p11 COLOC/SUSIE\n",
  "====================================================\n",
  "Region: chr",
  CHR,
  ":",
  REGION_START,
  "-",
  REGION_STOP,
  "\n",
  "Lead SNP: ",
  FUMA_LEAD_SNP,
  "\n",
  "PLINK: ",
  PLINK_EXE,
  "\n",
  "LD reference: ",
  REF_PREFIX,
  "\n",
  "Reference N: ",
  N_REF,
  "\n",
  "Reference selection: ",
  REF_INFO$selection,
  "\n",
  sep = ""
)

# ============================================================
# 6. READ REFERENCE BIM AND DEFINE REGION SNP UNIVERSE
# ============================================================

BIM <- fread(
  paste0(
    REF_PREFIX,
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

BIM[
  ,
  BIM_ORDER :=
    seq_len(.N)
]

REGION_BIM <- BIM[
  BIM$CHR == CHR &
  BIM$BP >= REGION_START &
  BIM$BP <= REGION_STOP
]

REGION_BIM <- REGION_BIM[
  !is.na(SNP) &
  nzchar(SNP)
]

REGION_BIM <- REGION_BIM[
  !duplicated(SNP)
]

if (
  nrow(REGION_BIM) <
  MIN_COMMON_SNPS
) {
  stop(
    "Too few reference SNPs in secondary locus: ",
    nrow(REGION_BIM)
  )
}

REGION_SNPS <- REGION_BIM$SNP

# ============================================================
# 7. MEMORY-SAFE REGION EXTRACTION FROM FINAL .ma FILES
# ============================================================

read_ma_subset <- function(
  ma_file,
  target_snps,
  chunk_n = 500000L
) {

  con <- file(
    ma_file,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  header <- readLines(
    con,
    n = 1L,
    warn = FALSE
  )

  if (!length(header)) {
    stop(
      "Empty .ma file: ",
      ma_file
    )
  }

  out_lines <- list()
  k <- 0L

  target_env <- new.env(
    hash = TRUE,
    parent = emptyenv()
  )

  for (s in target_snps) {
    assign(
      s,
      TRUE,
      envir = target_env
    )
  }

  total_rows <- 0L
  kept_rows <- 0L

  repeat {

    lines <- readLines(
      con,
      n = chunk_n,
      warn = FALSE
    )

    if (!length(lines)) {
      break
    }

    total_rows <- total_rows +
      length(lines)

    snp <- sub(
      "[[:space:]].*$",
      "",
      lines,
      perl = TRUE,
      useBytes = TRUE
    )

    keep <- vapply(
      snp,
      exists,
      logical(1),
      envir = target_env,
      inherits = FALSE
    )

    if (any(keep)) {

      k <- k + 1L

      out_lines[[k]] <- lines[
        keep
      ]

      kept_rows <- kept_rows +
        sum(keep)
    }

    cat(
      "\rReading ",
      basename(ma_file),
      ": ",
      format(
        total_rows,
        big.mark = ","
      ),
      " rows; matched ",
      kept_rows,
      sep = ""
    )
  }

  cat("\n")

  if (!length(out_lines)) {
    return(
      data.table()
    )
  }

  text_blob <- paste(
    c(
      header,
      unlist(
        out_lines,
        use.names = FALSE
      )
    ),
    collapse = "\n"
  )

  fread(
    text = text_blob,
    showProgress = FALSE,
    na.strings = c(
      "NA",
      "NaN",
      "nan",
      "."
    )
  )
}

AF_REGION_RAW <- read_ma_subset(
  AF_MA,
  REGION_SNPS
)

BMI_REGION_RAW <- read_ma_subset(
  BMI_MA,
  REGION_SNPS
)

required_ma_cols <- c(
  "SNP",
  "A1",
  "A2",
  "freq",
  "b",
  "se",
  "p",
  "n"
)

if (!all(
  required_ma_cols %in%
  names(AF_REGION_RAW)
)) {
  stop(
    "AF.ma subset lacks required columns."
  )
}

if (!all(
  required_ma_cols %in%
  names(BMI_REGION_RAW)
)) {
  stop(
    "BMI.ma subset lacks required columns."
  )
}

# Deduplicate conservatively using smallest P.
AF_REGION_RAW[
  ,
  p_num := safe_num(p)
]

BMI_REGION_RAW[
  ,
  p_num := safe_num(p)
]

setorder(
  AF_REGION_RAW,
  SNP,
  p_num
)

setorder(
  BMI_REGION_RAW,
  SNP,
  p_num
)

AF_REGION_RAW <- AF_REGION_RAW[
  !duplicated(SNP)
]

BMI_REGION_RAW <- BMI_REGION_RAW[
  !duplicated(SNP)
]

# ============================================================
# 8. ALLELE HARMONIZATION TO ORIGINAL REFERENCE BIM
# ============================================================

harmonize_to_reference <- function(
  g,
  ref,
  trait_label
) {

  # The original full reference BIM uses BIM_ORDER, whereas the
  # final PLINK subset created later in the workflow uses ORDER_FINAL.
  # Normalize both layouts to a single internal BIM_ORDER column.
  ref_local <- copy(
    ref
  )

  if (
    !"BIM_ORDER" %in%
      names(
        ref_local
      )
  ) {

    if (
      "ORDER_FINAL" %in%
        names(
          ref_local
        )
    ) {
      ref_local[
        ,
        BIM_ORDER :=
          ORDER_FINAL
      ]
    } else {
      ref_local[
        ,
        BIM_ORDER :=
          seq_len(
            .N
          )
      ]
    }
  }

  z <- merge(
    g,
    ref_local[
      ,
      .(
        SNP,
        CHR,
        BP,
        REF_A1 = toupper(REF_A1),
        REF_A2 = toupper(REF_A2),
        BIM_ORDER
      )
    ],
    by = "SNP"
  )

  z[
    ,
    `:=`(
      A1 = toupper(
        as.character(A1)
      ),
      A2 = toupper(
        as.character(A2)
      ),
      b_num = safe_num(b),
      se_num = safe_num(se),
      p_num = safe_num(p),
      freq_num = safe_num(freq),
      n_num = safe_num(n)
    )
  ]

  z[
    ,
    palindromic :=
      is_palindromic(
        A1,
        A2
      )
  ]

  z[
    ,
    `:=`(
      CA1 =
        complement_allele(
          A1
        ),
      CA2 =
        complement_allele(
          A2
        )
    )
  ]

  z[
    ,
    orient :=
      fifelse(
        A1 == REF_A1 &
        A2 == REF_A2,
        1L,

        fifelse(
          A1 == REF_A2 &
          A2 == REF_A1,
          -1L,

          fifelse(
            CA1 == REF_A1 &
            CA2 == REF_A2,
            1L,

            fifelse(
              CA1 == REF_A2 &
              CA2 == REF_A1,
              -1L,
              0L
            )
          )
        )
      )
  ]

  z[
    ,
    valid_numeric :=
      is.finite(b_num) &
      is.finite(se_num) &
      se_num > 0 &
      is.finite(p_num) &
      p_num >= 0 &
      p_num <= 1
  ]

  z[
    ,
    keep_harmonized :=
      !palindromic &
      orient != 0L &
      valid_numeric
  ]

  z[
    ,
    b_ref :=
      fifelse(
        orient == 1L,
        b_num,
        -b_num
      )
  ]

  z[
    ,
    freq_refA1 :=
      fifelse(
        is.finite(
          freq_num
        ),
        fifelse(
          orient == 1L,
          freq_num,
          1 -
            freq_num
        ),
        NA_real_
      )
  ]

  qc <- data.table(
    trait = trait_label,
    n_raw_region = nrow(g),
    n_with_reference = nrow(z),
    n_palindromic =
      sum(
        z$palindromic,
        na.rm = TRUE
      ),
    n_allele_unresolved =
      sum(
        z$orient == 0L,
        na.rm = TRUE
      ),
    n_invalid_numeric =
      sum(
        !z$valid_numeric,
        na.rm = TRUE
      ),
    n_harmonized =
      sum(
        z$keep_harmonized,
        na.rm = TRUE
      )
  )

  list(
    data = z[
      keep_harmonized == TRUE
    ],
    qc = qc
  )
}

AF_H0 <- harmonize_to_reference(
  AF_REGION_RAW,
  REGION_BIM,
  "AF"
)

BMI_H0 <- harmonize_to_reference(
  BMI_REGION_RAW,
  REGION_BIM,
  "BMI"
)

AF_H <- AF_H0$data
BMI_H <- BMI_H0$data

HARM_QC_INITIAL <- rbind(
  AF_H0$qc,
  BMI_H0$qc
)

# ============================================================
# 9. COMMON SNPs BEFORE LD SUBSETTING
# ============================================================

COMMON0 <- intersect(
  AF_H$SNP,
  BMI_H$SNP
)

COMMON0 <- intersect(
  COMMON0,
  REGION_BIM$SNP
)

COMMON0 <- REGION_BIM[
  SNP %in%
    COMMON0
][
  order(
    BIM_ORDER
  ),
  SNP
]

if (
  length(COMMON0) <
  MIN_COMMON_SNPS
) {
  stop(
    "Too few AF-BMI harmonized common SNPs before LD QC: ",
    length(COMMON0)
  )
}

EXTRACT0 <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "secondary_AF_BMI_common_preLD.snplist"
)

writeLines(
  COMMON0,
  EXTRACT0
)

# ============================================================
# 10. CREATE REFERENCE SUBSET WITH PLINK
# ============================================================

REF_SUBSET0 <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "ref_secondary_16p11_preQC"
)

PLINK_LOG0 <- file.path(
  OUT_ROOT,
  "07_LOGS",
  "plink_make_subset.log"
)

args_subset0 <- c(
  "--bfile",
  shQuote(
    REF_PREFIX
  ),
  "--chr",
  as.character(CHR),
  "--from-bp",
  as.character(REGION_START),
  "--to-bp",
  as.character(REGION_STOP),
  "--extract",
  shQuote(
    norm_path(
      EXTRACT0,
      TRUE
    )
  ),
  "--maf",
  as.character(REF_MAF_MIN),
  "--geno",
  as.character(PLINK_GENO_MAX),
  "--keep-allele-order",
  "--make-bed",
  "--out",
  shQuote(
    paste0(
      norm_path(
        dirname(
          REF_SUBSET0
        ),
        TRUE
      ),
      "/",
      basename(
        REF_SUBSET0
      )
    )
  )
)

status0 <- system2(
  PLINK_EXE,
  args = args_subset0,
  stdout = PLINK_LOG0,
  stderr = PLINK_LOG0
)

if (
  status0 != 0L ||
  !prefix_ok(
    REF_SUBSET0
  )
) {
  stop(
    "PLINK reference subsetting failed. Inspect:\n",
    PLINK_LOG0
  )
}

SUB0_BIM <- fread(
  paste0(
    REF_SUBSET0,
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

SUB0_BIM[
  ,
  BIM_ORDER :=
    seq_len(.N)
]

# Re-harmonize to the actual subset BIM.
AF_H1 <- harmonize_to_reference(
  AF_REGION_RAW,
  SUB0_BIM,
  "AF"
)

BMI_H1 <- harmonize_to_reference(
  BMI_REGION_RAW,
  SUB0_BIM,
  "BMI"
)

AF1 <- AF_H1$data
BMI1 <- BMI_H1$data

COMMON1 <- intersect(
  AF1$SNP,
  BMI1$SNP
)

COMMON1 <- SUB0_BIM[
  SNP %in%
    COMMON1
][
  order(
    BIM_ORDER
  ),
  SNP
]

if (
  length(COMMON1) <
  MIN_COMMON_SNPS
) {
  stop(
    "Too few common SNPs after PLINK MAF/geno QC: ",
    length(COMMON1)
  )
}

# ============================================================
# 11. FINAL LD BED SUBSET
# ============================================================

FINAL_SNPLIST <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "secondary_AF_BMI_final.snplist"
)

writeLines(
  COMMON1,
  FINAL_SNPLIST
)

REF_FINAL <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "ref_secondary_16p11_FINAL"
)

PLINK_LOG1 <- file.path(
  OUT_ROOT,
  "07_LOGS",
  "plink_make_final_subset.log"
)

args_final <- c(
  "--bfile",
  shQuote(
    plink_bfile_prefix(
      REF_SUBSET0
    )
  ),
  "--extract",
  shQuote(
    norm_path(
      FINAL_SNPLIST,
      TRUE
    )
  ),
  "--keep-allele-order",
  "--make-bed",
  "--out",
  shQuote(
    paste0(
      norm_path(
        dirname(
          REF_FINAL
        ),
        TRUE
      ),
      "/",
      basename(
        REF_FINAL
      )
    )
  )
)

status1 <- system2(
  PLINK_EXE,
  args = args_final,
  stdout = PLINK_LOG1,
  stderr = PLINK_LOG1
)

if (
  status1 != 0L ||
  !prefix_ok(
    REF_FINAL
  )
) {
  stop(
    "PLINK final reference subset failed. Inspect:\n",
    PLINK_LOG1
  )
}

FINAL_BIM <- fread(
  paste0(
    REF_FINAL,
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

FINAL_BIM[
  ,
  ORDER_FINAL :=
    seq_len(.N)
]

# ============================================================
# 12. FINAL HARMONIZATION AND ORDER
# ============================================================

AF_HF <- harmonize_to_reference(
  AF_REGION_RAW,
  FINAL_BIM,
  "AF"
)$data

BMI_HF <- harmonize_to_reference(
  BMI_REGION_RAW,
  FINAL_BIM,
  "BMI"
)$data

FINAL_COMMON <- intersect(
  AF_HF$SNP,
  BMI_HF$SNP
)

FINAL_COMMON <- FINAL_BIM[
  SNP %in%
    FINAL_COMMON
][
  order(
    ORDER_FINAL
  ),
  SNP
]

if (
  length(FINAL_COMMON) <
  MIN_COMMON_SNPS
) {
  stop(
    "Too few final common SNPs: ",
    length(FINAL_COMMON)
  )
}

# If final harmonization dropped anything, regenerate once.
if (
  length(FINAL_COMMON) !=
  nrow(FINAL_BIM)
) {

  writeLines(
    FINAL_COMMON,
    FINAL_SNPLIST
  )

  unlink(
    paste0(
      REF_FINAL,
      c(
        ".bed",
        ".bim",
        ".fam"
      )
    )
  )

  status1b <- system2(
    PLINK_EXE,
    args = c(
      "--bfile",
      shQuote(
        plink_bfile_prefix(
          REF_SUBSET0
        )
      ),
      "--extract",
      shQuote(
        norm_path(
          FINAL_SNPLIST,
          TRUE
        )
      ),
      "--keep-allele-order",
      "--make-bed",
      "--out",
      shQuote(
        paste0(
      norm_path(
        dirname(
          REF_FINAL
        ),
        TRUE
      ),
      "/",
      basename(
        REF_FINAL
      )
    )
      )
    ),
    stdout = PLINK_LOG1,
    stderr = PLINK_LOG1
  )

  if (
    status1b != 0L ||
    !prefix_ok(
      REF_FINAL
    )
  ) {
    stop(
      "Final PLINK regeneration failed."
    )
  }

  FINAL_BIM <- fread(
    paste0(
      REF_FINAL,
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

  FINAL_BIM[
    ,
    ORDER_FINAL :=
      seq_len(.N)
  ]

  AF_HF <- harmonize_to_reference(
    AF_REGION_RAW,
    FINAL_BIM,
    "AF"
  )$data

  BMI_HF <- harmonize_to_reference(
    BMI_REGION_RAW,
    FINAL_BIM,
    "BMI"
  )$data

  FINAL_COMMON <- FINAL_BIM$SNP
}

setkey(
  AF_HF,
  SNP
)

setkey(
  BMI_HF,
  SNP
)

AF_FINAL <- AF_HF[
  FINAL_COMMON
]

BMI_FINAL <- BMI_HF[
  FINAL_COMMON
]

if (
  any(
    AF_FINAL$SNP !=
      FINAL_COMMON
  ) ||
  any(
    BMI_FINAL$SNP !=
      FINAL_COMMON
  )
) {
  stop(
    "Final SNP ordering mismatch after harmonization."
  )
}

# ============================================================
# 13. PLINK REFERENCE MAF
# ============================================================

FREQ_PREFIX <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "ref_secondary_16p11_freq"
)

FREQ_LOG <- file.path(
  OUT_ROOT,
  "07_LOGS",
  "plink_freq.log"
)

freq_status <- system2(
  PLINK_EXE,
  args = c(
    "--bfile",
    shQuote(
      plink_bfile_prefix(
        REF_FINAL
      )
    ),
    "--freq",
    "--keep-allele-order",
    "--out",
    shQuote(
      paste0(
      norm_path(
        dirname(
          FREQ_PREFIX
        ),
        TRUE
      ),
      "/",
      basename(
        FREQ_PREFIX
      )
    )
    )
  ),
  stdout = FREQ_LOG,
  stderr = FREQ_LOG
)

if (
  freq_status != 0L ||
  !file.exists(
    paste0(
      FREQ_PREFIX,
      ".frq"
    )
  )
) {
  stop(
    "PLINK --freq failed."
  )
}

REF_FREQ <- fread(
  paste0(
    FREQ_PREFIX,
    ".frq"
  ),
  showProgress = FALSE
)

setnames(
  REF_FREQ,
  old = names(REF_FREQ),
  new = toupper(
    names(REF_FREQ)
  )
)

if (!all(
  c(
    "SNP",
    "MAF"
  ) %in%
  names(REF_FREQ)
)) {
  stop(
    "Unexpected PLINK .frq format."
  )
}

REF_FREQ[
  ,
  MAF_REF :=
    safe_num(
      MAF
    )
]

REF_FREQ <- REF_FREQ[
  ,
  .(
    SNP,
    MAF_REF
  )
]

# ============================================================
# 14. SIGNED LD MATRIX
# ============================================================

LD_PREFIX <- file.path(
  OUT_ROOT,
  "06_LD_WORK",
  "secondary_AF_BMI_LD"
)

LD_LOG <- file.path(
  OUT_ROOT,
  "07_LOGS",
  "plink_LD.log"
)

ld_status <- system2(
  PLINK_EXE,
  args = c(
    "--bfile",
    shQuote(
      plink_bfile_prefix(
        REF_FINAL
      )
    ),
    "--r",
    "square",
    "--keep-allele-order",
    "--out",
    shQuote(
      paste0(
      norm_path(
        dirname(
          LD_PREFIX
        ),
        TRUE
      ),
      "/",
      basename(
        LD_PREFIX
      )
    )
    )
  ),
  stdout = LD_LOG,
  stderr = LD_LOG
)

LD_FILE <- paste0(
  LD_PREFIX,
  ".ld"
)

if (
  ld_status != 0L ||
  !file.exists(
    LD_FILE
  )
) {
  stop(
    "PLINK signed LD calculation failed. Inspect:\n",
    LD_LOG
  )
}

R <- as.matrix(
  fread(
    LD_FILE,
    header = FALSE,
    showProgress = FALSE
  )
)

storage.mode(R) <- "double"

if (
  nrow(R) !=
    length(FINAL_COMMON) ||
  ncol(R) !=
    length(FINAL_COMMON)
) {
  stop(
    "LD matrix dimension mismatch: ",
    nrow(R),
    " x ",
    ncol(R),
    " vs ",
    length(FINAL_COMMON),
    " SNPs."
  )
}

if (any(!is.finite(R))) {
  stop(
    "LD matrix contains non-finite values."
  )
}

# Numerical cleanup only.
R <- (
  R +
  t(R)
) / 2

diag(R) <- 1

rownames(R) <- FINAL_COMMON
colnames(R) <- FINAL_COMMON

# ============================================================
# 15. BUILD FINAL COMMON VARIANT TABLE
# ============================================================

COMMON <- merge(
  FINAL_BIM[
    ,
    .(
      SNP,
      CHR,
      BP,
      REF_A1,
      REF_A2,
      ORDER_FINAL
    )
  ],
  REF_FREQ,
  by = "SNP",
  all.x = TRUE
)

COMMON <- merge(
  COMMON,
  AF_FINAL[
    ,
    .(
      SNP,
      AF_beta = b_ref,
      AF_se = se_num,
      AF_p = p_num,
      AF_freq_refA1 = freq_refA1,
      AF_n = n_num
    )
  ],
  by = "SNP"
)

COMMON <- merge(
  COMMON,
  BMI_FINAL[
    ,
    .(
      SNP,
      BMI_beta = b_ref,
      BMI_se = se_num,
      BMI_p = p_num,
      BMI_freq_refA1 = freq_refA1,
      BMI_n = n_num
    )
  ],
  by = "SNP"
)

setorder(
  COMMON,
  ORDER_FINAL
)

if (
  any(
    COMMON$SNP !=
      FINAL_COMMON
  )
) {
  stop(
    "COMMON table order differs from LD matrix."
  )
}

COMMON[
  ,
  AF_MAF :=
    pmin(
      AF_freq_refA1,
      1 -
        AF_freq_refA1
    )
]

COMMON[
  ,
  BMI_MAF :=
    pmin(
      BMI_freq_refA1,
      1 -
        BMI_freq_refA1
    )
]

COMMON[
  !is.finite(AF_MAF),
  AF_MAF := MAF_REF
]

COMMON[
  !is.finite(BMI_MAF),
  BMI_MAF := MAF_REF
]

# Final MAF sanity.
maf_ok <- is.finite(
  COMMON$AF_MAF
) &
  is.finite(
    COMMON$BMI_MAF
) &
  COMMON$AF_MAF > 0 &
  COMMON$AF_MAF < 0.5 &
  COMMON$BMI_MAF > 0 &
  COMMON$BMI_MAF < 0.5

if (!all(maf_ok)) {

  keep_snps <- COMMON$SNP[
    maf_ok
  ]

  if (
    length(keep_snps) <
    MIN_COMMON_SNPS
  ) {
    stop(
      "Too few SNPs after final MAF sanity filter."
    )
  }

  idx <- match(
    keep_snps,
    COMMON$SNP
  )

  COMMON <- COMMON[
    idx
  ]

  R <- R[
    idx,
    idx,
    drop = FALSE
  ]

  FINAL_COMMON <- COMMON$SNP

  rownames(R) <- FINAL_COMMON
  colnames(R) <- FINAL_COMMON
}

# ============================================================
# 16. SAMPLE SIZE QC
# ============================================================

AF_N_REGION <- median(
  COMMON$AF_n,
  na.rm = TRUE
)

BMI_N_REGION <- median(
  COMMON$BMI_n,
  na.rm = TRUE
)

if (
  !is.finite(AF_N_REGION)
) {
  AF_N_REGION <- AF_N_TOTAL
}

if (
  !is.finite(BMI_N_REGION)
) {
  stop(
    "BMI sample size is unavailable in the region."
  )
}

if (
  abs(
    AF_N_REGION -
      AF_N_TOTAL
  ) /
    AF_N_TOTAL >
    0.05
) {
  warning(
    paste0(
      "AF region median N (",
      AF_N_REGION,
      ") differs >5% from the published European total N (",
      AF_N_TOTAL,
      "). The coloc dataset will use published N."
    )
  )
}

# AF uses the published European total N because s must correspond
# to the same case/control cohort definition.
AF_N_USE <- AF_N_TOTAL
BMI_N_USE <- round(
  BMI_N_REGION
)

# ============================================================
# 17. COLOC DATASETS
# ============================================================

D_AF <- list(
  beta = COMMON$AF_beta,
  varbeta = COMMON$AF_se^2,
  snp = COMMON$SNP,
  position = COMMON$BP,
  MAF = COMMON$AF_MAF,
  type = "cc",
  s = AF_CASE_FRACTION,
  N = AF_N_USE,
  LD = R
)

D_BMI <- list(
  beta = COMMON$BMI_beta,
  varbeta = COMMON$BMI_se^2,
  snp = COMMON$SNP,
  position = COMMON$BP,
  MAF = COMMON$BMI_MAF,
  type = "quant",
  N = BMI_N_USE,
  LD = R
)

AF_CHECK <- capture.output(
  coloc::check_dataset(
    D_AF,
    req = "LD"
  )
)

BMI_CHECK <- capture.output(
  coloc::check_dataset(
    D_BMI,
    req = "LD"
  )
)

writeLines(
  AF_CHECK,
  file.path(
    OUT_ROOT,
    "00_QC",
    "check_dataset_AF.txt"
  )
)

writeLines(
  BMI_CHECK,
  file.path(
    OUT_ROOT,
    "00_QC",
    "check_dataset_BMI.txt"
  )
)

# ============================================================
# 18. COLOC.ABF — DEFAULT + CONSERVATIVE PRIORS
# ============================================================

ABF_DEFAULT <- coloc::coloc.abf(
  D_AF,
  D_BMI,
  p1 = P1,
  p2 = P2,
  p12 = P12_DEFAULT
)

ABF_CONSERVATIVE <- coloc::coloc.abf(
  D_AF,
  D_BMI,
  p1 = P1,
  p2 = P2,
  p12 = P12_CONSERVATIVE
)

abf_summary_to_dt <- function(
  obj,
  prior_label,
  p12
) {

  s <- obj$summary

  data.table(
    prior = prior_label,
    p1 = P1,
    p2 = P2,
    p12 = p12,
    nsnps =
      if (
        "nsnps" %in%
        names(s)
      ) {
        safe_num(
          s[["nsnps"]]
        )
      } else {
        length(
          D_AF$snp
        )
      },
    PP.H0 =
      safe_num(
        s[
          grep(
            "PP.H0",
            names(s)
          )[1]
        ]
      ),
    PP.H1 =
      safe_num(
        s[
          grep(
            "PP.H1",
            names(s)
          )[1]
        ]
      ),
    PP.H2 =
      safe_num(
        s[
          grep(
            "PP.H2",
            names(s)
          )[1]
        ]
      ),
    PP.H3 =
      safe_num(
        s[
          grep(
            "PP.H3",
            names(s)
          )[1]
        ]
      ),
    PP.H4 =
      safe_num(
        s[
          grep(
            "PP.H4",
            names(s)
          )[1]
        ]
      )
  )
}

ABF_SUMMARY <- rbind(
  abf_summary_to_dt(
    ABF_DEFAULT,
    "default",
    P12_DEFAULT
  ),
  abf_summary_to_dt(
    ABF_CONSERVATIVE,
    "conservative",
    P12_CONSERVATIVE
  )
)

ABF_SUMMARY[
  ,
  interpretation :=
    fifelse(
      PP.H4 >= 0.80 &
      PP.H4 > PP.H3,
      "H4_support",

      fifelse(
        PP.H3 >= 0.80 &
        PP.H3 > PP.H4,
        "H3_support",
        "non_dominant"
      )
    )
]

# ============================================================
# 19. ABF PRIOR-SENSITIVITY GRID
# ============================================================

P12_GRID <- c(
  1e-6,
  3e-6,
  1e-5
)

ABF_PRIOR_GRID <- rbindlist(
  lapply(
    P12_GRID,
    function(px) {

      obj <- coloc::coloc.abf(
        D_AF,
        D_BMI,
        p1 = P1,
        p2 = P2,
        p12 = px
      )

      z <- abf_summary_to_dt(
        obj,
        paste0(
          "p12_",
          format(
            px,
            scientific = TRUE
          )
        ),
        px
      )

      z
    }
  )
)

# ============================================================
# 20. LD MISMATCH DIAGNOSTIC BEFORE FORMAL SUSIE
# ============================================================

Z_AF <- COMMON$AF_beta /
  COMMON$AF_se

Z_BMI <- COMMON$BMI_beta /
  COMMON$BMI_se

S_EST_AF <- tryCatch(
  susieR::estimate_s_rss(
    z = Z_AF,
    R = R,
    n = AF_N_USE,
    method = "null-mle"
  ),
  error = function(e) NA_real_
)

S_EST_BMI <- tryCatch(
  susieR::estimate_s_rss(
    z = Z_BMI,
    R = R,
    n = BMI_N_USE,
    method = "null-mle"
  ),
  error = function(e) NA_real_
)

# ============================================================
# 21. BOUNDED MISMATCH-AWARE SUSIE
# ============================================================
#
# IMPORTANT:
# coloc::runsusie(repeat_until_convergence=TRUE) can keep expanding
# max iterations dramatically when the fit does not converge.
# For this secondary/post hoc analysis we do NOT force convergence.
#
# Frozen rule:
#   - one bounded formal attempt (maxit = 500)
#   - no automatic iteration escalation
#   - if not converged, preserve a direct susie_rss diagnostic fit
#     and classify the fine-mapping result as unresolved
#   - do NOT run coloc.susie on a non-converged fit
# ============================================================

SUSIE_MAXIT_BOUNDED <- 500L
SUSIE_DIAG_MAXIT <- 200L

run_bounded_susie <- function(
  d,
  z,
  n_use,
  trait_label
) {

  cat(
    "\nRunning bounded mismatch-aware SuSiE for ",
    trait_label,
    " (maxit=",
    SUSIE_MAXIT_BOUNDED,
    "; no automatic escalation)...\n",
    sep = ""
  )

  formal_fit <- tryCatch(
    coloc::runsusie(
      d,
      suffix = paste0(
        trait_label,
        "_secondary_16p11"
      ),
      maxit = SUSIE_MAXIT_BOUNDED,
      repeat_until_convergence = FALSE,
      L = min(
        SUSIE_L,
        length(
          FINAL_COMMON
        )
      ),
      coverage = SUSIE_COVERAGE,
      R_finite = N_REF,
      R_mismatch = "eb"
    ),
    error = function(e) {
      structure(
        list(
          message = conditionMessage(e)
        ),
        class = "bounded_runsusie_error"
      )
    }
  )

  if (
    !inherits(
      formal_fit,
      "bounded_runsusie_error"
    ) &&
    isTRUE(
      formal_fit$converged
    )
  ) {
    return(
      list(
        fit = formal_fit,
        formal_valid = TRUE,
        status = "converged",
        error_message = NA_character_
      )
    )
  }

  formal_error <- if (
    inherits(
      formal_fit,
      "bounded_runsusie_error"
    )
  ) {
    formal_fit$message
  } else {
    "runsusie returned a non-converged fit"
  }

  cat(
    "\n",
    trait_label,
    " did not converge within the bounded formal attempt.\n",
    "Running a bounded direct susie_rss diagnostic fit only.\n",
    sep = ""
  )

  diagnostic_fit <- susieR::susie_rss(
    z = z,
    R = R,
    n = n_use,
    L = min(
      SUSIE_L,
      length(
        FINAL_COMMON
      )
    ),
    coverage = SUSIE_COVERAGE,
    R_finite = N_REF,
    R_mismatch = "eb",
    max_iter = SUSIE_DIAG_MAXIT,
    verbose = FALSE
  )

  list(
    fit = diagnostic_fit,
    formal_valid = FALSE,
    status = if (
      isTRUE(
        diagnostic_fit$converged
      )
    ) {
      "runsusie_failed_but_diagnostic_fit_converged"
    } else {
      "unresolved_nonconvergence"
    },
    error_message = formal_error
  )
}

SUSIE_AF_RUN <- run_bounded_susie(
  D_AF,
  Z_AF,
  AF_N_USE,
  "AF"
)

SUSIE_BMI_RUN <- run_bounded_susie(
  D_BMI,
  Z_BMI,
  BMI_N_USE,
  "BMI"
)

SUSIE_AF <- SUSIE_AF_RUN$fit
SUSIE_BMI <- SUSIE_BMI_RUN$fit

AF_SUSIE_FORMAL_VALID <- isTRUE(
  SUSIE_AF_RUN$formal_valid
)

BMI_SUSIE_FORMAL_VALID <- isTRUE(
  SUSIE_BMI_RUN$formal_valid
)

SUSIE_RUN_STATUS <- data.table(
  trait = c(
    "AF",
    "BMI"
  ),
  formal_valid = c(
    AF_SUSIE_FORMAL_VALID,
    BMI_SUSIE_FORMAL_VALID
  ),
  status = c(
    SUSIE_AF_RUN$status,
    SUSIE_BMI_RUN$status
  ),
  bounded_maxit = SUSIE_MAXIT_BOUNDED,
  diagnostic_maxit = SUSIE_DIAG_MAXIT,
  error_message = c(
    SUSIE_AF_RUN$error_message,
    SUSIE_BMI_RUN$error_message
  )
)

# ============================================================
# 22. FORMAL COLOC.SUSIE — ONLY IF BOTH TRAITS CONVERGED
# ============================================================

if (
  AF_SUSIE_FORMAL_VALID &&
  BMI_SUSIE_FORMAL_VALID
) {

  SUSIE_COLOC_DEFAULT <- coloc::coloc.susie(
    SUSIE_AF,
    SUSIE_BMI,
    p1 = P1,
    p2 = P2,
    p12 = P12_DEFAULT
  )

  SUSIE_COLOC_CONS <- coloc::coloc.susie(
    SUSIE_AF,
    SUSIE_BMI,
    p1 = P1,
    p2 = P2,
    p12 = P12_CONSERVATIVE
  )

} else {

  cat(
    "\nFormal coloc.susie skipped because at least one trait did not ",
    "produce a valid converged SuSiE fit.\n",
    "This is recorded as fine-mapping unresolved, not as H3 or H4.\n",
    sep = ""
  )

  SUSIE_COLOC_DEFAULT <- list(
    summary = data.frame()
  )

  SUSIE_COLOC_CONS <- list(
    summary = data.frame()
  )
}

standardize_susie_summary <- function(
  obj,
  prior_label,
  p12
) {

  # Empty/absent coloc.susie output = no formal signal-pair table.
  if (
    is.null(obj$summary) ||
    nrow(as.data.frame(obj$summary)) == 0L
  ) {
    return(
      data.table(
        prior = prior_label,
        p12 = p12,
        signal_pair = NA_integer_,
        hit1 = NA_character_,
        hit2 = NA_character_,
        PP.H0 = NA_real_,
        PP.H1 = NA_real_,
        PP.H2 = NA_real_,
        PP.H3 = NA_real_,
        PP.H4 = NA_real_
      )
    )
  }

  x <- as.data.table(
    obj$summary
  )

  first_match <- function(pattern) {
    z <- grep(
      pattern,
      names(x),
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(z)) z[1] else NA_character_
  }

  numeric_col_or_na <- function(colname) {
    if (
      is.na(colname) ||
      !colname %in% names(x)
    ) {
      return(
        rep(
          NA_real_,
          nrow(x)
        )
      )
    }
    safe_num(
      x[[colname]]
    )
  }

  char_col_or_na <- function(colname) {
    if (
      is.na(colname) ||
      !colname %in% names(x)
    ) {
      return(
        rep(
          NA_character_,
          nrow(x)
        )
      )
    }
    as.character(
      x[[colname]]
    )
  }

  h0 <- first_match("PP[.]H0")
  h1 <- first_match("PP[.]H1")
  h2 <- first_match("PP[.]H2")
  h3 <- first_match("PP[.]H3")
  h4 <- first_match("PP[.]H4")
  hit1_col <- first_match("^hit1$|hit1")
  hit2_col <- first_match("^hit2$|hit2")

  data.table(
    prior = prior_label,
    p12 = p12,
    signal_pair = seq_len(
      nrow(x)
    ),
    hit1 = char_col_or_na(
      hit1_col
    ),
    hit2 = char_col_or_na(
      hit2_col
    ),
    PP.H0 = numeric_col_or_na(
      h0
    ),
    PP.H1 = numeric_col_or_na(
      h1
    ),
    PP.H2 = numeric_col_or_na(
      h2
    ),
    PP.H3 = numeric_col_or_na(
      h3
    ),
    PP.H4 = numeric_col_or_na(
      h4
    )
  )
}

SUSIE_SUM_DEFAULT <- standardize_susie_summary(
  SUSIE_COLOC_DEFAULT,
  "default",
  P12_DEFAULT
)

SUSIE_SUM_CONS <- standardize_susie_summary(
  SUSIE_COLOC_CONS,
  "conservative",
  P12_CONSERVATIVE
)

SUSIE_SUMMARY_ALL <- rbind(
  SUSIE_SUM_DEFAULT,
  SUSIE_SUM_CONS,
  fill = TRUE
)

max_or_na <- function(x) {

  z <- x[
    is.finite(x)
  ]

  if (!length(z)) {
    return(NA_real_)
  }

  max(z)
}

n_signal_pairs_default <- sum(
  is.finite(
    SUSIE_SUM_DEFAULT$PP.H4
  ) |
  is.finite(
    SUSIE_SUM_DEFAULT$PP.H3
  )
)

n_signal_pairs_cons <- sum(
  is.finite(
    SUSIE_SUM_CONS$PP.H4
  ) |
  is.finite(
    SUSIE_SUM_CONS$PP.H3
  )
)

MAX_H4_DEFAULT <- max_or_na(
  SUSIE_SUM_DEFAULT$PP.H4
)

MAX_H3_DEFAULT <- max_or_na(
  SUSIE_SUM_DEFAULT$PP.H3
)

MAX_H4_CONS <- max_or_na(
  SUSIE_SUM_CONS$PP.H4
)

MAX_H3_CONS <- max_or_na(
  SUSIE_SUM_CONS$PP.H3
)

FORMAL_CLASS <-
  if (
    !AF_SUSIE_FORMAL_VALID ||
    !BMI_SUSIE_FORMAL_VALID
  ) {

    "Fine_mapping_unresolved_nonconvergence"

  } else if (
    n_signal_pairs_default == 0L
  ) {

    "No_signal_pair"

  } else if (
    is.finite(
      MAX_H3_DEFAULT
    ) &&
    MAX_H3_DEFAULT >= 0.80 &&
    (
      !is.finite(
        MAX_H4_DEFAULT
      ) ||
      MAX_H4_DEFAULT < 0.20
    )
  ) {

    "H3_dominant"

  } else if (
    is.finite(
      MAX_H4_DEFAULT
    ) &&
    MAX_H4_DEFAULT >= 0.80 &&
    is.finite(
      MAX_H4_CONS
    ) &&
    MAX_H4_CONS >= 0.80
  ) {

    "H4_robust_both_priors"

  } else if (
    is.finite(
      MAX_H4_DEFAULT
    ) &&
    MAX_H4_DEFAULT >= 0.80
  ) {

    "H4_prior_sensitive"

  } else {

    "No_dominant_H3_H4"
  }

# ============================================================
# 23. SUSIE CREDIBLE SET EXTRACTION
# ============================================================

extract_cs <- function(
  fit,
  trait_label,
  common_table
) {

  cs <- fit$sets$cs

  if (
    is.null(cs) ||
    !length(cs)
  ) {
    return(
      data.table()
    )
  }

  out <- list()

  for (nm in names(cs)) {

    idx <- cs[[nm]]

    idx <- idx[
      idx >= 1L &
      idx <= nrow(
        common_table
      )
    ]

    if (!length(idx)) {
      next
    }

    pip <- fit$pip[
      idx
    ]

    z <- common_table[
      idx,
      .(
        SNP,
        BP
      )
    ]

    z[
      ,
      `:=`(
        trait = trait_label,
        CS = nm,
        PIP = pip
      )
    ]

    out[[
      length(out) + 1L
    ]] <- z
  }

  if (!length(out)) {
    return(
      data.table()
    )
  }

  rbindlist(
    out,
    fill = TRUE
  )
}

CS_AF <- if (
  AF_SUSIE_FORMAL_VALID
) {
  extract_cs(
    SUSIE_AF,
    "AF",
    COMMON
  )
} else {
  data.table()
}

CS_BMI <- if (
  BMI_SUSIE_FORMAL_VALID
) {
  extract_cs(
    SUSIE_BMI,
    "BMI",
    COMMON
  )
} else {
  data.table()
}

# ============================================================
# 24. MISMATCH DIAGNOSTIC EXTRACTION
# ============================================================

diag_scalar <- function(
  fit,
  nm
) {

  d <- fit$R_finite_diagnostics

  if (
    is.null(d) ||
    !nm %in%
      names(d)
  ) {
    return(NA)
  }

  x <- d[[nm]]

  if (
    length(x) == 1L
  ) {
    return(x)
  }

  NA
}

SUSIE_DIAGNOSTICS <- data.table(
  trait = c(
    "AF",
    "BMI"
  ),
  n_reference = N_REF,
  estimate_s_rss = c(
    safe_num(
      S_EST_AF
    ),
    safe_num(
      S_EST_BMI
    )
  ),
  converged = c(
    isTRUE(
      SUSIE_AF$converged
    ),
    isTRUE(
      SUSIE_BMI$converged
    )
  ),
  n_credible_sets = c(
    if (
      is.null(
        SUSIE_AF$sets$cs
      )
    ) {
      0L
    } else {
      length(
        SUSIE_AF$sets$cs
      )
    },
    if (
      is.null(
        SUSIE_BMI$sets$cs
      )
    ) {
      0L
    } else {
      length(
        SUSIE_BMI$sets$cs
      )
    }
  ),
  R_reliability_flag = c(
    as.character(
      diag_scalar(
        SUSIE_AF,
        "R_reliability_flag"
      )
    ),
    as.character(
      diag_scalar(
        SUSIE_BMI,
        "R_reliability_flag"
      )
    )
  ),
  R_sensitivity_flag = c(
    as.character(
      diag_scalar(
        SUSIE_AF,
        "R_sensitivity_flag"
      )
    ),
    as.character(
      diag_scalar(
        SUSIE_BMI,
        "R_sensitivity_flag"
      )
    )
  ),
  Q_art = c(
    safe_num(
      diag_scalar(
        SUSIE_AF,
        "Q_art"
      )
    ),
    safe_num(
      diag_scalar(
        SUSIE_BMI,
        "Q_art"
      )
    )
  ),
  r_over_B = c(
    safe_num(
      diag_scalar(
        SUSIE_AF,
        "r_over_B"
      )
    ),
    safe_num(
      diag_scalar(
        SUSIE_BMI,
        "r_over_B"
      )
    )
  ),
  lambda_bias = c(
    safe_num(
      diag_scalar(
        SUSIE_AF,
        "lambda_bias"
      )
    ),
    safe_num(
      diag_scalar(
        SUSIE_BMI,
        "lambda_bias"
      )
    )
  ),
  B_corrected = c(
    safe_num(
      diag_scalar(
        SUSIE_AF,
        "B_corrected"
      )
    ),
    safe_num(
      diag_scalar(
        SUSIE_BMI,
        "B_corrected"
      )
    )
  )
)

# ============================================================
# 25. SER SENSITIVITY AUDIT
# ============================================================

run_ser <- function(
  z,
  n,
  trait_label
) {

  fit <- susieR::susie_rss(
    z = z,
    R = R,
    n = n,
    L = 1L,
    coverage = 0.95,
    R_finite = N_REF,
    R_mismatch = "eb",
    max_iter = 100,
    verbose = FALSE
  )

  alpha <- fit$alpha

  if (
    is.null(
      alpha
    )
  ) {
    return(
      list(
        fit = fit,
        summary = data.table(
          trait = trait_label,
          SER_top_snp = NA_character_,
          SER_top_PIP = NA_real_,
          SER_CS_size = 0L,
          SER_CS_snps = ""
        ),
        cs = character()
      )
    )
  }

  if (
    is.matrix(alpha)
  ) {
    pip <- as.numeric(
      alpha[1, ]
    )
  } else {
    pip <- as.numeric(
      alpha
    )
  }

  ord <- order(
    pip,
    decreasing = TRUE
  )

  cum <- cumsum(
    pip[
      ord
    ]
  )

  crossing <- which(
    cum >= 0.95
  )

  k <- if (
    length(crossing)
  ) {
    crossing[1]
  } else {
    length(ord)
  }

  cs_idx <- ord[
    seq_len(k)
  ]

  cs_snps <- COMMON$SNP[
    cs_idx
  ]

  top_idx <- ord[1]

  list(
    fit = fit,
    summary = data.table(
      trait = trait_label,
      SER_top_snp =
        COMMON$SNP[
          top_idx
        ],
      SER_top_PIP =
        pip[
          top_idx
        ],
      SER_CS_size =
        length(
          cs_snps
        ),
      SER_CS_snps =
        paste(
          cs_snps,
          collapse = ";"
        )
    ),
    cs = cs_snps
  )
}

SER_AF <- run_ser(
  Z_AF,
  AF_N_USE,
  "AF"
)

SER_BMI <- run_ser(
  Z_BMI,
  BMI_N_USE,
  "BMI"
)

cs_union <- function(
  fit
) {

  cs <- fit$sets$cs

  if (
    is.null(cs) ||
    !length(cs)
  ) {
    return(
      character()
    )
  }

  unique(
    COMMON$SNP[
      unlist(
        cs,
        use.names = FALSE
      )
    ]
  )
}

EB_CS_AF <- if (
  AF_SUSIE_FORMAL_VALID
) {
  cs_union(
    SUSIE_AF
  )
} else {
  character()
}

EB_CS_BMI <- if (
  BMI_SUSIE_FORMAL_VALID
) {
  cs_union(
    SUSIE_BMI
  )
} else {
  character()
}

jaccard <- function(
  a,
  b
) {

  if (
    !length(a) ||
    !length(b)
  ) {
    return(NA_real_)
  }

  length(
    intersect(
      a,
      b
    )
  ) /
    length(
      union(
        a,
        b
      )
    )
}

SER_AUDIT <- rbind(
  cbind(
    SER_AF$summary,
    data.table(
      partner_trait = "BMI",
      SER_top_in_partner_EB_CS =
        SER_AF$summary$SER_top_snp %in%
          EB_CS_BMI,
      SER_partner_Jaccard =
        jaccard(
          SER_AF$cs,
          EB_CS_BMI
        )
    )
  ),
  cbind(
    SER_BMI$summary,
    data.table(
      partner_trait = "AF",
      SER_top_in_partner_EB_CS =
        SER_BMI$summary$SER_top_snp %in%
          EB_CS_AF,
      SER_partner_Jaccard =
        jaccard(
          SER_BMI$cs,
          EB_CS_AF
        )
    )
  )
)

# ============================================================
# 26. FINAL FORMAL SUMMARY
# ============================================================

FINAL_SUMMARY <- data.table(
  analysis = "secondary_AF_BMI_chr16_28538336_29008079",
  status = "secondary_hypothesis_generating_followup",
  chromosome = CHR,
  start = REGION_START,
  stop = REGION_STOP,
  FUMA_lead_SNP = FUMA_LEAD_SNP,
  n_common_snps = nrow(COMMON),
  AF_N = AF_N_USE,
  AF_cases = AF_CASES,
  AF_controls = AF_CONTROLS,
  AF_case_fraction = AF_CASE_FRACTION,
  BMI_N_region_median = BMI_N_USE,
  LD_reference_prefix = REF_PREFIX,
  LD_reference_N = N_REF,
  LD_reference_selection = REF_INFO$selection,
  AF_SuSiE_formal_valid =
    AF_SUSIE_FORMAL_VALID,
  BMI_SuSiE_formal_valid =
    BMI_SUSIE_FORMAL_VALID,
  AF_SuSiE_status =
    SUSIE_AF_RUN$status,
  BMI_SuSiE_status =
    SUSIE_BMI_RUN$status,
  ABF_H3_default =
    ABF_SUMMARY[
      prior ==
        "default",
      PP.H3
    ],
  ABF_H4_default =
    ABF_SUMMARY[
      prior ==
        "default",
      PP.H4
    ],
  ABF_H3_conservative =
    ABF_SUMMARY[
      prior ==
        "conservative",
      PP.H3
    ],
  ABF_H4_conservative =
    ABF_SUMMARY[
      prior ==
        "conservative",
      PP.H4
    ],
  SuSiE_signal_pairs_default =
    n_signal_pairs_default,
  SuSiE_max_H3_default =
    MAX_H3_DEFAULT,
  SuSiE_max_H4_default =
    MAX_H4_DEFAULT,
  SuSiE_signal_pairs_conservative =
    n_signal_pairs_cons,
  SuSiE_max_H3_conservative =
    MAX_H3_CONS,
  SuSiE_max_H4_conservative =
    MAX_H4_CONS,
  SuSiE_formal_class =
    FORMAL_CLASS
)

# ============================================================
# 27. QC SUMMARY
# ============================================================

QC_SUMMARY <- data.table(
  item = c(
    "reference_region_snps",
    "AF_ma_region_rows",
    "BMI_ma_region_rows",
    "AF_initial_harmonized",
    "BMI_initial_harmonized",
    "common_preLD",
    "common_after_ref_MAF_geno_QC",
    "final_common_snps",
    "lead_SNP_present_final",
    "LD_matrix_rows",
    "LD_reference_N",
    "AF_N_used",
    "BMI_N_used",
    "AF_s_estimate",
    "BMI_s_estimate",
    "AF_susie_converged",
    "BMI_susie_converged",
    "AF_CS_n",
    "BMI_CS_n",
    "formal_SuSiE_class"
  ),
  value = c(
    nrow(REGION_BIM),
    nrow(AF_REGION_RAW),
    nrow(BMI_REGION_RAW),
    HARM_QC_INITIAL[
      trait == "AF",
      n_harmonized
    ],
    HARM_QC_INITIAL[
      trait == "BMI",
      n_harmonized
    ],
    length(COMMON0),
    nrow(SUB0_BIM),
    nrow(COMMON),
    FUMA_LEAD_SNP %in%
      COMMON$SNP,
    nrow(R),
    N_REF,
    AF_N_USE,
    BMI_N_USE,
    S_EST_AF,
    S_EST_BMI,
    isTRUE(
      SUSIE_AF$converged
    ),
    isTRUE(
      SUSIE_BMI$converged
    ),
    if (
      is.null(
        SUSIE_AF$sets$cs
      )
    ) {
      0L
    } else {
      length(
        SUSIE_AF$sets$cs
      )
    },
    if (
      is.null(
        SUSIE_BMI$sets$cs
      )
    ) {
      0L
    } else {
      length(
        SUSIE_BMI$sets$cs
      )
    },
    FORMAL_CLASS
  )
)

# ============================================================
# 28. EXPORT TABLES
# ============================================================

COMMON_FILE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_secondary_AF_BMI_common_harmonized_variants.csv.gz"
)

fwrite(
  COMMON,
  COMMON_FILE
)

HARM_QC_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V4_harmonization_QC.csv"
)

fwrite(
  HARM_QC_INITIAL,
  HARM_QC_FILE
)

QC_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V4_QC_summary.csv"
)

fwrite(
  QC_SUMMARY,
  QC_FILE
)

RUN_STATUS_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V4_SuSiE_bounded_run_status.csv"
)

fwrite(
  SUSIE_RUN_STATUS,
  RUN_STATUS_FILE
)

DIAG_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V4_SuSiE_LD_mismatch_diagnostics.csv"
)

fwrite(
  SUSIE_DIAGNOSTICS,
  DIAG_FILE
)

ABF_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table_secondary_AF_BMI_coloc_ABF_summary.csv"
)

fwrite(
  ABF_SUMMARY,
  ABF_FILE
)

SUSIE_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table_secondary_AF_BMI_coloc_SuSiE_signal_pairs.csv"
)

fwrite(
  SUSIE_SUMMARY_ALL,
  SUSIE_FILE
)

FINAL_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table_secondary_AF_BMI_FINAL_summary.csv"
)

fwrite(
  FINAL_SUMMARY,
  FINAL_FILE
)

PRIOR_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_secondary_AF_BMI_ABF_prior_sensitivity.csv"
)

fwrite(
  ABF_PRIOR_GRID,
  PRIOR_FILE
)

CS_AF_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_secondary_AF_SuSiE_credible_sets.csv"
)

fwrite(
  CS_AF,
  CS_AF_FILE
)

CS_BMI_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_secondary_BMI_SuSiE_credible_sets.csv"
)

fwrite(
  CS_BMI,
  CS_BMI_FILE
)

SER_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_secondary_AF_BMI_SER_audit.csv"
)

fwrite(
  SER_AUDIT,
  SER_FILE
)

# Save formal objects for reproducibility.
saveRDS(
  SUSIE_AF,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SUSIE_AF_formal_fit.rds"
  )
)

saveRDS(
  SUSIE_BMI,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SUSIE_BMI_formal_fit.rds"
  )
)

saveRDS(
  SUSIE_COLOC_DEFAULT,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "COLOC_SUSIE_default.rds"
  )
)

saveRDS(
  SUSIE_COLOC_CONS,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "COLOC_SUSIE_conservative.rds"
  )
)

# ============================================================
# 29. FIGURE SOURCE DATA
# ============================================================

FIG_AF_SOURCE <- COMMON[
  ,
  .(
    SNP,
    BP,
    Mb = BP / 1e6,
    p = AF_p,
    minus_log10_p =
      -log10(
        pmax(
          AF_p,
          1e-300
        )
      ),
    PIP =
      if (
        AF_SUSIE_FORMAL_VALID
      ) {
        SUSIE_AF$pip
      } else {
        rep(
          NA_real_,
          .N
        )
      },
    formal_SuSiE_valid =
      AF_SUSIE_FORMAL_VALID,
    in_CS =
      SNP %in%
        EB_CS_AF,
    lead_snp =
      SNP ==
        FUMA_LEAD_SNP
  )
]

FIG_BMI_SOURCE <- COMMON[
  ,
  .(
    SNP,
    BP,
    Mb = BP / 1e6,
    p = BMI_p,
    minus_log10_p =
      -log10(
        pmax(
          BMI_p,
          1e-300
        )
      ),
    PIP =
      if (
        BMI_SUSIE_FORMAL_VALID
      ) {
        SUSIE_BMI$pip
      } else {
        rep(
          NA_real_,
          .N
        )
      },
    formal_SuSiE_valid =
      BMI_SUSIE_FORMAL_VALID,
    in_CS =
      SNP %in%
        EB_CS_BMI,
    lead_snp =
      SNP ==
        FUMA_LEAD_SNP
  )
]

FIG_AF_FILE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_secondary_AF_region_and_PIP.csv"
)

FIG_BMI_FILE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_secondary_BMI_region_and_PIP.csv"
)

fwrite(
  FIG_AF_SOURCE,
  FIG_AF_FILE
)

fwrite(
  FIG_BMI_SOURCE,
  FIG_BMI_FILE
)

# ============================================================
# 30. MAIN FIGURES — REGIONAL ASSOCIATION
# ============================================================

make_assoc_plot <- function(
  dat,
  trait_label,
  filename
) {

  p <- ggplot(
    dat,
    aes(
      x = Mb,
      y = minus_log10_p
    )
  ) +
    geom_point(
      shape = 16,
      size = 1.45,
      alpha = 0.62,
      colour = "#9EA4AA"
    ) +
    geom_point(
      data = dat[
        in_CS == TRUE
      ],
      shape = 21,
      size = 2.1,
      stroke = 0.35,
      fill = "#6E919C",
      colour = "#2D454E"
    ) +
    geom_point(
      data = dat[
        lead_snp == TRUE
      ],
      shape = 23,
      size = 3.2,
      stroke = 0.55,
      fill = "#B77A70",
      colour = "#5D3630"
    ) +
    labs(
      x = "Chromosome 16 position (Mb)",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic(
      base_size = 11,
      base_family = "Arial"
    ) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      legend.position = "none",
      panel.border = element_blank(),
      plot.margin = margin(
        10,
        12,
        10,
        10
      )
    )

  ggsave(
    filename,
    p,
    width = 6.8,
    height = 4.6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

AF_ASSOC_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_secondary_AF_regional_association.tiff"
)

BMI_ASSOC_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_secondary_BMI_regional_association.tiff"
)

make_assoc_plot(
  FIG_AF_SOURCE,
  "AF",
  AF_ASSOC_TIFF
)

make_assoc_plot(
  FIG_BMI_SOURCE,
  "BMI",
  BMI_ASSOC_TIFF
)

# ============================================================
# 31. MAIN FIGURES — SUSIE PIP
# ============================================================

make_pip_plot <- function(
  dat,
  filename
) {

  p <- ggplot(
    dat,
    aes(
      x = Mb,
      y = PIP
    )
  ) +
    geom_segment(
      aes(
        xend = Mb,
        y = 0,
        yend = PIP
      ),
      linewidth = 0.35,
      colour = "#A9AEB4"
    ) +
    geom_point(
      shape = 21,
      size = 1.8,
      stroke = 0.3,
      fill = "#7B98A2",
      colour = "#31454C"
    ) +
    geom_point(
      data = dat[
        lead_snp == TRUE
      ],
      shape = 23,
      size = 3.3,
      stroke = 0.55,
      fill = "#B77A70",
      colour = "#5D3630"
    ) +
    labs(
      x = "Chromosome 16 position (Mb)",
      y = "Posterior inclusion probability"
    ) +
    coord_cartesian(
      ylim = c(
        0,
        1
      )
    ) +
    theme_classic(
      base_size = 11,
      base_family = "Arial"
    ) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      legend.position = "none",
      panel.border = element_blank()
    )

  ggsave(
    filename,
    p,
    width = 6.8,
    height = 4.6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

AF_PIP_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_secondary_AF_SuSiE_PIP.tiff"
)

BMI_PIP_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_secondary_BMI_SuSiE_PIP.tiff"
)

if (
  AF_SUSIE_FORMAL_VALID
) {
  make_pip_plot(
    FIG_AF_SOURCE,
    AF_PIP_TIFF
  )
}

if (
  BMI_SUSIE_FORMAL_VALID
) {
  make_pip_plot(
    FIG_BMI_SOURCE,
    BMI_PIP_TIFF
  )
}

# ============================================================
# 32. SUPPLEMENTARY SIGNAL-PAIR HEATMAP
# ============================================================

PAIR_DEFAULT <- SUSIE_SUM_DEFAULT[
  is.finite(PP.H4)
]

if (nrow(PAIR_DEFAULT)) {

  PAIR_DEFAULT[
    ,
    hit1_plot :=
      fifelse(
        is.na(hit1) |
        !nzchar(hit1),
        paste0(
          "AF signal ",
          signal_pair
        ),
        hit1
      )
  ]

  PAIR_DEFAULT[
    ,
    hit2_plot :=
      fifelse(
        is.na(hit2) |
        !nzchar(hit2),
        paste0(
          "BMI signal ",
          signal_pair
        ),
        hit2
      )
  ]

  p_pair <- ggplot(
    PAIR_DEFAULT,
    aes(
      x = hit2_plot,
      y = hit1_plot,
      fill = PP.H4
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.2f",
          PP.H4
        )
      ),
      family = "Arial",
      size = 3.4,
      colour = "black"
    ) +
    scale_fill_gradient(
      low = "#F0F1F2",
      high = "#6E919C",
      limits = c(
        0,
        1
      )
    ) +
    labs(
      x = "BMI SuSiE signal",
      y = "AF SuSiE signal"
    ) +
    theme_classic(
      base_size = 11,
      base_family = "Arial"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      legend.position = "none",
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )

  ggsave(
    file.path(
      OUT_ROOT,
      "05_SUPPLEMENTARY_FIGURES",
      "FigureS_secondary_AF_BMI_SuSiE_signal_pair_H4.tiff"
    ),
    p_pair,
    width = 5.8,
    height = 5.0,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

# ============================================================
# 33. METHOD PROVENANCE
# ============================================================

PROVENANCE <- data.table(
  item = c(
    "analysis_status",
    "reason_for_secondary_followup",
    "region",
    "genome_build",
    "AF_source",
    "AF_type",
    "AF_cases",
    "AF_controls",
    "BMI_type",
    "LD_reference",
    "LD_reference_N",
    "allele_QC",
    "palindromic_rule",
    "SuSiE_R_finite",
    "SuSiE_R_mismatch",
    "ABF_default_p12",
    "ABF_conservative_p12",
    "SuSiE_default_p12",
    "SuSiE_conservative_p12",
    "formal_interpretation_rule"
  ),
  value = c(
    "Secondary hypothesis-generating locus-level follow-up",
    "Convergent conjFDR + LAVA + FUMA + eQTL/sQTL evidence after the pre-specified STEP8 analyses",
    paste0(
      "chr16:",
      REGION_START,
      "-",
      REGION_STOP
    ),
    "GRCh37/hg19",
    "GCST90624412 European AF GWAS",
    "case-control",
    as.character(
      AF_CASES
    ),
    as.character(
      AF_CONTROLS
    ),
    "quantitative",
    REF_PREFIX,
    as.character(
      N_REF
    ),
    "GWAS effect alleles oriented to PLINK reference A1; unresolved allele pairs removed",
    "All A/T and C/G palindromic SNPs excluded conservatively",
    as.character(
      N_REF
    ),
    "eb",
    as.character(
      P12_DEFAULT
    ),
    as.character(
      P12_CONSERVATIVE
    ),
    as.character(
      P12_DEFAULT
    ),
    as.character(
      P12_CONSERVATIVE
    ),
    "H4 robust if max PP.H4 >=0.80 under both priors; H4 prior-sensitive if >=0.80 only under default; H3-dominant if default PP.H3>=0.80 and PP.H4<0.20"
  )
)

PROV_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V4_method_provenance.csv"
)

fwrite(
  PROVENANCE,
  PROV_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_ROOT,
    "07_LOGS",
    "STEP10B_V4_sessionInfo.txt"
  )
)

# ============================================================
# 34. FINAL CONSOLE REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP10B_V4 SECONDARY AF-BMI 16p11 COLOC/SUSIE COMPLETE\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nFinal summary:\n"
)

print(
  FINAL_SUMMARY
)

cat(
  "\nABF prior sensitivity:\n"
)

print(
  ABF_SUMMARY
)

cat(
  "\nFormal SuSiE signal-pair summary:\n"
)

print(
  SUSIE_SUMMARY_ALL
)

cat(
  "\nLD mismatch diagnostics:\n"
)

print(
  SUSIE_DIAGNOSTICS
)

cat(
  "\nSER audit:\n"
)

print(
  SER_AUDIT
)

cat(
  "\nINTERPRETATION RULE:\n",
  "This is a SECONDARY follow-up. ",
  "It must not be presented as one of the eight pre-specified STEP8 tests.\n",
  "H4 supports shared causal-signal evidence; H3 supports distinct causal signals; ",
  "No signal pair means unresolved fine-mapping, not proof of no sharing.\n",
  sep = ""
)

cat(
  "\nUPLOAD THESE CORE OUTPUTS:\n",
  "1) ", QC_FILE, "\n",
  "2) ", DIAG_FILE, "\n",
  "3) ", FINAL_FILE, "\n",
  "4) ", ABF_FILE, "\n",
  "5) ", SUSIE_FILE, "\n",
  "6) ", SER_FILE, "\n",
  "7) ", CS_AF_FILE, "\n",
  "8) ", CS_BMI_FILE, "\n",
  "9) ", FIG_AF_FILE, "\n",
  "10) ", FIG_BMI_FILE, "\n",
  sep = ""
)

cat(
  "====================================================\n"
)
