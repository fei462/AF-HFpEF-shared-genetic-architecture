# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP3B_V2_MR_MVMR_22CHR_FIXED(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP 3B V2 — Bidirectional UVMR + MVMR (22-chromosome LD-reference fix)
# Project: AF – HFpEF – BMI – OSA shared genetic architecture
#
# PRIMARY PRINCIPLES
#   1) MiXeR is NOT required for this script. MR/MVMR can run independently.
#   2) Primary IV selection follows the reference article:
#        P < 5e-8
#        F >= 10
#        LD clumping r2 < 0.001, 10,000 kb, EUR 1000G
#   3) If the PRIMARY threshold yields too few independent IVs (<3),
#      the script automatically performs PRE-SPECIFIED exploratory relaxation:
#        Tier 1: P < 5e-8   (PRIMARY)
#        Tier 2: P < 1e-6   (EXPLORATORY)
#        Tier 3: P < 5e-6   (EXPLORATORY)
#
#      IMPORTANT:
#      - F >= 10 is NEVER relaxed.
#      - r2 < 0.001 and 10,000 kb are NEVER relaxed.
#      - The first relaxed tier reaching >=3 IVs is used for the exploratory MR.
#      - Primary and relaxed analyses are always kept separate.
#      - Relaxed results are NEVER silently substituted for the primary results.
#
#   4) UVMR:
#        Primary estimator = IVW multiplicative random effects (>=2 IVs)
#        1 IV              = Wald ratio
#        Sensitivity       = weighted median + MR-Egger (>=3 IVs)
#        Heterogeneity     = Cochran Q
#        Pleiotropy        = MR-Egger intercept
#        MR-PRESSO         = run when >=4 IVs and package available
#        Leave-one-out     = run when >=3 IVs
#
#   5) MVMR:
#        Exposures = BMI + OSA
#        Outcomes  = HFpEF, AF
#        IVW MVMR + conditional F + Q-statistic
#        gencov = 0 (appropriate when exposure GWAS are treated as non-overlapping)
#
#   6) Main visualization mirrors the reference article Figure 5:
#        squares = IVW/Wald point estimate
#        horizontal lines = 95% CI
#        muted red/pink = Bonferroni significant
#        blue            = nominal P < 0.05
#        gray            = non-significant
#        no legend inside figure; legend saved separately
#
# REQUIRED INPUTS (auto-discovered under D:/A/data)
#   AF    : GCST90624412.tsv.gz
#   HFpEF : FORMAT-METAL_Pheno4_EUR.tsv.gz
#   BMI   : SNP_gwas_mc_merge_nogc.tbl.uniq.gz
#   OSA   : finngen_R9_G6_SLEEPAPNO.gz
#
# OUTPUT ROOT
#   D:/A/data/STEP3_MR/
#
# NOTE
#   This script is intentionally conservative. If even P<5e-6 does not provide
#   enough independent strong IVs, the pair is reported as sparse/not estimable;
#   the script does NOT keep relaxing until a significant result appears.
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

DATA_DIR <- "D:/A/data"
OUT_ROOT <- file.path(DATA_DIR, "STEP3_MR_V2")

# Pre-specified instrument thresholds
P_TIERS <- c(
  primary_GWS = 5e-8,
  exploratory_P1e_6 = 1e-6,
  exploratory_P5e_6 = 5e-6
)

MIN_IV_TARGET <- 3L
F_MIN <- 10
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000

# Planned UVMR family = 6 tests.
# Keep denominator fixed at 6 even if one direction is not estimable.
UVMR_PLANNED_TESTS <- 6L

# Planned MVMR direct effects = BMI and OSA for two outcomes = 4.
MVMR_PLANNED_TESTS <- 4L

# Harmonisation:
# action=2 = infer strand using EAF for palindromic variants; conservative default.
HARMONISE_ACTION <- 2L

# MR-PRESSO permutations
MRPRESSO_NB_DISTRIBUTION <- 1000L

# Local LD reference auto-download if not already present.
ALLOW_REFERENCE_DOWNLOAD <- TRUE

# V2 LD-reference safety controls
# A candidate instrument must be found in the EUR reference.  If less than
# this fraction can be mapped, the run stops instead of silently producing
# a misleading clumped set.
MIN_REF_MATCH_RATE <- 0.50

# The local reference may be either:
#   (A) one whole-genome EUR.bed/bim/fam prefix; or
#   (B) a complete 22-chromosome family such as
#       1000G.EUR.QC.1 ... 1000G.EUR.QC.22
# V2 explicitly forbids treating one .1-.22 chromosome prefix as the
# whole-genome reference.

# Figure settings
FIG_DPI <- 600
BASE_FAMILY <- "Arial"


# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

DIR_REF   <- file.path(OUT_ROOT, "00_reference")
DIR_INST  <- file.path(OUT_ROOT, "01_instruments")
DIR_HARM  <- file.path(OUT_ROOT, "02_harmonized")
DIR_UVMR  <- file.path(OUT_ROOT, "03_uvmr")
DIR_MVMR  <- file.path(OUT_ROOT, "04_mvmr")
DIR_SENS  <- file.path(OUT_ROOT, "05_sensitivity")
DIR_FIG   <- file.path(OUT_ROOT, "06_figures")
DIR_REP   <- file.path(OUT_ROOT, "07_reports")
DIR_CACHE <- file.path(OUT_ROOT, "08_cache")

for (d in c(
  OUT_ROOT, DIR_REF, DIR_INST, DIR_HARM,
  DIR_UVMR, DIR_MVMR, DIR_SENS, DIR_FIG,
  DIR_REP, DIR_CACHE
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# 2. PACKAGE INSTALLATION
# ==============================================================================

options(repos = c(
  MRCIEU = "https://mrcieu.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg)
  }
}

for (p in c(
  "data.table", "DBI", "duckdb", "ggplot2", "scales",
  "ieugwasr", "remotes", "tidyr"
)) {
  install_cran_if_missing(p)
}

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  message("Installing TwoSampleMR...")
  try(
    install.packages(
      "TwoSampleMR",
      repos = c(
        "https://mrcieu.r-universe.dev",
        "https://cloud.r-project.org"
      )
    ),
    silent = TRUE
  )
}
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop(
    "TwoSampleMR installation failed. Please run:\n",
    "install.packages('TwoSampleMR', repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))"
  )
}

# MVMR is installed from GitHub if necessary.
MVMR_AVAILABLE <- requireNamespace("MVMR", quietly = TRUE)
if (!MVMR_AVAILABLE) {
  message("Installing MVMR from GitHub...")
  try(
    remotes::install_github("WSpiller/MVMR", upgrade = "never"),
    silent = TRUE
  )
  MVMR_AVAILABLE <- requireNamespace("MVMR", quietly = TRUE)
}

# MR-PRESSO is optional; UVMR continues if installation fails.
MRPRESSO_AVAILABLE <- requireNamespace("MRPRESSO", quietly = TRUE)
if (!MRPRESSO_AVAILABLE) {
  message("Trying to install MR-PRESSO...")
  try(
    remotes::install_github("rondolab/MR-PRESSO", upgrade = "never"),
    silent = TRUE
  )
  MRPRESSO_AVAILABLE <- requireNamespace("MRPRESSO", quietly = TRUE)
}

library(data.table)
library(DBI)
library(duckdb)
library(ggplot2)


# ==============================================================================
# 3. LOGGING
# ==============================================================================

LOG_FILE <- file.path(DIR_REP, "STEP3B_master_console.log")

log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", txt, "\n"),
    file = LOG_FILE,
    append = TRUE
  )
}

if (file.exists(LOG_FILE)) unlink(LOG_FILE)

log_msg("============================================================")
log_msg("STEP 3B V2 START — 22-chromosome-safe clumping")
log_msg("============================================================")


# ==============================================================================
# 4. FILE DISCOVERY
# ==============================================================================

all_files <- list.files(
  DATA_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

choose_shortest <- function(x) {
  if (length(x) == 0) return(NA_character_)
  x[which.min(nchar(x))]
}

find_exact_file <- function(fname, prefer_regex = NULL) {
  hit <- all_files[tolower(basename(all_files)) == tolower(fname)]

  if (!is.null(prefer_regex) && length(hit) > 0) {
    hp <- hit[grepl(prefer_regex, hit, ignore.case = TRUE)]
    if (length(hp) > 0) return(choose_shortest(hp))
  }

  choose_shortest(hit)
}

AF_PATH <- find_exact_file("GCST90624412.tsv.gz")

HF_PATH <- find_exact_file(
  "FORMAT-METAL_Pheno4_EUR.tsv.gz",
  prefer_regex = "STEP1B_HFpEF_PREP|HERMES"
)

BMI_PATH <- find_exact_file("SNP_gwas_mc_merge_nogc.tbl.uniq.gz")
OSA_PATH <- find_exact_file("finngen_R9_G6_SLEEPAPNO.gz")

GWAS_PATHS <- c(
  AF = AF_PATH,
  HFpEF = HF_PATH,
  BMI = BMI_PATH,
  OSA = OSA_PATH
)

log_msg("Resolved GWAS files:")
for (nm in names(GWAS_PATHS)) {
  log_msg("  ", nm, " = ", GWAS_PATHS[[nm]])
}

if (any(is.na(GWAS_PATHS)) || any(!file.exists(GWAS_PATHS))) {
  stop(
    "One or more required GWAS files were not found:\n",
    paste(names(GWAS_PATHS), GWAS_PATHS, sep = " = ", collapse = "\n")
  )
}


# ==============================================================================
# 5. LOCAL PLINK + EUR LD REFERENCE — V2 SAFE DETECTION
# ==============================================================================

has_triplet <- function(prefix) {
  all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
}

is_single_chr_prefix <- function(prefix) {
  grepl("\\.([1-9]|1[0-9]|2[0-2])$", prefix)
}

find_all_plink_prefixes <- function(root_dir) {

  beds <- list.files(
    root_dir,
    pattern = "\\.bed$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(beds) == 0) return(character(0))

  pref <- sub("\\.bed$", "", beds, ignore.case = TRUE)
  pref <- unique(pref[vapply(pref, has_triplet, logical(1))])
  pref
}

score_eur_prefix <- function(x) {
  s <- 0
  s <- s + 100 * grepl("EUR", x, ignore.case = TRUE)
  s <- s + 40 * grepl("1000G|1KG|Phase.?3", x, ignore.case = TRUE)
  s <- s + 20 * grepl("QC", x, ignore.case = TRUE)
  s
}

find_complete_split_eur_family <- function(prefixes) {

  # User's existing reference follows:
  #   <ROOT>.1, <ROOT>.2, ... <ROOT>.22
  split_pref <- prefixes[
    grepl("\\.([1-9]|1[0-9]|2[0-2])$", prefixes)
  ]

  if (length(split_pref) == 0) return(NULL)

  roots <- unique(
    sub(
      "\\.([1-9]|1[0-9]|2[0-2])$",
      "",
      split_pref
    )
  )

  good <- list()

  for (root in roots) {

    chr_pref <- paste0(root, ".", 1:22)

    if (all(vapply(chr_pref, has_triplet, logical(1)))) {

      good[[length(good) + 1]] <- list(
        root = root,
        chr_prefix = chr_pref,
        score = score_eur_prefix(root)
      )
    }
  }

  if (length(good) == 0) return(NULL)

  sc <- vapply(good, function(z) z$score, numeric(1))
  good[[which.max(sc)]]
}

find_whole_genome_eur_prefix <- function(prefixes) {

  # A whole-genome prefix must NOT end in .1 ... .22.
  whole <- prefixes[!is_single_chr_prefix(prefixes)]

  if (length(whole) == 0) return(NA_character_)

  score <- vapply(whole, score_eur_prefix, numeric(1))

  # Prefer the conventional OpenGWAS prefix exactly named EUR.
  score <- score + 200 * (tolower(basename(whole)) == "eur")

  whole[order(score, nchar(whole), decreasing = TRUE)][1]
}

find_plink_binary <- function() {

  roots <- unique(c(DATA_DIR, dirname(DATA_DIR)))

  hits <- character(0)

  for (r in roots) {

    if (dir.exists(r)) {

      hits <- c(
        hits,
        list.files(
          r,
          pattern = "^plink(\\.exe)?$",
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        )
      )
    }
  }

  if (length(hits) > 0) {
    return(hits[which.min(nchar(hits))])
  }

  if (!requireNamespace("genetics.binaRies", quietly = TRUE)) {

    try(
      remotes::install_github(
        "MRCIEU/genetics.binaRies",
        upgrade = "never"
      ),
      silent = TRUE
    )
  }

  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    return(genetics.binaRies::get_plink_binary())
  }

  NA_character_
}

# --------------------------
# Detect existing reference
# --------------------------

PLINK_PREFIXES <- find_all_plink_prefixes(DATA_DIR)

SPLIT_REF <- find_complete_split_eur_family(PLINK_PREFIXES)
WHOLE_REF <- find_whole_genome_eur_prefix(PLINK_PREFIXES)

LD_REF_MODE <- NA_character_
BFILE_WHOLE <- NA_character_
BFILE_CHR <- rep(NA_character_, 22)

# Prefer a complete 22-chromosome family if it exists.
if (!is.null(SPLIT_REF)) {

  LD_REF_MODE <- "split_22chr"
  BFILE_CHR <- SPLIT_REF$chr_prefix

} else if (
  !is.na(WHOLE_REF) &&
  has_triplet(WHOLE_REF)
) {

  LD_REF_MODE <- "whole_genome"
  BFILE_WHOLE <- WHOLE_REF
}

# --------------------------
# Download fallback only if no valid local reference exists
# --------------------------

if (is.na(LD_REF_MODE) && ALLOW_REFERENCE_DOWNLOAD) {

  log_msg("No complete EUR PLINK reference detected locally.")
  log_msg("Downloading official OpenGWAS/1000G EUR LD reference...")

  TGZ <- file.path(DIR_REF, "1kg.v3.tgz")

  if (!file.exists(TGZ) || file.info(TGZ)$size < 1e6) {

    urls <- c(
      "https://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz",
      "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
    )

    success <- FALSE

    for (u in urls) {

      z <- try(
        utils::download.file(
          u,
          TGZ,
          mode = "wb",
          quiet = FALSE
        ),
        silent = TRUE
      )

      if (
        !inherits(z, "try-error") &&
        file.exists(TGZ) &&
        file.info(TGZ)$size > 1e6
      ) {
        success <- TRUE
        break
      }
    }

    if (!success) {
      stop(
        "Could not download the EUR LD reference automatically.\n",
        "You already appear to have chromosome-specific 1000G files in most ",
        "setups; if so, verify that all 1-22 BED/BIM/FAM triplets exist."
      )
    }
  }

  EXTRACT_DIR <- file.path(DIR_REF, "1kg.v3")
  dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

  try(
    utils::untar(
      TGZ,
      exdir = EXTRACT_DIR
    ),
    silent = TRUE
  )

  new_pref <- find_all_plink_prefixes(EXTRACT_DIR)
  new_split <- find_complete_split_eur_family(new_pref)
  new_whole <- find_whole_genome_eur_prefix(new_pref)

  if (!is.null(new_split)) {

    LD_REF_MODE <- "split_22chr"
    BFILE_CHR <- new_split$chr_prefix

  } else if (
    !is.na(new_whole) &&
    has_triplet(new_whole)
  ) {

    LD_REF_MODE <- "whole_genome"
    BFILE_WHOLE <- new_whole
  }
}

if (is.na(LD_REF_MODE)) {
  stop(
    "No valid EUR LD reference was found.\n",
    "V2 requires either:\n",
    "  1) one whole-genome EUR.bed/EUR.bim/EUR.fam, OR\n",
    "  2) a COMPLETE chromosome family *.1 ... *.22, each with BED/BIM/FAM."
  )
}

# HARD PROTECTION against the exact V1 bug.
if (
  identical(LD_REF_MODE, "whole_genome") &&
  is_single_chr_prefix(BFILE_WHOLE)
) {
  stop(
    "SAFETY STOP: a chromosome-specific prefix was selected as a whole-genome ",
    "reference: ", BFILE_WHOLE,
    "\nThis is the exact V1 error and V2 refuses to continue."
  )
}

if (identical(LD_REF_MODE, "split_22chr")) {

  if (length(BFILE_CHR) != 22L) {
    stop("SAFETY STOP: split reference does not contain exactly 22 prefixes.")
  }

  ok22 <- vapply(BFILE_CHR, has_triplet, logical(1))

  if (!all(ok22)) {
    stop(
      "SAFETY STOP: missing BED/BIM/FAM for chromosomes: ",
      paste(which(!ok22), collapse = ", ")
    )
  }
}

PLINK_BIN <- find_plink_binary()

if (is.na(PLINK_BIN) || !file.exists(PLINK_BIN)) {
  stop("PLINK binary could not be located or installed.")
}

log_msg("PLINK = ", PLINK_BIN)
log_msg("LD reference mode = ", LD_REF_MODE)

if (LD_REF_MODE == "whole_genome") {

  log_msg("Whole-genome EUR bfile = ", BFILE_WHOLE)

} else {

  log_msg("22-chromosome EUR reference root = ", SPLIT_REF$root)

  for (cc in 1:22) {
    log_msg("  chr", cc, " = ", BFILE_CHR[cc])
  }
}

# --------------------------
# Build rsID -> reference chromosome index for split reference
# --------------------------

REF_SNP_INDEX <- NULL

if (LD_REF_MODE == "split_22chr") {

  log_msg("Building rsID-to-chromosome index from 22 BIM files...")

  ref_list <- vector("list", 22)

  for (cc in 1:22) {

    bim_file <- paste0(BFILE_CHR[cc], ".bim")

    b <- data.table::fread(
      bim_file,
      header = FALSE,
      select = c(1, 2, 4),
      col.names = c("CHR", "SNP", "BP"),
      showProgress = FALSE
    )

    b[, CHR := as.integer(CHR)]
    b[, BP := as.integer(BP)]

    # Defensive check: chromosome file should actually represent cc.
    chr_values <- unique(b$CHR[is.finite(b$CHR)])

    if (
      length(chr_values) > 0 &&
      !all(chr_values == cc)
    ) {
      warning(
        "Reference file ", bim_file,
        " contains chromosome labels other than ", cc,
        ". The BIM chromosome labels will be used as authoritative."
      )
    }

    ref_list[[cc]] <- b[, .(SNP, REF_CHR = CHR, REF_BP = BP)]
  }

  REF_SNP_INDEX <- rbindlist(ref_list, use.names = TRUE)

  # rsIDs should be unique across autosomes. Keep first occurrence defensively.
  REF_SNP_INDEX <- unique(
    REF_SNP_INDEX,
    by = "SNP"
  )

  log_msg(
    "Reference index unique autosomal SNPs = ",
    nrow(REF_SNP_INDEX)
  )

  # Reference QC
  ref_qc <- rbindlist(
    lapply(
      1:22,
      function(cc) {
        data.table(
          chromosome = cc,
          bfile_prefix = BFILE_CHR[cc],
          bed_exists = file.exists(paste0(BFILE_CHR[cc], ".bed")),
          bim_exists = file.exists(paste0(BFILE_CHR[cc], ".bim")),
          fam_exists = file.exists(paste0(BFILE_CHR[cc], ".fam")),
          n_bim_snps = sum(REF_SNP_INDEX$REF_CHR == cc, na.rm = TRUE)
        )
      }
    )
  )

  fwrite(
    ref_qc,
    file.path(
      DIR_REP,
      "STEP3B_V2_LD_reference_QC.csv"
    )
  )

} else {

  bim_file <- paste0(BFILE_WHOLE, ".bim")

  b <- data.table::fread(
    bim_file,
    header = FALSE,
    select = c(1, 2, 4),
    col.names = c("CHR", "SNP", "BP"),
    showProgress = FALSE
  )

  b[, CHR := as.integer(CHR)]
  b[, BP := as.integer(BP)]

  REF_SNP_INDEX <- unique(
    b[
      CHR %in% 1:22,
      .(
        SNP,
        REF_CHR = CHR,
        REF_BP = BP
      )
    ],
    by = "SNP"
  )

  fwrite(
    data.table(
      mode = "whole_genome",
      bfile_prefix = BFILE_WHOLE,
      n_autosomal_snps = nrow(REF_SNP_INDEX)
    ),
    file.path(
      DIR_REP,
      "STEP3B_V2_LD_reference_QC.csv"
    )
  )
}

if (
  is.null(REF_SNP_INDEX) ||
  nrow(REF_SNP_INDEX) == 0
) {
  stop("SAFETY STOP: the EUR reference SNP index is empty.")
}


# ==============================================================================
# 6. DUCKDB CONNECTION
# ==============================================================================

DB_FILE <- file.path(DIR_CACHE, "STEP3B_temp.duckdb")

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_FILE,
  read_only = FALSE
)

dbExecute(con, "SET threads=4;")

dbExecute(
  con,
  sprintf(
    "SET temp_directory='%s';",
    gsub("'", "''", normalizePath(DIR_CACHE, winslash = "/", mustWork = FALSE))
  )
)


# ==============================================================================
# 7. GWAS COLUMN MAPPING
# ==============================================================================

sql_path <- function(x) {
  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(x, winslash = "/", mustWork = TRUE)
    ),
    "'"
  )
}

qid <- function(x) {
  paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
}

get_schema <- function(path) {
  dbGetQuery(
    con,
    sprintf(
      "DESCRIBE SELECT * FROM read_csv_auto(%s, header=true, sample_size=100000);",
      sql_path(path)
    )
  )
}

pick_col <- function(cols, candidates, required = TRUE) {
  low <- tolower(cols)

  for (cand in candidates) {
    i <- which(low == tolower(cand))
    if (length(i) > 0) return(cols[i[1]])
  }

  if (required) {
    stop(
      "Required column not found. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable:\n",
      paste(cols, collapse = " | ")
    )
  }

  NA_character_
}

is_avail <- function(x) {
  length(x) == 1 && !is.na(x) && nzchar(x)
}

make_spec <- function(trait, path) {

  sch <- get_schema(path)
  cols <- sch$column_name

  SNP <- pick_col(
    cols,
    c("rs_id", "rsID", "rsids", "SNP", "MarkerName")
  )

  A1 <- pick_col(
    cols,
    c("effect_allele", "A1", "alt", "allele1")
  )

  A2 <- pick_col(
    cols,
    c("other_allele", "A2", "ref", "allele2")
  )

  BETA <- pick_col(
    cols,
    c("beta", "A1_beta", "b", "effect", "estimate")
  )

  SE <- pick_col(
    cols,
    c("standard_error", "A1_se", "sebeta", "se", "stderr")
  )

  P <- pick_col(
    cols,
    c("p_value", "A1_p", "pval", "p", "pvalue")
  )

  EAF <- pick_col(
    cols,
    c(
      "effect_allele_frequency",
      "A1_freq",
      "af_alt",
      "Freq1.Hapmap",
      "Freq1",
      "eaf",
      "af"
    )
  )

  CHR <- pick_col(
    cols,
    c("chromosome", "#chrom", "chrom", "chr"),
    required = FALSE
  )

  BP <- pick_col(
    cols,
    c("base_pair_location", "pos_b37", "pos", "position", "bp"),
    required = FALSE
  )

  N <- pick_col(
    cols,
    c("n", "N", "N_total", "n_total", "samplesize"),
    required = FALSE
  )

  NCASE <- pick_col(
    cols,
    c("N_case", "n_case", "ncase"),
    required = FALSE
  )

  NTOTAL <- pick_col(
    cols,
    c("N_total", "n_total", "n"),
    required = FALSE
  )

  list(
    trait = trait,
    path = path,
    cols = cols,
    SNP = SNP,
    A1 = A1,
    A2 = A2,
    BETA = BETA,
    SE = SE,
    P = P,
    EAF = EAF,
    CHR = CHR,
    BP = BP,
    N = N,
    NCASE = NCASE,
    NTOTAL = NTOTAL
  )
}

SPECS <- list(
  AF = make_spec("AF", AF_PATH),
  HFpEF = make_spec("HFpEF", HF_PATH),
  BMI = make_spec("BMI", BMI_PATH),
  OSA = make_spec("OSA", OSA_PATH)
)

header_audit <- rbindlist(
  lapply(
    SPECS,
    function(s) {
      data.table(
        trait = s$trait,
        path = s$path,
        SNP = s$SNP,
        A1 = s$A1,
        A2 = s$A2,
        BETA = s$BETA,
        SE = s$SE,
        P = s$P,
        EAF = s$EAF,
        CHR = s$CHR,
        BP = s$BP,
        N = s$N,
        NCASE = s$NCASE,
        NTOTAL = s$NTOTAL
      )
    }
  ),
  fill = TRUE
)

fwrite(
  header_audit,
  file.path(DIR_REP, "STEP3B_header_mapping.csv")
)


# ==============================================================================
# 8. EFFECTIVE SAMPLE SIZE EXPRESSIONS
# ==============================================================================

# AF study counts used in this project
AF_CASE <- 228926
AF_CTRL <- 1611415
AF_TOTAL <- AF_CASE + AF_CTRL
AF_CASE_FRAC <- AF_CASE / AF_TOTAL
AF_NEFF_FACTOR <- 4 * AF_CASE_FRAC * (1 - AF_CASE_FRAC)
AF_NEFF_CONST <- 4 / (1 / AF_CASE + 1 / AF_CTRL)

# OSA FinnGen R9 counts
OSA_CASE <- 38998
OSA_CTRL <- 336659
OSA_NEFF_CONST <- 4 / (1 / OSA_CASE + 1 / OSA_CTRL)

make_n_expr <- function(spec) {

  tr <- spec$trait

  if (tr == "AF") {
    if (is_avail(spec$N)) {
      return(
        sprintf(
          "(TRY_CAST(%s AS DOUBLE) * %.15f)",
          qid(spec$N),
          AF_NEFF_FACTOR
        )
      )
    } else {
      return(sprintf("%.15f", AF_NEFF_CONST))
    }
  }

  if (tr == "HFpEF") {
    if (is_avail(spec$NCASE) && is_avail(spec$NTOTAL)) {
      return(
        sprintf(
          "(4.0 / (1.0 / TRY_CAST(%s AS DOUBLE) + 1.0 / (TRY_CAST(%s AS DOUBLE) - TRY_CAST(%s AS DOUBLE))))",
          qid(spec$NCASE),
          qid(spec$NTOTAL),
          qid(spec$NCASE)
        )
      )
    } else if (is_avail(spec$N)) {
      return(sprintf("TRY_CAST(%s AS DOUBLE)", qid(spec$N)))
    } else {
      return("NULL::DOUBLE")
    }
  }

  if (tr == "OSA") {
    return(sprintf("%.15f", OSA_NEFF_CONST))
  }

  if (tr == "BMI") {
    if (is_avail(spec$N)) {
      return(sprintf("TRY_CAST(%s AS DOUBLE)", qid(spec$N)))
    } else {
      return("NULL::DOUBLE")
    }
  }

  "NULL::DOUBLE"
}


# ==============================================================================
# 9. STANDARDISED GWAS QUERY
# ==============================================================================

standardised_base_sql <- function(spec) {

  snp_expr <- sprintf(
    "regexp_extract(CAST(%s AS VARCHAR), 'rs[0-9]+')",
    qid(spec$SNP)
  )

  chr_expr <- if (is_avail(spec$CHR)) {
    sprintf("TRY_CAST(%s AS INTEGER)", qid(spec$CHR))
  } else {
    "NULL::INTEGER"
  }

  bp_expr <- if (is_avail(spec$BP)) {
    sprintf("TRY_CAST(%s AS BIGINT)", qid(spec$BP))
  } else {
    "NULL::BIGINT"
  }

  n_expr <- make_n_expr(spec)

  sprintf(
"
SELECT
  %s AS SNP,
  upper(CAST(%s AS VARCHAR)) AS A1,
  upper(CAST(%s AS VARCHAR)) AS A2,
  TRY_CAST(%s AS DOUBLE) AS BETA,
  TRY_CAST(%s AS DOUBLE) AS SE,
  TRY_CAST(%s AS DOUBLE) AS P,
  TRY_CAST(%s AS DOUBLE) AS EAF,
  %s AS N,
  %s AS CHR,
  %s AS BP
FROM read_csv_auto(%s, header=true, sample_size=100000)
",
    snp_expr,
    qid(spec$A1),
    qid(spec$A2),
    qid(spec$BETA),
    qid(spec$SE),
    qid(spec$P),
    qid(spec$EAF),
    n_expr,
    chr_expr,
    bp_expr,
    sql_path(spec$path)
  )
}

query_trait <- function(spec, p_max = NULL, snps = NULL) {

  base <- standardised_base_sql(spec)

  temp_name <- NULL
  join_sql <- ""

  if (!is.null(snps)) {

    snps <- unique(as.character(snps))
    snps <- snps[grepl("^rs[0-9]+$", snps)]

    if (length(snps) == 0) return(data.table())

    temp_name <- paste0(
      "wanted_",
      gsub("[^A-Za-z0-9]", "_", spec$trait),
      "_",
      as.integer(runif(1, 1, 1e8))
    )

    dbWriteTable(
      con,
      temp_name,
      data.frame(SNP = snps),
      overwrite = TRUE,
      temporary = TRUE
    )

    join_sql <- sprintf(
      "INNER JOIN %s w ON x.SNP = w.SNP",
      qid(temp_name)
    )
  }

  p_filter <- if (!is.null(p_max)) {
    sprintf("AND x.P <= %.17g", p_max)
  } else {
    ""
  }

  sql <- sprintf(
"
WITH x AS (
  %s
),
q AS (
  SELECT x.*
  FROM x
  %s
  WHERE
    x.SNP IS NOT NULL
    AND regexp_matches(x.SNP, '^rs[0-9]+$')
    AND x.A1 IN ('A','C','G','T')
    AND x.A2 IN ('A','C','G','T')
    AND x.A1 <> x.A2
    AND x.BETA IS NOT NULL
    AND x.SE IS NOT NULL
    AND x.SE > 0
    AND x.P IS NOT NULL
    AND x.P > 0
    AND x.P <= 1
    AND x.EAF IS NOT NULL
    AND x.EAF > 0
    AND x.EAF < 1
    %s
),
d AS (
  SELECT *,
         row_number() OVER (
           PARTITION BY SNP
           ORDER BY P ASC
         ) AS rn
  FROM q
)
SELECT SNP, A1, A2, BETA, SE, P, EAF, N, CHR, BP
FROM d
WHERE rn = 1
",
    base,
    join_sql,
    p_filter
  )

  out <- as.data.table(dbGetQuery(con, sql))

  if (!is.null(temp_name)) {
    try(dbRemoveTable(con, temp_name), silent = TRUE)
  }

  out
}


# ==============================================================================
# 10. LOAD ONLY POTENTIAL INSTRUMENTS (MAX P = 5e-6)
# ==============================================================================

MAX_P <- max(P_TIERS)

candidate_tables <- list()

for (tr in names(SPECS)) {

  log_msg("Reading potential instruments for ", tr, " (P <= ", MAX_P, ") ...")

  cand <- query_trait(
    SPECS[[tr]],
    p_max = MAX_P
  )

  cand[, Fstat := (BETA / SE)^2]

  setorder(cand, P)

  candidate_tables[[tr]] <- cand

  fwrite(
    cand,
    file.path(
      DIR_CACHE,
      paste0("STEP3B_", tr, "_candidate_Ple5e-6.csv.gz")
    )
  )

  log_msg(
    "  ", tr,
    ": candidate SNPs = ", nrow(cand),
    "; F>=10 = ", sum(cand$Fstat >= F_MIN, na.rm = TRUE)
  )
}


# ==============================================================================
# 11. LOCAL LD CLUMPING — V2 CHROMOSOME-SAFE IMPLEMENTATION
# ==============================================================================

CLUMP_QC_ACCUM <- list()

local_clump <- function(dt, trait, clump_p = 1) {

  if (nrow(dt) == 0) return(dt[0])

  x <- copy(dt)

  x <- x[
    !is.na(SNP) &
    grepl("^rs[0-9]+$", SNP) &
    is.finite(P)
  ]

  x <- unique(
    x,
    by = "SNP"
  )

  if (nrow(x) == 0) return(x)

  n_input <- nrow(x)

  # --------------------------------------------------------------
  # Map every candidate rsID to its chromosome in the actual
  # 1000G EUR reference. This is authoritative for clumping.
  # --------------------------------------------------------------

  map <- REF_SNP_INDEX[
    SNP %in% x$SNP
  ]

  xmap <- merge(
    x,
    map,
    by = "SNP",
    all = FALSE
  )

  n_mapped <- nrow(xmap)
  match_rate <- n_mapped / n_input

  missing_snps <- x[
    !SNP %in% xmap$SNP
  ]

  if (nrow(missing_snps) > 0) {

    fwrite(
      missing_snps,
      file.path(
        DIR_INST,
        paste0(
          gsub("[^A-Za-z0-9_.-]", "_", trait),
          "__not_in_1000G_EUR_reference.csv"
        )
      )
    )
  }

  log_msg(
    "CLUMP ", trait,
    ": input=", n_input,
    "; mapped_to_ref=", n_mapped,
    "; match_rate=", sprintf("%.3f", match_rate)
  )

  if (
    n_input >= 10 &&
    match_rate < MIN_REF_MATCH_RATE
  ) {
    stop(
      "SAFETY STOP for ", trait, ": only ",
      sprintf("%.1f%%", 100 * match_rate),
      " of candidate SNPs were found in the EUR LD reference.\n",
      "This suggests an rsID/build/reference mismatch. ",
      "The script will NOT continue with a biased subset."
    )
  }

  if (n_mapped == 0) {
    stop(
      "No candidate SNPs for ", trait,
      " are present in the EUR LD reference."
    )
  }

  kept_all <- character(0)
  per_chr_qc <- list()

  # --------------------------------------------------------------
  # Whole-genome reference:
  # one ordinary local PLINK clump is valid.
  # --------------------------------------------------------------

  if (LD_REF_MODE == "whole_genome") {

    if (is_single_chr_prefix(BFILE_WHOLE)) {
      stop(
        "SAFETY STOP: single-chromosome reference reached ",
        "whole-genome clumping code."
      )
    }

    inp <- unique(
      data.table(
        rsid = xmap$SNP,
        pval = xmap$P,
        id = trait
      ),
      by = "rsid"
    )

    out <- try(
      ieugwasr::ld_clump(
        dat = as.data.frame(inp),
        clump_kb = CLUMP_KB,
        clump_r2 = CLUMP_R2,
        clump_p = clump_p,
        pop = "EUR",
        bfile = BFILE_WHOLE,
        plink_bin = PLINK_BIN
      ),
      silent = TRUE
    )

    if (inherits(out, "try-error")) {
      stop(
        "Whole-genome local LD clumping failed for ",
        trait,
        ".\nError:\n",
        as.character(out)
      )
    }

    snp_col <- if ("rsid" %in% names(out)) {
      "rsid"
    } else if ("variant" %in% names(out)) {
      "variant"
    } else {
      stop("Could not identify SNP column in ld_clump output.")
    }

    kept_all <- unique(as.character(out[[snp_col]]))

    per_chr_qc[[1]] <- data.table(
      trait = trait,
      chromosome = NA_integer_,
      n_input = n_input,
      n_mapped = n_mapped,
      n_kept = length(kept_all),
      reference_mode = "whole_genome",
      bfile_prefix = BFILE_WHOLE
    )

  } else {

    # ------------------------------------------------------------
    # SPLIT 22-CHR REFERENCE:
    # clump chromosome-by-chromosome using the matching prefix.
    # Cross-chromosome LD does not need to be modeled.
    # ------------------------------------------------------------

    chr_present <- sort(
      unique(
        xmap$REF_CHR[
          xmap$REF_CHR %in% 1:22
        ]
      )
    )

    if (length(chr_present) == 0) {
      stop(
        "Mapped SNPs for ", trait,
        " do not have valid autosomal chromosome labels."
      )
    }

    for (cc in chr_present) {

      dchr <- xmap[
        REF_CHR == cc
      ]

      if (nrow(dchr) == 0) next

      bfile_cc <- BFILE_CHR[cc]

      # HARD PROTECTION: exact chromosome prefix must match cc.
      if (!grepl(
        paste0("\\.", cc, "$"),
        bfile_cc
      )) {
        stop(
          "SAFETY STOP: chromosome ", cc,
          " is mapped to unexpected bfile prefix:\n",
          bfile_cc
        )
      }

      if (!has_triplet(bfile_cc)) {
        stop(
          "Missing BED/BIM/FAM for chromosome ",
          cc, ": ", bfile_cc
        )
      }

      inp <- unique(
        data.table(
          rsid = dchr$SNP,
          pval = dchr$P,
          id = paste0(trait, "_chr", cc)
        ),
        by = "rsid"
      )

      out <- try(
        ieugwasr::ld_clump(
          dat = as.data.frame(inp),
          clump_kb = CLUMP_KB,
          clump_r2 = CLUMP_R2,
          clump_p = clump_p,
          pop = "EUR",
          bfile = bfile_cc,
          plink_bin = PLINK_BIN
        ),
        silent = TRUE
      )

      if (inherits(out, "try-error")) {
        stop(
          "Chromosome-specific LD clumping failed.\n",
          "Trait: ", trait,
          "\nChromosome: ", cc,
          "\nBFILE: ", bfile_cc,
          "\nError:\n",
          as.character(out)
        )
      }

      snp_col <- if ("rsid" %in% names(out)) {
        "rsid"
      } else if ("variant" %in% names(out)) {
        "variant"
      } else {
        stop(
          "Could not identify SNP column in ld_clump output ",
          "for chromosome ", cc, "."
        )
      }

      kept_chr <- unique(
        as.character(
          out[[snp_col]]
        )
      )

      kept_all <- c(
        kept_all,
        kept_chr
      )

      per_chr_qc[[length(per_chr_qc) + 1]] <- data.table(
        trait = trait,
        chromosome = cc,
        n_input = nrow(dchr),
        n_mapped = nrow(dchr),
        n_kept = length(kept_chr),
        reference_mode = "split_22chr",
        bfile_prefix = bfile_cc
      )

      log_msg(
        "  chr", cc,
        ": input=", nrow(dchr),
        "; kept=", length(kept_chr)
      )
    }
  }

  kept_all <- unique(kept_all)

  ans <- x[
    SNP %in% kept_all
  ]

  setorder(ans, P)

  # --------------------------------------------------------------
  # QC summaries
  # --------------------------------------------------------------

  qc_chr <- rbindlist(
    per_chr_qc,
    fill = TRUE
  )

  qc_summary <- data.table(
    trait = trait,
    clump_p = clump_p,
    clump_r2 = CLUMP_R2,
    clump_kb = CLUMP_KB,
    reference_mode = LD_REF_MODE,
    n_input = n_input,
    n_mapped_reference = n_mapped,
    reference_match_rate = match_rate,
    n_chromosomes_with_candidates =
      length(unique(xmap$REF_CHR)),
    n_kept = nrow(ans)
  )

  fwrite(
    qc_chr,
    file.path(
      DIR_INST,
      paste0(
        gsub("[^A-Za-z0-9_.-]", "_", trait),
        "__clump_by_chromosome_QC.csv"
      )
    )
  )

  fwrite(
    qc_summary,
    file.path(
      DIR_INST,
      paste0(
        gsub("[^A-Za-z0-9_.-]", "_", trait),
        "__clump_summary_QC.csv"
      )
    )
  )

  log_msg(
    "CLUMP COMPLETE ", trait,
    ": mapped=", n_mapped,
    "; kept=", nrow(ans),
    "; chromosomes=", length(unique(xmap$REF_CHR))
  )

  # Very important sanity warning: if input candidates span many chromosomes
  # but clumped SNPs all come from one chromosome, flag it.
  if (
    length(unique(xmap$REF_CHR)) >= 5 &&
    nrow(ans) >= 5
  ) {

    kept_chr_n <- unique(
      xmap[
        SNP %in% ans$SNP,
        REF_CHR
      ]
    )

    if (length(kept_chr_n) == 1) {
      stop(
        "SAFETY STOP: candidates span ",
        length(unique(xmap$REF_CHR)),
        " chromosomes but every retained SNP is from chromosome ",
        kept_chr_n,
        ". This is suspicious for an LD-reference/clumping error."
      )
    }
  }

  ans
}


# ==============================================================================
# 12. SCIENTIFIC AUTOMATIC THRESHOLD RELAXATION
# ==============================================================================

selection_qc <- list()
instrument_sets <- list()

for (tr in names(candidate_tables)) {

  cand <- candidate_tables[[tr]]

  tier_sets <- list()

  for (tier_nm in names(P_TIERS)) {

    pthr <- P_TIERS[[tier_nm]]

    x <- cand[
      P <= pthr &
      is.finite(Fstat) &
      Fstat >= F_MIN
    ]

    n_pre <- nrow(x)

    xc <- local_clump(
      x,
      trait = paste0(tr, "_", tier_nm),
      clump_p = pthr
    )

    n_post <- nrow(xc)

    tier_sets[[tier_nm]] <- xc

    selection_qc[[length(selection_qc) + 1]] <- data.table(
      trait = tr,
      tier = tier_nm,
      p_threshold = pthr,
      F_threshold = F_MIN,
      clump_r2 = CLUMP_R2,
      clump_kb = CLUMP_KB,
      n_before_clump = n_pre,
      n_after_clump = n_post,
      reaches_target_n3 = n_post >= MIN_IV_TARGET,
      analysis_class = ifelse(
        tier_nm == "primary_GWS",
        "PRIMARY",
        "EXPLORATORY"
      )
    )

    fwrite(
      xc,
      file.path(
        DIR_INST,
        paste0(
          tr, "__", tier_nm,
          "__clumped_instruments.csv"
        )
      )
    )
  }

  strict <- tier_sets[["primary_GWS"]]

  # Selected exploratory tier:
  # - if strict already has >=3 IVs, selected = strict
  # - otherwise first relaxed tier reaching >=3
  # - if none reaches >=3, choose tier with maximum available IVs
  if (nrow(strict) >= MIN_IV_TARGET) {

    selected_nm <- "primary_GWS"

  } else {

    relaxed_names <- setdiff(names(P_TIERS), "primary_GWS")
    reached <- relaxed_names[
      vapply(
        tier_sets[relaxed_names],
        nrow,
        integer(1)
      ) >= MIN_IV_TARGET
    ]

    if (length(reached) > 0) {

      selected_nm <- reached[1]

    } else {

      ns <- vapply(tier_sets, nrow, integer(1))
      selected_nm <- names(ns)[which.max(ns)]
    }
  }

  instrument_sets[[tr]] <- list(
    strict = strict,
    selected = tier_sets[[selected_nm]],
    selected_tier = selected_nm,
    selected_p = P_TIERS[[selected_nm]],
    all_tiers = tier_sets
  )

  log_msg(
    "Instrument selection ", tr,
    ": strict=", nrow(strict),
    "; selected tier=", selected_nm,
    "; selected n=", nrow(tier_sets[[selected_nm]])
  )
}

selection_qc_dt <- rbindlist(selection_qc, fill = TRUE)

fwrite(
  selection_qc_dt,
  file.path(
    DIR_REP,
    "STEP3B_instrument_selection_QC.csv"
  )
)


# ==============================================================================
# 13. QUERY ALL REQUIRED SNP ASSOCIATIONS ONCE PER TRAIT
# ==============================================================================

all_needed_snps <- unique(
  unlist(
    lapply(
      instrument_sets,
      function(z) unique(c(z$strict$SNP, z$selected$SNP))
    )
  )
)

log_msg("Unique SNPs requiring cross-trait lookup = ", length(all_needed_snps))

assoc_tables <- list()

for (tr in names(SPECS)) {

  log_msg("Extracting cross-trait associations for ", tr, "...")

  aa <- query_trait(
    SPECS[[tr]],
    snps = all_needed_snps
  )

  assoc_tables[[tr]] <- aa

  fwrite(
    aa,
    file.path(
      DIR_CACHE,
      paste0("STEP3B_", tr, "_needed_SNP_associations.csv.gz")
    )
  )

  log_msg("  found = ", nrow(aa), " SNPs")
}


# ==============================================================================
# 14. TWOSAMPLEMR FORMAT HELPERS
# ==============================================================================

to_tsmr_exposure <- function(dt, trait) {

  if (nrow(dt) == 0) return(data.frame())

  d <- data.frame(
    SNP = dt$SNP,
    beta = dt$BETA,
    se = dt$SE,
    effect_allele = dt$A1,
    other_allele = dt$A2,
    eaf = dt$EAF,
    pval = dt$P,
    samplesize = dt$N,
    Phenotype = trait,
    id = trait,
    stringsAsFactors = FALSE
  )

  x <- TwoSampleMR::format_data(
    d,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    id_col = "id"
  )

  x
}

to_tsmr_outcome <- function(dt, trait) {

  if (nrow(dt) == 0) return(data.frame())

  d <- data.frame(
    SNP = dt$SNP,
    beta = dt$BETA,
    se = dt$SE,
    effect_allele = dt$A1,
    other_allele = dt$A2,
    eaf = dt$EAF,
    pval = dt$P,
    samplesize = dt$N,
    Phenotype = trait,
    id = trait,
    stringsAsFactors = FALSE
  )

  x <- TwoSampleMR::format_data(
    d,
    type = "outcome",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    id_col = "id"
  )

  x
}


# ==============================================================================
# 15. MR-PRESSO HELPER
# ==============================================================================

run_presso <- function(hdat, tag) {

  out_txt <- file.path(
    DIR_SENS,
    paste0(tag, "__MRPRESSO.txt")
  )

  if (!MRPRESSO_AVAILABLE) {
    writeLines(
      "MR-PRESSO package unavailable; analysis skipped.",
      out_txt
    )
    return(
      data.table(
        tag = tag,
        MRPRESSO_ran = FALSE,
        MRPRESSO_note = "package unavailable"
      )
    )
  }

  if (nrow(hdat) < 4) {
    writeLines(
      "MR-PRESSO not run because <4 harmonised IVs.",
      out_txt
    )
    return(
      data.table(
        tag = tag,
        MRPRESSO_ran = FALSE,
        MRPRESSO_note = "<4 IVs"
      )
    )
  }

  pr <- try(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = as.data.frame(hdat),
      NbDistribution = MRPRESSO_NB_DISTRIBUTION,
      SignifThreshold = 0.05
    ),
    silent = TRUE
  )

  capture.output(
    pr,
    file = out_txt
  )

  if (inherits(pr, "try-error")) {
    return(
      data.table(
        tag = tag,
        MRPRESSO_ran = FALSE,
        MRPRESSO_note = paste0("error: ", as.character(pr))
      )
    )
  }

  # Extract global P if structure is available; otherwise leave NA.
  gp <- NA_real_

  gp_try <- try(
    pr[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
    silent = TRUE
  )

  if (!inherits(gp_try, "try-error") && length(gp_try) > 0) {
    gp_chr <- as.character(gp_try[1])
    gp_num <- suppressWarnings(as.numeric(gsub("[^0-9.eE-]", "", gp_chr)))
    if (is.finite(gp_num)) gp <- gp_num
  }

  data.table(
    tag = tag,
    MRPRESSO_ran = TRUE,
    MRPRESSO_global_p = gp,
    MRPRESSO_note = "see txt output"
  )
}


# ==============================================================================
# 16. RUN ONE UVMR ANALYSIS
# ==============================================================================

run_uvmr_one <- function(
  exposure,
  outcome,
  inst_dt,
  tier_name,
  p_threshold,
  analysis_class
) {

  tag <- paste0(
    exposure, "__to__", outcome,
    "__", tier_name
  )

  log_msg(
    "UVMR: ", exposure, " -> ", outcome,
    " | ", tier_name,
    " | instruments before harmonisation = ", nrow(inst_dt)
  )

  if (nrow(inst_dt) == 0) {

    return(
      list(
        summary = data.table(
          exposure = exposure,
          outcome = outcome,
          tier = tier_name,
          p_threshold = p_threshold,
          analysis_class = analysis_class,
          method = NA_character_,
          nsnp = 0L,
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = "NOT_ESTIMABLE_NO_IV"
        ),
        harmonised = data.table()
      )
    )
  }

  exp_assoc <- merge(
    inst_dt[, .(SNP)],
    assoc_tables[[exposure]],
    by = "SNP",
    all.x = FALSE,
    all.y = FALSE
  )

  out_assoc <- assoc_tables[[outcome]][
    SNP %in% exp_assoc$SNP
  ]

  if (nrow(exp_assoc) == 0 || nrow(out_assoc) == 0) {

    return(
      list(
        summary = data.table(
          exposure = exposure,
          outcome = outcome,
          tier = tier_name,
          p_threshold = p_threshold,
          analysis_class = analysis_class,
          method = NA_character_,
          nsnp = 0L,
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = "NOT_ESTIMABLE_NO_OVERLAP"
        ),
        harmonised = data.table()
      )
    )
  }

  exp_fmt <- to_tsmr_exposure(exp_assoc, exposure)
  out_fmt <- to_tsmr_outcome(out_assoc, outcome)

  h <- TwoSampleMR::harmonise_data(
    exp_fmt,
    out_fmt,
    action = HARMONISE_ACTION
  )

  h <- as.data.table(h)
  h <- h[mr_keep %in% TRUE]

  # Final strength gate
  h[, Fstat_exposure := (beta.exposure / se.exposure)^2]
  h <- h[
    is.finite(Fstat_exposure) &
    Fstat_exposure >= F_MIN
  ]

  fwrite(
    h,
    file.path(
      DIR_HARM,
      paste0(tag, "__harmonised.csv")
    )
  )

  n_iv <- nrow(h)

  if (n_iv == 0) {

    return(
      list(
        summary = data.table(
          exposure = exposure,
          outcome = outcome,
          tier = tier_name,
          p_threshold = p_threshold,
          analysis_class = analysis_class,
          method = NA_character_,
          nsnp = 0L,
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = "NOT_ESTIMABLE_AFTER_HARMONISATION"
        ),
        harmonised = h
      )
    )
  }

  if (n_iv == 1) {

    methods <- c("mr_wald_ratio")

  } else if (n_iv == 2) {

    methods <- c("mr_ivw_mre")

  } else {

    methods <- c(
      "mr_ivw_mre",
      "mr_weighted_median",
      "mr_egger_regression"
    )
  }

  mrres <- TwoSampleMR::mr(
    as.data.frame(h),
    method_list = methods
  )

  mrres <- as.data.table(mrres)

  mrres[, `:=`(
    exposure_trait = exposure,
    outcome_trait = outcome,
    tier = tier_name,
    p_threshold = p_threshold,
    analysis_class = analysis_class,
    status = "OK",
    mean_F = mean(h$Fstat_exposure, na.rm = TRUE),
    min_F = min(h$Fstat_exposure, na.rm = TRUE)
  )]

  fwrite(
    mrres,
    file.path(
      DIR_UVMR,
      paste0(tag, "__MR_results.csv")
    )
  )

  # Heterogeneity
  if (n_iv >= 2) {

    het <- try(
      TwoSampleMR::mr_heterogeneity(
        as.data.frame(h),
        method_list = c("mr_ivw", "mr_egger_regression")
      ),
      silent = TRUE
    )

    if (!inherits(het, "try-error")) {
      fwrite(
        as.data.table(het),
        file.path(
          DIR_SENS,
          paste0(tag, "__heterogeneity.csv")
        )
      )
    }
  }

  # Egger intercept
  if (n_iv >= 3) {

    pleio <- try(
      TwoSampleMR::mr_pleiotropy_test(
        as.data.frame(h)
      ),
      silent = TRUE
    )

    if (!inherits(pleio, "try-error")) {
      fwrite(
        as.data.table(pleio),
        file.path(
          DIR_SENS,
          paste0(tag, "__Egger_intercept.csv")
        )
      )
    }
  }

  # Leave-one-out
  if (n_iv >= 3) {

    loo <- try(
      TwoSampleMR::mr_leaveoneout(
        as.data.frame(h),
        method = TwoSampleMR::mr_ivw
      ),
      silent = TRUE
    )

    if (!inherits(loo, "try-error")) {
      fwrite(
        as.data.table(loo),
        file.path(
          DIR_SENS,
          paste0(tag, "__leave_one_out.csv")
        )
      )
    }
  }

  # MR-PRESSO
  presso_qc <- run_presso(h, tag)

  fwrite(
    presso_qc,
    file.path(
      DIR_SENS,
      paste0(tag, "__MRPRESSO_QC.csv")
    )
  )

  list(
    summary = mrres,
    harmonised = h
  )
}


# ==============================================================================
# 17. PRE-SPECIFIED UVMR PAIRS
# ==============================================================================

uvmr_pairs <- data.table(
  exposure = c(
    "AF",
    "HFpEF",
    "BMI",
    "BMI",
    "OSA",
    "OSA"
  ),
  outcome = c(
    "HFpEF",
    "AF",
    "HFpEF",
    "AF",
    "HFpEF",
    "AF"
  )
)

fwrite(
  uvmr_pairs,
  file.path(
    DIR_REP,
    "STEP3B_predefined_UVMR_pairs.csv"
  )
)


# ==============================================================================
# 18. RUN PRIMARY + AUTO-RELAXED UVMR
# ==============================================================================

uvmr_results <- list()
uvmr_harmonised <- list()

for (i in seq_len(nrow(uvmr_pairs))) {

  ex <- uvmr_pairs$exposure[i]
  oy <- uvmr_pairs$outcome[i]

  z <- instrument_sets[[ex]]

  # ALWAYS run primary strict analysis, even if only one IV.
  r1 <- run_uvmr_one(
    exposure = ex,
    outcome = oy,
    inst_dt = z$strict,
    tier_name = "primary_GWS",
    p_threshold = P_TIERS[["primary_GWS"]],
    analysis_class = "PRIMARY"
  )

  uvmr_results[[length(uvmr_results) + 1]] <- r1$summary
  uvmr_harmonised[[paste0(ex, "__", oy, "__PRIMARY")]] <- r1$harmonised

  # Automatic exploratory relaxed analysis ONLY if primary has <3 IVs
  # and selected tier is actually relaxed.
  if (
    nrow(z$strict) < MIN_IV_TARGET &&
    z$selected_tier != "primary_GWS"
  ) {

    r2 <- run_uvmr_one(
      exposure = ex,
      outcome = oy,
      inst_dt = z$selected,
      tier_name = z$selected_tier,
      p_threshold = z$selected_p,
      analysis_class = "EXPLORATORY_RELAXED"
    )

    uvmr_results[[length(uvmr_results) + 1]] <- r2$summary
    uvmr_harmonised[[paste0(ex, "__", oy, "__RELAXED")]] <- r2$harmonised
  }
}

uvmr_all <- rbindlist(uvmr_results, fill = TRUE)

fwrite(
  uvmr_all,
  file.path(
    DIR_UVMR,
    "STEP3B_UVMR_all_methods_all_tiers.csv"
  )
)


# ==============================================================================
# 19. IDENTIFY PRIMARY ESTIMATOR + MULTIPLE-TESTING CORRECTION
# ==============================================================================

is_primary_method <- function(method, nsnp) {

  if (is.na(method) || is.na(nsnp)) return(FALSE)

  if (nsnp == 1) {
    return(grepl("Wald", method, ignore.case = TRUE))
  }

  grepl(
    "multiplicative random effects",
    method,
    ignore.case = TRUE
  )
}

primary_rows <- uvmr_all[
  analysis_class == "PRIMARY" &
  status == "OK"
]

if (nrow(primary_rows) > 0) {

  primary_rows[
    ,
    primary_estimator := mapply(
      is_primary_method,
      method,
      nsnp
    )
  ]

  primary_rows <- primary_rows[
    primary_estimator == TRUE
  ]
}

# Add one row for pre-specified pairs that were not estimable.
pair_key <- paste(
  uvmr_pairs$exposure,
  uvmr_pairs$outcome,
  sep = " -> "
)

primary_present_key <- if (nrow(primary_rows) > 0) {
  paste(
    primary_rows$exposure_trait,
    primary_rows$outcome_trait,
    sep = " -> "
  )
} else {
  character(0)
}

missing_keys <- setdiff(pair_key, primary_present_key)

if (length(missing_keys) > 0) {

  missing_rows <- rbindlist(
    lapply(
      missing_keys,
      function(k) {
        sp <- strsplit(k, " -> ", fixed = TRUE)[[1]]
        data.table(
          exposure_trait = sp[1],
          outcome_trait = sp[2],
          method = NA_character_,
          nsnp = 0L,
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          analysis_class = "PRIMARY",
          tier = "primary_GWS",
          status = "NOT_ESTIMABLE_AT_P5e-8"
        )
      }
    ),
    fill = TRUE
  )

  primary_rows <- rbindlist(
    list(primary_rows, missing_rows),
    fill = TRUE
  )
}

# Fixed Bonferroni family of 6 planned primary tests
primary_rows[, bonferroni_threshold := 0.05 / UVMR_PLANNED_TESTS]

nonmiss <- which(is.finite(primary_rows$pval))

primary_rows[, BH_FDR_q := NA_real_]

if (length(nonmiss) > 0) {
  primary_rows$BH_FDR_q[nonmiss] <- p.adjust(
    primary_rows$pval[nonmiss],
    method = "BH",
    n = UVMR_PLANNED_TESTS
  )
}

primary_rows[
  ,
  bonferroni_significant :=
    is.finite(pval) &
    pval < bonferroni_threshold
]

primary_rows[
  ,
  nominal_significant :=
    is.finite(pval) &
    pval < 0.05
]

primary_rows[
  ,
  significance_class := fifelse(
    bonferroni_significant,
    "Bonferroni significant",
    fifelse(
      nominal_significant,
      "Nominal",
      "Not significant"
    )
  )
]

primary_rows[
  ,
  `:=`(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se
  )
]

fwrite(
  primary_rows,
  file.path(
    DIR_UVMR,
    "STEP3B_UVMR_PRIMARY_final.csv"
  )
)

# Exploratory relaxed primary estimator table
relaxed_rows <- uvmr_all[
  analysis_class == "EXPLORATORY_RELAXED" &
  status == "OK"
]

if (nrow(relaxed_rows) > 0) {

  relaxed_rows[
    ,
    primary_estimator := mapply(
      is_primary_method,
      method,
      nsnp
    )
  ]

  relaxed_rows <- relaxed_rows[
    primary_estimator == TRUE
  ]

  relaxed_rows[
    ,
    `:=`(
      ci_low = b - 1.96 * se,
      ci_high = b + 1.96 * se,
      significance_class = fifelse(
        pval < 0.05,
        "Nominal exploratory",
        "Not significant"
      )
    )
  ]

  fwrite(
    relaxed_rows,
    file.path(
      DIR_UVMR,
      "STEP3B_UVMR_EXPLORATORY_RELAXED_final.csv"
    )
  )
}


# ==============================================================================
# 20. MAIN UVMR FIGURES — REFERENCE ARTICLE FIGURE 5 STYLE
# ==============================================================================

COL_BONF <- "#C65353"
COL_NOM  <- "#3F78A8"
COL_NS   <- "#8B9098"

make_primary_forest <- function(
  result_dt,
  outcome_name,
  filename_stub
) {

  d <- copy(result_dt[outcome_trait == outcome_name])

  wanted_order <- if (outcome_name == "HFpEF") {
    c("AF", "BMI", "OSA")
  } else {
    c("HFpEF", "BMI", "OSA")
  }

  d[, exposure_trait := factor(
    exposure_trait,
    levels = rev(wanted_order)
  )]

  # Binary outcomes AF/HFpEF -> OR scale
  d[, `:=`(
    est_plot = exp(b),
    low_plot = exp(ci_low),
    high_plot = exp(ci_high)
  )]

  d[
    ,
    point_colour := fifelse(
      significance_class == "Bonferroni significant",
      COL_BONF,
      fifelse(
        significance_class == "Nominal",
        COL_NOM,
        COL_NS
      )
    )
  ]

  d[
    ,
    effect_text := fifelse(
      is.finite(est_plot),
      sprintf(
        "%.2f (%.2f–%.2f)",
        est_plot,
        low_plot,
        high_plot
      ),
      "Not estimable"
    )
  ]

  p <- ggplot(
    d,
    aes(
      y = exposure_trait,
      x = est_plot
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "#6E6E6E"
    ) +
    geom_errorbarh(
      aes(
        xmin = low_plot,
        xmax = high_plot,
        colour = point_colour
      ),
      height = 0,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    geom_point(
      aes(colour = point_colour),
      shape = 15,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_colour_identity() +
    scale_x_log10() +
    labs(
      x = paste0(
        "Odds ratio for ",
        outcome_name,
        " (95% CI)"
      ),
      y = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family = BASE_FAMILY
    ) +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(8, 18, 8, 8)
    )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(filename_stub, ".pdf")
    ),
    p,
    width = 6.4,
    height = 3.4,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(filename_stub, ".tiff")
    ),
    p,
    width = 6.4,
    height = 3.4,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  fwrite(
    d,
    file.path(
      DIR_FIG,
      paste0(filename_stub, "_source_data.csv")
    )
  )
}

make_primary_forest(
  primary_rows,
  outcome_name = "HFpEF",
  filename_stub = "Figure5A_UVMR_PRIMARY_HFpEF_outcome"
)

make_primary_forest(
  primary_rows,
  outcome_name = "AF",
  filename_stub = "Figure5B_UVMR_PRIMARY_AF_outcome"
)

writeLines(
  c(
    "Figure 5 UVMR legend",
    "",
    "Squares indicate primary MR point estimates and horizontal lines indicate 95% confidence intervals.",
    "Primary analysis uses genome-wide significant instruments (P < 5×10^-8), F >= 10, LD r2 < 0.001 within 10,000 kb.",
    "For >=2 IVs, the primary estimator is inverse-variance weighted multiplicative random effects; for a single IV, the Wald ratio is used.",
    paste0(
      "Muted red/pink: Bonferroni significant at P < ",
      signif(0.05 / UVMR_PLANNED_TESTS, 4),
      " (0.05/6 pre-specified primary tests)."
    ),
    "Blue: nominal P < 0.05 but not Bonferroni significant.",
    "Gray: non-significant.",
    "Any P-value-threshold-relaxed analyses are exploratory and are excluded from the primary Figure 5 panels."
  ),
  file.path(
    DIR_FIG,
    "Figure5_UVMR_PRIMARY_legend.txt"
  )
)


# ==============================================================================
# 21. EXPLORATORY RELAXED FOREST PLOT
# ==============================================================================

if (exists("relaxed_rows") && nrow(relaxed_rows) > 0) {

  d <- copy(relaxed_rows)

  d[, pair := paste0(exposure_trait, " → ", outcome_trait)]
  d[, pair := factor(pair, levels = rev(unique(pair)))]

  d[, `:=`(
    est_plot = exp(b),
    low_plot = exp(ci_low),
    high_plot = exp(ci_high),
    col = fifelse(
      pval < 0.05,
      COL_NOM,
      COL_NS
    )
  )]

  p <- ggplot(
    d,
    aes(
      y = pair,
      x = est_plot
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "#6E6E6E"
    ) +
    geom_errorbarh(
      aes(
        xmin = low_plot,
        xmax = high_plot,
        colour = col
      ),
      height = 0,
      linewidth = 0.7
    ) +
    geom_point(
      aes(colour = col),
      shape = 15,
      size = 3.0
    ) +
    scale_colour_identity() +
    scale_x_log10() +
    labs(
      x = "Odds ratio (95% CI)",
      y = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family = BASE_FAMILY
    ) +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(
    file.path(
      DIR_FIG,
      "Supplementary_UVMR_EXPLORATORY_RELAXED.pdf"
    ),
    p,
    width = 7.0,
    height = max(3.0, 0.55 * nrow(d) + 1.6),
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    file.path(
      DIR_FIG,
      "Supplementary_UVMR_EXPLORATORY_RELAXED.tiff"
    ),
    p,
    width = 7.0,
    height = max(3.0, 0.55 * nrow(d) + 1.6),
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  writeLines(
    c(
      "Exploratory relaxed-IV MR legend",
      "",
      "These analyses were triggered only when the genome-wide-significant primary instrument set contained fewer than 3 independent IVs.",
      "Relaxation sequence was pre-specified: P<5×10^-8 → P<1×10^-6 → P<5×10^-6.",
      "F>=10 and LD clumping r2<0.001 within 10,000 kb were never relaxed.",
      "These results are exploratory and must not replace the primary genome-wide-significant MR results."
    ),
    file.path(
      DIR_FIG,
      "Supplementary_UVMR_EXPLORATORY_RELAXED_legend.txt"
    )
  )
}


# ==============================================================================
# 22. MVMR HELPERS
# ==============================================================================

# Use strict BMI/OSA instruments whenever each has >=3 strict IVs.
# Otherwise use the pre-specified selected relaxed tier and mark entire MVMR
# analysis as exploratory.

get_mvmr_source <- function(trait) {

  z <- instrument_sets[[trait]]

  if (nrow(z$strict) >= MIN_IV_TARGET) {
    list(
      dt = z$strict,
      tier = "primary_GWS",
      p = P_TIERS[["primary_GWS"]],
      class = "PRIMARY"
    )
  } else {
    list(
      dt = z$selected,
      tier = z$selected_tier,
      p = z$selected_p,
      class = "EXPLORATORY_RELAXED"
    )
  }
}

BMI_SRC <- get_mvmr_source("BMI")
OSA_SRC <- get_mvmr_source("OSA")

MVMR_CLASS <- if (
  BMI_SRC$class == "PRIMARY" &&
  OSA_SRC$class == "PRIMARY"
) {
  "PRIMARY"
} else {
  "EXPLORATORY_RELAXED"
}

log_msg(
  "MVMR instrument class = ", MVMR_CLASS,
  "; BMI tier=", BMI_SRC$tier,
  "; OSA tier=", OSA_SRC$tier
)


# ==============================================================================
# 23. BUILD UNION MVMR INSTRUMENT SET + RE-CLUMP
# ==============================================================================

mvmr_union_raw <- rbindlist(
  list(
    BMI_SRC$dt[, .(
      SNP,
      P,
      source_exposure = "BMI"
    )],
    OSA_SRC$dt[, .(
      SNP,
      P,
      source_exposure = "OSA"
    )]
  ),
  fill = TRUE
)

mvmr_union_for_clump <- mvmr_union_raw[
  ,
  .(
    P = min(P, na.rm = TRUE)
  ),
  by = SNP
]

# Re-clump union so instruments remain independent across exposures.
mvmr_union_clumped <- local_clump(
  mvmr_union_for_clump,
  trait = "BMI_OSA_MVMR_union",
  clump_p = 1
)

MVMR_SNPS <- mvmr_union_clumped$SNP

fwrite(
  mvmr_union_clumped,
  file.path(
    DIR_INST,
    "MVMR_BMI_OSA_union_reclumped.csv"
  )
)


# ==============================================================================
# 24. PREPARE MVMR MATRICES BY HARMONISING TO BMI
# ==============================================================================

prepare_mvmr <- function(outcome_trait) {

  tag <- paste0(
    "BMI_OSA__to__", outcome_trait
  )

  bmi <- assoc_tables[["BMI"]][SNP %in% MVMR_SNPS]
  osa <- assoc_tables[["OSA"]][SNP %in% MVMR_SNPS]
  out <- assoc_tables[[outcome_trait]][SNP %in% MVMR_SNPS]

  if (
    nrow(bmi) < 3 ||
    nrow(osa) < 3 ||
    nrow(out) < 3
  ) {
    return(
      list(
        status = "NOT_ESTIMABLE_TOO_FEW_OVERLAPPING_SNPS"
      )
    )
  }

  # Align OSA to BMI allele coding
  h_bmi_osa <- TwoSampleMR::harmonise_data(
    to_tsmr_exposure(bmi, "BMI"),
    to_tsmr_outcome(osa, "OSA"),
    action = HARMONISE_ACTION
  )

  h_bmi_osa <- as.data.table(h_bmi_osa)[mr_keep %in% TRUE]

  if (nrow(h_bmi_osa) < 3) {
    return(
      list(
        status = "NOT_ESTIMABLE_AFTER_BMI_OSA_HARMONISATION"
      )
    )
  }

  # Align outcome to the same BMI reference allele coding
  out2 <- out[SNP %in% h_bmi_osa$SNP]
  bmi2 <- bmi[SNP %in% h_bmi_osa$SNP]

  h_bmi_out <- TwoSampleMR::harmonise_data(
    to_tsmr_exposure(bmi2, "BMI"),
    to_tsmr_outcome(out2, outcome_trait),
    action = HARMONISE_ACTION
  )

  h_bmi_out <- as.data.table(h_bmi_out)[mr_keep %in% TRUE]

  # Merge aligned OSA and aligned outcome by SNP.
  d <- merge(
    h_bmi_osa[
      ,
      .(
        SNP,
        beta_BMI = beta.exposure,
        se_BMI = se.exposure,
        p_BMI = pval.exposure,
        beta_OSA = beta.outcome,
        se_OSA = se.outcome,
        p_OSA = pval.outcome
      )
    ],
    h_bmi_out[
      ,
      .(
        SNP,
        beta_OUT = beta.outcome,
        se_OUT = se.outcome,
        p_OUT = pval.outcome
      )
    ],
    by = "SNP"
  )

  d <- d[
    complete.cases(d)
  ]

  fwrite(
    d,
    file.path(
      DIR_HARM,
      paste0(tag, "__MVMR_harmonised.csv")
    )
  )

  if (nrow(d) < 3) {
    return(
      list(
        status = "NOT_ESTIMABLE_TOO_FEW_FINAL_SNPS",
        data = d
      )
    )
  }

  list(
    status = "OK",
    data = d
  )
}


# ==============================================================================
# 25. RUN MVMR
# ==============================================================================

mvmr_final_results <- list()
mvmr_qc_results <- list()

run_mvmr_one <- function(outcome_trait) {

  tag <- paste0(
    "BMI_OSA__to__", outcome_trait
  )

  prep <- prepare_mvmr(outcome_trait)

  if (prep$status != "OK") {

    log_msg(
      "MVMR ", tag, ": ", prep$status
    )

    return(
      list(
        result = data.table(
          exposure = c("BMI", "OSA"),
          outcome = outcome_trait,
          analysis_class = MVMR_CLASS,
          nsnp = if (!is.null(prep$data)) nrow(prep$data) else 0L,
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = prep$status
        ),
        qc = data.table(
          outcome = outcome_trait,
          status = prep$status
        )
      )
    )
  }

  d <- prep$data

  if (!MVMR_AVAILABLE) {

    return(
      list(
        result = data.table(
          exposure = c("BMI", "OSA"),
          outcome = outcome_trait,
          analysis_class = MVMR_CLASS,
          nsnp = nrow(d),
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = "MVMR_PACKAGE_UNAVAILABLE"
        ),
        qc = data.table(
          outcome = outcome_trait,
          status = "MVMR_PACKAGE_UNAVAILABLE"
        )
      )
    )
  }

  BX <- as.matrix(
    d[, .(beta_BMI, beta_OSA)]
  )

  colnames(BX) <- c("BMI", "OSA")

  seBX <- as.matrix(
    d[, .(se_BMI, se_OSA)]
  )

  colnames(seBX) <- c("BMI", "OSA")

  rin <- MVMR::format_mvmr(
    BXGs = BX,
    BYG = d$beta_OUT,
    seBXGs = seBX,
    seBYG = d$se_OUT,
    RSID = d$SNP
  )

  ivw <- try(
    MVMR::ivw_mvmr(
      r_input = rin,
      gencov = 0
    ),
    silent = TRUE
  )

  strength <- try(
    MVMR::strength_mvmr(
      r_input = rin,
      gencov = 0
    ),
    silent = TRUE
  )

  pleio <- try(
    MVMR::pleiotropy_mvmr(
      r_input = rin,
      gencov = 0
    ),
    silent = TRUE
  )

  capture.output(
    ivw,
    file = file.path(
      DIR_MVMR,
      paste0(tag, "__ivw_mvmr.txt")
    )
  )

  capture.output(
    strength,
    file = file.path(
      DIR_MVMR,
      paste0(tag, "__conditional_F.txt")
    )
  )

  capture.output(
    pleio,
    file = file.path(
      DIR_MVMR,
      paste0(tag, "__pleiotropy_Q.txt")
    )
  )

  if (inherits(ivw, "try-error")) {

    return(
      list(
        result = data.table(
          exposure = c("BMI", "OSA"),
          outcome = outcome_trait,
          analysis_class = MVMR_CLASS,
          nsnp = nrow(d),
          b = NA_real_,
          se = NA_real_,
          pval = NA_real_,
          status = paste0(
            "MVMR_IVW_ERROR: ",
            as.character(ivw)
          )
        ),
        qc = data.table(
          outcome = outcome_trait,
          status = "MVMR_IVW_ERROR"
        )
      )
    )
  }

  ivw_df <- as.data.frame(ivw)

  # Robustly locate columns
  est_col <- grep(
    "^Estimate$",
    names(ivw_df),
    value = TRUE
  )[1]

  se_col <- grep(
    "Std\\.? ?Error",
    names(ivw_df),
    value = TRUE,
    ignore.case = TRUE
  )[1]

  p_col <- grep(
    "Pr\\(|p",
    names(ivw_df),
    value = TRUE,
    ignore.case = TRUE
  )

  # Prefer explicit Pr(>|t|) if present
  p_col2 <- names(ivw_df)[
    grepl(
      "Pr",
      names(ivw_df),
      fixed = TRUE
    )
  ]

  if (length(p_col2) > 0) {
    p_use <- p_col2[1]
  } else {
    # fall back to last numeric column
    numeric_cols <- names(ivw_df)[
      vapply(ivw_df, is.numeric, logical(1))
    ]
    p_use <- tail(numeric_cols, 1)
  }

  if (
    is.na(est_col) ||
    is.na(se_col) ||
    length(p_use) == 0
  ) {
    stop(
      "Could not parse ivw_mvmr result columns for ",
      outcome_trait,
      ". See saved txt file."
    )
  }

  res <- data.table(
    exposure = c("BMI", "OSA")[seq_len(nrow(ivw_df))],
    outcome = outcome_trait,
    analysis_class = MVMR_CLASS,
    nsnp = nrow(d),
    b = as.numeric(ivw_df[[est_col]]),
    se = as.numeric(ivw_df[[se_col]]),
    pval = as.numeric(ivw_df[[p_use]]),
    status = "OK"
  )

  # Extract conditional F values conservatively.
  condF <- c(NA_real_, NA_real_)

  if (!inherits(strength, "try-error")) {
    s_df <- as.data.frame(strength)
    vals <- unlist(
      lapply(
        s_df,
        function(x) {
          if (is.numeric(x)) x else numeric(0)
        }
      ),
      use.names = FALSE
    )
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 2) {
      condF <- vals[1:2]
    }
  }

  res[, conditional_F := condF]

  qc <- data.table(
    outcome = outcome_trait,
    analysis_class = MVMR_CLASS,
    n_final_snps = nrow(d),
    BMI_conditional_F = condF[1],
    OSA_conditional_F = condF[2],
    covariance_assumption = "gencov=0",
    status = "OK"
  )

  list(
    result = res,
    qc = qc
  )
}

for (oy in c("HFpEF", "AF")) {

  z <- run_mvmr_one(oy)

  mvmr_final_results[[length(mvmr_final_results) + 1]] <- z$result
  mvmr_qc_results[[length(mvmr_qc_results) + 1]] <- z$qc
}

mvmr_res <- rbindlist(
  mvmr_final_results,
  fill = TRUE
)

mvmr_qc <- rbindlist(
  mvmr_qc_results,
  fill = TRUE
)


# ==============================================================================
# 26. MVMR MULTIPLE TESTING
# ==============================================================================

mvmr_res[
  ,
  bonferroni_threshold := 0.05 / MVMR_PLANNED_TESTS
]

idx <- which(is.finite(mvmr_res$pval))

mvmr_res[, BH_FDR_q := NA_real_]

if (length(idx) > 0) {
  mvmr_res$BH_FDR_q[idx] <- p.adjust(
    mvmr_res$pval[idx],
    method = "BH",
    n = MVMR_PLANNED_TESTS
  )
}

mvmr_res[
  ,
  significance_class := fifelse(
    is.finite(pval) &
      pval < bonferroni_threshold,
    "Bonferroni significant",
    fifelse(
      is.finite(pval) &
        pval < 0.05,
      "Nominal",
      "Not significant"
    )
  )
]

mvmr_res[
  ,
  `:=`(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se
  )
]

fwrite(
  mvmr_res,
  file.path(
    DIR_MVMR,
    "STEP3B_MVMR_BMI_OSA_final.csv"
  )
)

fwrite(
  mvmr_qc,
  file.path(
    DIR_MVMR,
    "STEP3B_MVMR_QC.csv"
  )
)


# ==============================================================================
# 27. MVMR FOREST FIGURE
# ==============================================================================

if (nrow(mvmr_res[status == "OK"]) > 0) {

  d <- copy(mvmr_res)

  d[, pair := paste0(
    exposure,
    " → ",
    outcome
  )]

  d[, pair := factor(
    pair,
    levels = rev(c(
      "BMI → HFpEF",
      "OSA → HFpEF",
      "BMI → AF",
      "OSA → AF"
    ))
  )]

  d[, `:=`(
    est_plot = exp(b),
    low_plot = exp(ci_low),
    high_plot = exp(ci_high)
  )]

  d[
    ,
    col := fifelse(
      significance_class == "Bonferroni significant",
      COL_BONF,
      fifelse(
        significance_class == "Nominal",
        COL_NOM,
        COL_NS
      )
    )
  ]

  p <- ggplot(
    d,
    aes(
      y = pair,
      x = est_plot
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "#6E6E6E"
    ) +
    geom_errorbarh(
      aes(
        xmin = low_plot,
        xmax = high_plot,
        colour = col
      ),
      height = 0,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    geom_point(
      aes(colour = col),
      shape = 15,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_colour_identity() +
    scale_x_log10() +
    labs(
      x = "Direct effect OR (95% CI)",
      y = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family = BASE_FAMILY
    ) +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  stub <- if (MVMR_CLASS == "PRIMARY") {
    "Figure5C_MVMR_PRIMARY_BMI_OSA"
  } else {
    "Supplementary_MVMR_EXPLORATORY_BMI_OSA"
  }

  ggsave(
    file.path(
      DIR_FIG,
      paste0(stub, ".pdf")
    ),
    p,
    width = 6.6,
    height = 4.0,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(stub, ".tiff")
    ),
    p,
    width = 6.6,
    height = 4.0,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  fwrite(
    d,
    file.path(
      DIR_FIG,
      paste0(stub, "_source_data.csv")
    )
  )

  writeLines(
    c(
      "MVMR figure legend",
      "",
      "Squares indicate IVW multivariable MR direct-effect estimates and horizontal lines indicate 95% confidence intervals.",
      "BMI and OSA are included simultaneously as exposures.",
      paste0("Analysis class: ", MVMR_CLASS, "."),
      "Conditional F-statistics were estimated using the MVMR package.",
      "The exposure-effect covariance was fixed at zero (gencov=0); this assumes no material sample overlap between the BMI and OSA exposure GWAS.",
      paste0(
        "Bonferroni threshold = ",
        signif(0.05 / MVMR_PLANNED_TESTS, 4),
        " across four pre-specified direct-effect tests."
      )
    ),
    file.path(
      DIR_FIG,
      paste0(stub, "_legend.txt")
    )
  )
}


# ==============================================================================
# 28. FINAL READINESS / INTERPRETATION TABLE
# ==============================================================================

final_status <- merge(
  uvmr_pairs,
  primary_rows[
    ,
    .(
      exposure = exposure_trait,
      outcome = outcome_trait,
      primary_nsnp = nsnp,
      primary_method = method,
      primary_b = b,
      primary_se = se,
      primary_p = pval,
      primary_q = BH_FDR_q,
      primary_status = status
    )
  ],
  by = c("exposure", "outcome"),
  all.x = TRUE
)

# Add selected IV tier used for potential exploratory rescue.
final_status[
  ,
  selected_IV_tier := vapply(
    exposure,
    function(x) instrument_sets[[x]]$selected_tier,
    character(1)
  )
]

final_status[
  ,
  selected_IV_p_threshold := vapply(
    exposure,
    function(x) instrument_sets[[x]]$selected_p,
    numeric(1)
  )
]

final_status[
  ,
  selected_IV_n := vapply(
    exposure,
    function(x) nrow(instrument_sets[[x]]$selected),
    integer(1)
  )
]

fwrite(
  final_status,
  file.path(
    DIR_REP,
    "STEP3B_FINAL_analysis_status.csv"
  )
)


# ==============================================================================
# 29. METHODS / REPORTING NOTES
# ==============================================================================

methods_txt <- c(
  "STEP 3B METHODS NOTES",
  "",
  "Primary UVMR",
  "- Genome-wide significant IVs: P < 5×10^-8.",
  "- Weak-instrument filter: F statistic >=10.",
  "- LD clumping: r2 <0.001 within 10,000 kb using 1000 Genomes Phase 3 EUR reference.",
  "- Exposure and outcome alleles harmonised with TwoSampleMR action=2.",
  "- IVW multiplicative random effects is the primary estimator for >=2 IVs; Wald ratio is used for a single IV.",
  "- Weighted median and MR-Egger are sensitivity estimators when >=3 IVs.",
  "- Cochran Q, Egger intercept, leave-one-out and MR-PRESSO are saved where estimable.",
  "- Bonferroni correction uses the six pre-specified primary UVMR tests (0.05/6), with BH-FDR also reported.",
  "",
  "Automatic sparse-IV handling",
  "- If the P<5×10^-8 set contains fewer than 3 independent IVs, the script automatically evaluates P<1×10^-6, then P<5×10^-6.",
  "- The first relaxed tier reaching >=3 independent IVs is used for a separate exploratory analysis.",
  "- F>=10 and LD clumping thresholds are NOT relaxed.",
  "- Exploratory relaxed results are never substituted into the primary Figure 5 results.",
  "- If even P<5×10^-6 is insufficient, the analysis remains sparse/not estimable.",
  "",
  "MVMR",
  "- Exposures: BMI and OSA; outcomes: HFpEF and AF.",
  "- A union of exposure instruments is re-clumped at r2<0.001, 10,000 kb.",
  "- Direct effects are estimated with IVW MVMR.",
  "- Conditional F-statistics and the MVMR Q statistic are evaluated.",
  "- gencov is fixed at 0; interpretation assumes no material overlap between the BMI and OSA exposure GWAS.",
  "",
  "Steiger note",
  "- Automated Steiger filtering is intentionally not used as a primary filter here because AF, HFpEF and OSA are binary traits and reliable Steiger conversion requires appropriate binary-trait prevalence / liability-scale information.",
  "- Bidirectional MR itself provides the pre-specified directional comparison."
)

writeLines(
  methods_txt,
  file.path(
    DIR_REP,
    "STEP3B_METHODS_NOTES.txt"
  )
)


# ==============================================================================
# 30. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_REP,
    "STEP3B_sessionInfo.txt"
  )
)


# ==============================================================================
# 31. CLOSE DATABASE
# ==============================================================================

dbDisconnect(
  con,
  shutdown = TRUE
)


# ==============================================================================
# 32. FINAL CONSOLE
# ==============================================================================

log_msg("============================================================")
log_msg("STEP 3B V2 COMPLETE")
log_msg("============================================================")
log_msg("Output root: ", OUT_ROOT)
log_msg("")
log_msg("Please send back these files first:")
log_msg("1) 07_reports/STEP3B_instrument_selection_QC.csv")
log_msg("2) 03_uvmr/STEP3B_UVMR_PRIMARY_final.csv")
log_msg("3) 03_uvmr/STEP3B_UVMR_EXPLORATORY_RELAXED_final.csv (if generated)")
log_msg("4) 04_mvmr/STEP3B_MVMR_BMI_OSA_final.csv")
log_msg("5) 04_mvmr/STEP3B_MVMR_QC.csv")
log_msg("6) 07_reports/STEP3B_FINAL_analysis_status.csv")
log_msg("7) 07_reports/STEP3B_master_console.log")
log_msg("")
log_msg("Figures are in: 06_figures")
