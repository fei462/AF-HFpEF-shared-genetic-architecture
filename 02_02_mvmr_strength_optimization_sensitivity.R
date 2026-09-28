# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP3C_MVMR_STRENGTH_OPTIMIZATION(1).R
# See repository README.md for execution order and external dependencies.

# ==============================================================================
# STEP 3C — Strength-Optimized MVMR for BMI + OSA -> HFpEF / AF
# Project: AF–HFpEF–BMI–OSA shared genetic architecture
#
# WHY STEP 3C?
#   STEP3B V2 strict MVMR used genome-wide significant (P<5e-8) BMI/OSA IVs,
#   but conditional F-statistics were <10 (BMI ~7.24; OSA ~4.02), indicating
#   conditional weak-instrument bias.
#
# SCIENTIFIC PRINCIPLE
#   We DO NOT tune instrument thresholds against the outcome P-value.
#   We evaluate a pre-specified exposure-only threshold grid:
#
#      BMI tier: 5e-8, 1e-6, 5e-6
#      OSA tier: 5e-8, 1e-6, 5e-6
#
#   This creates 9 combinations. For each combination:
#     1) take previously clumped BMI and OSA IV sets from STEP3B V2,
#     2) take their union,
#     3) re-clump the union genome-wide at r2<0.001, 10,000 kb,
#     4) harmonize BMI + OSA + outcome,
#     5) calculate conditional F-statistics,
#     6) fit IVW MVMR and Q-statistic.
#
#   Selection of the "optimized" tier is based ONLY on instrument-strength QC:
#     - both conditional F >=10;
#     - choose the LEAST relaxed threshold combination;
#     - if tied, choose the combination with the largest minimum conditional F.
#
#   Outcome effect size/P-value is NEVER used to select the tier.
#
# IMPORTANT INTERPRETATION
#   If the selected tier contains P>5e-8 IVs, the optimized MVMR is an
#   EXPLORATORY / SENSITIVITY MVMR. It does NOT replace the strict GWS MVMR.
#
# OPTIONAL WEAK-INSTRUMENT ROBUST Q-MINIMIZATION
#   MVMR::qhet_mvmr() requires a phenotypic correlation matrix for BMI and OSA.
#   Because a genetic correlation is NOT a phenotypic correlation, this script
#   will NOT invent one. If you later have a defensible BMI–OSA phenotypic
#   correlation from individual-level data or a suitable external study, enter
#   it in BMI_OSA_PHENO_COR below and the script will run qhet_mvmr().
#
# INPUT ROOT
#   D:/A/data/STEP3_MR_V2
#
# OUTPUT ROOT
#   D:/A/data/STEP3_MR_STEP3C
#
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

DATA_DIR <- "D:/A/data"
STEP3B_DIR <- file.path(DATA_DIR, "STEP3_MR_V2")
OUT_DIR <- file.path(DATA_DIR, "STEP3_MR_STEP3C")

CLUMP_R2 <- 0.001
CLUMP_KB <- 10000
F_PASS <- 10
HARMONISE_ACTION <- 2L

# Optional robust qhet MVMR:
# Put a defensible phenotypic correlation here, e.g. 0.30.
# Leave NA to skip qhet_mvmr (recommended unless you have a justified value).
BMI_OSA_PHENO_COR <- NA_real_

# If qhet is enabled:
QHET_BOOTSTRAP_CI <- TRUE
QHET_ITERATIONS <- 1000L

FIG_DPI <- 600
BASE_FAMILY <- "Arial"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DIR_INST <- file.path(OUT_DIR, "01_instrument_grid")
DIR_HARM <- file.path(OUT_DIR, "02_harmonized")
DIR_MVMR <- file.path(OUT_DIR, "03_mvmr_results")
DIR_FIG  <- file.path(OUT_DIR, "04_figures")
DIR_REP  <- file.path(OUT_DIR, "05_reports")
DIR_CACHE <- file.path(OUT_DIR, "06_cache")

for (d in c(DIR_INST, DIR_HARM, DIR_MVMR, DIR_FIG, DIR_REP, DIR_CACHE)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

options(repos = c(
  MRCIEU = "https://mrcieu.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (p in c(
  "data.table", "DBI", "duckdb", "ggplot2",
  "TwoSampleMR", "MVMR", "ieugwasr", "remotes"
)) {
  try(install_if_missing(p), silent = TRUE)
}

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("TwoSampleMR is required.")
}

if (!requireNamespace("MVMR", quietly = TRUE)) {
  stop("MVMR is required.")
}

# STEP3C requires the corrected recent MVMR implementation.
# v0.4.7 corrected qhet_mvmr(); v0.4.8 corrected strhet_mvmr().
if (utils::packageVersion("MVMR") < "0.4.8") {
  message("Updating MVMR to >=0.4.8 ...")
  try(
    install.packages(
      "MVMR",
      repos = c(
        "https://mrcieu.r-universe.dev",
        "https://cloud.r-project.org"
      )
    ),
    silent = TRUE
  )
}

if (utils::packageVersion("MVMR") < "0.4.8") {
  stop(
    "STEP3C requires MVMR >=0.4.8. Current version: ",
    as.character(utils::packageVersion("MVMR")),
    "\nPlease restart R and rerun after updating MVMR."
  )
}

if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  stop("ieugwasr is required.")
}

library(data.table)
library(DBI)
library(duckdb)
library(ggplot2)


# ==============================================================================
# 2. LOGGING
# ==============================================================================

LOG_FILE <- file.path(DIR_REP, "STEP3C_master_console.log")

if (file.exists(LOG_FILE)) unlink(LOG_FILE)

log_msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x, "\n"),
    file = LOG_FILE,
    append = TRUE
  )
}

log_msg("============================================================")
log_msg("STEP 3C START")
log_msg("Strength-optimized MVMR")
log_msg("============================================================")


# ==============================================================================
# 3. CHECK STEP3B V2 INPUTS
# ==============================================================================

if (!dir.exists(STEP3B_DIR)) {
  stop(
    "STEP3B V2 directory not found: ",
    STEP3B_DIR
  )
}

REF_QC_FILE <- file.path(
  STEP3B_DIR,
  "07_reports",
  "STEP3B_V2_LD_reference_QC.csv"
)

if (!file.exists(REF_QC_FILE)) {
  stop(
    "Missing V2 LD reference QC file:\n",
    REF_QC_FILE
  )
}

ref_qc <- fread(REF_QC_FILE)

if (
  !"chromosome" %in% names(ref_qc) ||
  !"bfile_prefix" %in% names(ref_qc)
) {
  stop("Unexpected LD reference QC format.")
}

ref_qc <- ref_qc[chromosome %in% 1:22]
setorder(ref_qc, chromosome)

if (nrow(ref_qc) != 22) {
  stop(
    "Expected 22 chromosome-specific EUR reference entries, found ",
    nrow(ref_qc)
  )
}

BFILE_CHR <- ref_qc$bfile_prefix

has_triplet <- function(prefix) {
  all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
}

if (!all(vapply(BFILE_CHR, has_triplet, logical(1)))) {
  bad <- which(!vapply(BFILE_CHR, has_triplet, logical(1)))
  stop(
    "Missing reference BED/BIM/FAM for chromosome(s): ",
    paste(bad, collapse = ", ")
  )
}

log_msg("22-chromosome EUR reference: PASS")


# ==============================================================================
# 4. FIND PLINK
# ==============================================================================

find_plink <- function() {

  hits <- list.files(
    DATA_DIR,
    pattern = "^plink(\\.exe)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

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

PLINK_BIN <- find_plink()

if (is.na(PLINK_BIN) || !file.exists(PLINK_BIN)) {
  stop("PLINK binary could not be located.")
}

log_msg("PLINK = ", PLINK_BIN)


# ==============================================================================
# 5. LOAD PREVIOUSLY CLUMPED BMI / OSA INSTRUMENT TIERS
# ==============================================================================

INST_DIR <- file.path(
  STEP3B_DIR,
  "01_instruments"
)

tier_meta <- data.table(
  tier = c(
    "GWS",
    "P1e_6",
    "P5e_6"
  ),
  file_tag = c(
    "primary_GWS",
    "exploratory_P1e_6",
    "exploratory_P5e_6"
  ),
  p_threshold = c(
    5e-8,
    1e-6,
    5e-6
  ),
  relax_rank = c(
    0L,
    1L,
    2L
  )
)

load_tier <- function(trait, file_tag) {

  f <- file.path(
    INST_DIR,
    paste0(
      trait,
      "__",
      file_tag,
      "__clumped_instruments.csv"
    )
  )

  if (!file.exists(f)) {
    stop("Missing STEP3B V2 instrument file:\n", f)
  }

  d <- fread(f)

  req <- c("SNP", "P")
  miss <- setdiff(req, names(d))

  if (length(miss) > 0) {
    stop(
      "Instrument file missing columns ",
      paste(miss, collapse = ", "),
      ":\n", f
    )
  }

  d[, source_trait := trait]
  d
}

BMI_TIERS <- setNames(
  lapply(
    tier_meta$file_tag,
    function(x) load_tier("BMI", x)
  ),
  tier_meta$tier
)

OSA_TIERS <- setNames(
  lapply(
    tier_meta$file_tag,
    function(x) load_tier("OSA", x)
  ),
  tier_meta$tier
)

tier_counts <- rbindlist(
  lapply(
    tier_meta$tier,
    function(tt) {
      data.table(
        tier = tt,
        BMI_n = nrow(BMI_TIERS[[tt]]),
        OSA_n = nrow(OSA_TIERS[[tt]])
      )
    }
  )
)

fwrite(
  tier_counts,
  file.path(
    DIR_REP,
    "STEP3C_input_tier_counts.csv"
  )
)


# ==============================================================================
# 6. BUILD 1000G rsID -> CHROMOSOME INDEX
# ==============================================================================

log_msg("Building EUR reference SNP index...")

ref_index <- rbindlist(
  lapply(
    1:22,
    function(cc) {

      bim <- fread(
        paste0(BFILE_CHR[cc], ".bim"),
        header = FALSE,
        select = c(1, 2, 4),
        col.names = c("CHR", "SNP", "BP"),
        showProgress = FALSE
      )

      bim[
        ,
        .(
          SNP,
          REF_CHR = as.integer(CHR),
          REF_BP = as.integer(BP)
        )
      ]
    }
  )
)

ref_index <- unique(
  ref_index,
  by = "SNP"
)

log_msg(
  "Reference index SNPs = ",
  nrow(ref_index)
)


# ==============================================================================
# 7. CHROMOSOME-SAFE RE-CLUMPING OF BMI + OSA UNION
# ==============================================================================

reclump_union <- function(
  bmi_inst,
  osa_inst,
  combo_id
) {

  u <- rbindlist(
    list(
      bmi_inst[, .(
        SNP,
        P_union = P,
        source = "BMI"
      )],
      osa_inst[, .(
        SNP,
        P_union = P,
        source = "OSA"
      )]
    ),
    fill = TRUE
  )

  u <- u[
    ,
    .(
      P_union = min(P_union, na.rm = TRUE),
      source = paste(
        sort(unique(source)),
        collapse = "+"
      )
    ),
    by = SNP
  ]

  umap <- merge(
    u,
    ref_index,
    by = "SNP",
    all = FALSE
  )

  match_rate <- nrow(umap) / nrow(u)

  if (
    nrow(u) >= 10 &&
    match_rate < 0.50
  ) {
    stop(
      "Reference match rate <50% for ",
      combo_id
    )
  }

  kept <- character(0)
  qc_chr <- list()

  for (cc in sort(unique(umap$REF_CHR))) {

    if (!cc %in% 1:22) next

    d <- umap[REF_CHR == cc]

    inp <- unique(
      data.table(
        rsid = d$SNP,
        pval = d$P_union,
        id = paste0(combo_id, "_chr", cc)
      ),
      by = "rsid"
    )

    cl <- try(
      ieugwasr::ld_clump(
        dat = as.data.frame(inp),
        clump_kb = CLUMP_KB,
        clump_r2 = CLUMP_R2,
        clump_p = 1,
        pop = "EUR",
        bfile = BFILE_CHR[cc],
        plink_bin = PLINK_BIN
      ),
      silent = TRUE
    )

    if (inherits(cl, "try-error")) {
      stop(
        "Re-clumping failed: ",
        combo_id,
        " chr", cc,
        "\n",
        as.character(cl)
      )
    }

    snp_col <- if ("rsid" %in% names(cl)) {
      "rsid"
    } else if ("variant" %in% names(cl)) {
      "variant"
    } else {
      stop("Cannot identify clumped SNP column.")
    }

    kk <- unique(as.character(cl[[snp_col]]))

    kept <- c(kept, kk)

    qc_chr[[length(qc_chr) + 1]] <- data.table(
      combo = combo_id,
      chromosome = cc,
      n_input = nrow(d),
      n_kept = length(kk)
    )
  }

  kept <- unique(kept)

  out <- umap[SNP %in% kept]
  setorder(out, P_union)

  fwrite(
    out,
    file.path(
      DIR_INST,
      paste0(
        combo_id,
        "__union_reclumped.csv"
      )
    )
  )

  fwrite(
    rbindlist(qc_chr, fill = TRUE),
    file.path(
      DIR_INST,
      paste0(
        combo_id,
        "__reclump_chr_QC.csv"
      )
    )
  )

  list(
    data = out,
    n_union = nrow(u),
    n_mapped = nrow(umap),
    match_rate = match_rate,
    n_reclumped = nrow(out)
  )
}


# ==============================================================================
# 8. BUILD ALL 9 PRE-SPECIFIED EXPOSURE-TIER COMBINATIONS
# ==============================================================================

combo_grid <- CJ(
  BMI_tier = tier_meta$tier,
  OSA_tier = tier_meta$tier,
  unique = TRUE
)

combo_grid <- merge(
  combo_grid,
  tier_meta[
    ,
    .(
      BMI_tier = tier,
      BMI_p = p_threshold,
      BMI_rank = relax_rank
    )
  ],
  by = "BMI_tier"
)

combo_grid <- merge(
  combo_grid,
  tier_meta[
    ,
    .(
      OSA_tier = tier,
      OSA_p = p_threshold,
      OSA_rank = relax_rank
    )
  ],
  by = "OSA_tier"
)

combo_grid[
  ,
  relax_score := BMI_rank + OSA_rank
]

combo_grid[
  ,
  combo_id := paste0(
    "BMI_", BMI_tier,
    "__OSA_", OSA_tier
  )
]

setorder(
  combo_grid,
  relax_score,
  BMI_rank,
  OSA_rank
)

union_sets <- list()
union_qc <- list()

for (i in seq_len(nrow(combo_grid))) {

  rr <- combo_grid[i]

  log_msg(
    "Re-clumping union: ",
    rr$combo_id
  )

  z <- reclump_union(
    BMI_TIERS[[rr$BMI_tier]],
    OSA_TIERS[[rr$OSA_tier]],
    rr$combo_id
  )

  union_sets[[rr$combo_id]] <- z$data

  union_qc[[length(union_qc) + 1]] <- data.table(
    combo_id = rr$combo_id,
    BMI_tier = rr$BMI_tier,
    OSA_tier = rr$OSA_tier,
    BMI_p = rr$BMI_p,
    OSA_p = rr$OSA_p,
    BMI_rank = rr$BMI_rank,
    OSA_rank = rr$OSA_rank,
    relax_score = rr$relax_score,
    n_union_before_reclump = z$n_union,
    n_mapped_1000G = z$n_mapped,
    reference_match_rate = z$match_rate,
    n_union_after_reclump = z$n_reclumped
  )
}

union_qc <- rbindlist(union_qc)

fwrite(
  union_qc,
  file.path(
    DIR_REP,
    "STEP3C_union_grid_QC.csv"
  )
)


# ==============================================================================
# 9. GWAS FILE DISCOVERY
# ==============================================================================

all_files <- list.files(
  DATA_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

find_exact <- function(fname, prefer_regex = NULL) {

  x <- all_files[
    tolower(basename(all_files)) ==
      tolower(fname)
  ]

  if (
    !is.null(prefer_regex) &&
    length(x) > 0
  ) {

    y <- x[
      grepl(
        prefer_regex,
        x,
        ignore.case = TRUE
      )
    ]

    if (length(y) > 0) {
      return(y[which.min(nchar(y))])
    }
  }

  if (length(x) == 0) return(NA_character_)

  x[which.min(nchar(x))]
}

AF_PATH <- find_exact("GCST90624412.tsv.gz")
HF_PATH <- find_exact(
  "FORMAT-METAL_Pheno4_EUR.tsv.gz",
  prefer_regex = "HERMES|STEP1B_HFpEF_PREP"
)
BMI_PATH <- find_exact("SNP_gwas_mc_merge_nogc.tbl.uniq.gz")
OSA_PATH <- find_exact("finngen_R9_G6_SLEEPAPNO.gz")

paths <- c(
  BMI = BMI_PATH,
  OSA = OSA_PATH,
  AF = AF_PATH,
  HFpEF = HF_PATH
)

if (
  any(is.na(paths)) ||
  any(!file.exists(paths))
) {
  stop(
    "Missing original GWAS file(s):\n",
    paste(
      names(paths),
      paths,
      sep = " = ",
      collapse = "\n"
    )
  )
}

log_msg("GWAS paths resolved.")


# ==============================================================================
# 10. DUCKDB + EXPLICIT GWAS COLUMN SPECS
# ==============================================================================

DB_FILE <- file.path(
  DIR_CACHE,
  "STEP3C.duckdb"
)

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_FILE,
  read_only = FALSE
)

dbExecute(con, "SET threads=4;")

sql_path <- function(x) {
  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(
        x,
        winslash = "/",
        mustWork = TRUE
      )
    ),
    "'"
  )
}

qid <- function(x) {
  paste0(
    '"',
    gsub('"', '""', x, fixed = TRUE),
    '"'
  )
}

schema_cols <- function(path) {
  dbGetQuery(
    con,
    sprintf(
      "DESCRIBE SELECT * FROM read_csv_auto(%s, header=true, sample_size=100000);",
      sql_path(path)
    )
  )$column_name
}

pick <- function(cols, candidates) {

  for (x in candidates) {
    ii <- which(tolower(cols) == tolower(x))
    if (length(ii) > 0) {
      return(cols[ii[1]])
    }
  }

  stop(
    "Column not found. Tried: ",
    paste(candidates, collapse = ", "),
    "\nAvailable: ",
    paste(cols, collapse = " | ")
  )
}

make_spec <- function(trait, path) {

  cc <- schema_cols(path)

  list(
    trait = trait,
    path = path,
    SNP = pick(
      cc,
      c(
        "rs_id", "rsID", "rsids",
        "SNP", "MarkerName"
      )
    ),
    A1 = pick(
      cc,
      c(
        "effect_allele",
        "A1",
        "alt"
      )
    ),
    A2 = pick(
      cc,
      c(
        "other_allele",
        "A2",
        "ref"
      )
    ),
    BETA = pick(
      cc,
      c(
        "beta",
        "A1_beta",
        "b"
      )
    ),
    SE = pick(
      cc,
      c(
        "standard_error",
        "A1_se",
        "sebeta",
        "se"
      )
    ),
    P = pick(
      cc,
      c(
        "p_value",
        "A1_p",
        "pval",
        "p"
      )
    ),
    EAF = pick(
      cc,
      c(
        "effect_allele_frequency",
        "A1_freq",
        "af_alt",
        "Freq1.Hapmap",
        "eaf"
      )
    )
  )
}

specs <- list(
  BMI = make_spec("BMI", BMI_PATH),
  OSA = make_spec("OSA", OSA_PATH),
  AF = make_spec("AF", AF_PATH),
  HFpEF = make_spec("HFpEF", HF_PATH)
)

header_map <- rbindlist(
  lapply(
    specs,
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
        EAF = s$EAF
      )
    }
  )
)

fwrite(
  header_map,
  file.path(
    DIR_REP,
    "STEP3C_header_mapping.csv"
  )
)


# ==============================================================================
# 11. QUERY REQUIRED SNPs FROM ORIGINAL GWAS
# ==============================================================================

all_needed_snps <- unique(
  unlist(
    lapply(
      union_sets,
      function(x) x$SNP
    )
  )
)

wanted <- data.frame(
  SNP = all_needed_snps,
  stringsAsFactors = FALSE
)

dbWriteTable(
  con,
  "wanted_snps",
  wanted,
  overwrite = TRUE,
  temporary = TRUE
)

query_spec <- function(s) {

  snp_expr <- sprintf(
    "regexp_extract(CAST(%s AS VARCHAR), 'rs[0-9]+')",
    qid(s$SNP)
  )

  sql <- sprintf(
"
WITH x AS (
  SELECT
    %s AS SNP,
    upper(CAST(%s AS VARCHAR)) AS A1,
    upper(CAST(%s AS VARCHAR)) AS A2,
    TRY_CAST(%s AS DOUBLE) AS BETA,
    TRY_CAST(%s AS DOUBLE) AS SE,
    TRY_CAST(%s AS DOUBLE) AS P,
    TRY_CAST(%s AS DOUBLE) AS EAF
  FROM read_csv_auto(%s, header=true, sample_size=100000)
),
q AS (
  SELECT x.*
  FROM x
  INNER JOIN wanted_snps w
    ON x.SNP = w.SNP
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
),
d AS (
  SELECT *,
    row_number() OVER (
      PARTITION BY SNP
      ORDER BY P ASC
    ) AS rn
  FROM q
)
SELECT
  SNP, A1, A2, BETA, SE, P, EAF
FROM d
WHERE rn=1
",
    snp_expr,
    qid(s$A1),
    qid(s$A2),
    qid(s$BETA),
    qid(s$SE),
    qid(s$P),
    qid(s$EAF),
    sql_path(s$path)
  )

  as.data.table(
    dbGetQuery(
      con,
      sql
    )
  )
}

assoc <- list()

for (tr in names(specs)) {

  log_msg(
    "Extracting associations: ",
    tr
  )

  assoc[[tr]] <- query_spec(
    specs[[tr]]
  )

  fwrite(
    assoc[[tr]],
    file.path(
      DIR_CACHE,
      paste0(
        tr,
        "__required_SNP_associations.csv.gz"
      )
    )
  )

  log_msg(
    "  found ",
    nrow(assoc[[tr]]),
    " / ",
    length(all_needed_snps)
  )
}


# ==============================================================================
# 12. TwoSampleMR FORMAT HELPERS
# ==============================================================================

fmt_exp <- function(d, trait) {

  TwoSampleMR::format_data(
    data.frame(
      SNP = d$SNP,
      beta = d$BETA,
      se = d$SE,
      effect_allele = d$A1,
      other_allele = d$A2,
      eaf = d$EAF,
      pval = d$P,
      Phenotype = trait,
      id = trait,
      stringsAsFactors = FALSE
    ),
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    id_col = "id"
  )
}

fmt_out <- function(d, trait) {

  TwoSampleMR::format_data(
    data.frame(
      SNP = d$SNP,
      beta = d$BETA,
      se = d$SE,
      effect_allele = d$A1,
      other_allele = d$A2,
      eaf = d$EAF,
      pval = d$P,
      Phenotype = trait,
      id = trait,
      stringsAsFactors = FALSE
    ),
    type = "outcome",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    id_col = "id"
  )
}


# ==============================================================================
# 13. PREPARE HARMONIZED MVMR MATRIX
# ==============================================================================

prepare_mvmr <- function(
  snps,
  outcome_trait,
  combo_id
) {

  bmi <- assoc$BMI[SNP %in% snps]
  osa <- assoc$OSA[SNP %in% snps]
  out <- assoc[[outcome_trait]][SNP %in% snps]

  common0 <- Reduce(
    intersect,
    list(
      bmi$SNP,
      osa$SNP,
      out$SNP
    )
  )

  if (length(common0) < 3) {
    return(
      list(
        status = "TOO_FEW_COMMON_SNPS",
        data = data.table()
      )
    )
  }

  bmi <- bmi[SNP %in% common0]
  osa <- osa[SNP %in% common0]
  out <- out[SNP %in% common0]

  # Align OSA to BMI coding.
  h12 <- as.data.table(
    TwoSampleMR::harmonise_data(
      fmt_exp(bmi, "BMI"),
      fmt_out(osa, "OSA"),
      action = HARMONISE_ACTION
    )
  )

  h12 <- h12[mr_keep %in% TRUE]

  if (nrow(h12) < 3) {
    return(
      list(
        status = "TOO_FEW_AFTER_BMI_OSA_HARMONIZATION",
        data = data.table()
      )
    )
  }

  # Align outcome to the same BMI coding.
  bmi2 <- bmi[SNP %in% h12$SNP]
  out2 <- out[SNP %in% h12$SNP]

  h1y <- as.data.table(
    TwoSampleMR::harmonise_data(
      fmt_exp(bmi2, "BMI"),
      fmt_out(out2, outcome_trait),
      action = HARMONISE_ACTION
    )
  )

  h1y <- h1y[mr_keep %in% TRUE]

  d <- merge(
    h12[
      ,
      .(
        SNP,
        beta_BMI = beta.exposure,
        se_BMI = se.exposure,
        beta_OSA = beta.outcome,
        se_OSA = se.outcome
      )
    ],
    h1y[
      ,
      .(
        SNP,
        beta_OUT = beta.outcome,
        se_OUT = se.outcome
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
      paste0(
        combo_id,
        "__to__",
        outcome_trait,
        "__harmonized.csv"
      )
    )
  )

  if (nrow(d) < 3) {
    return(
      list(
        status = "TOO_FEW_FINAL_SNPS",
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
# 14. ROBUST RESULT PARSERS
# ==============================================================================

extract_two_F <- function(x) {

  if (inherits(x, "try-error")) {
    return(c(NA_real_, NA_real_))
  }

  z <- as.data.frame(x)

  vals <- unlist(
    lapply(
      z,
      function(v) {
        if (is.numeric(v)) v else numeric(0)
      }
    ),
    use.names = FALSE
  )

  vals <- vals[is.finite(vals)]

  if (length(vals) < 2) {
    return(c(NA_real_, NA_real_))
  }

  vals[1:2]
}

parse_ivw <- function(x) {

  if (inherits(x, "try-error")) {
    return(
      data.table(
        exposure = c("BMI", "OSA"),
        b = NA_real_,
        se = NA_real_,
        pval = NA_real_
      )
    )
  }

  z <- as.data.frame(x)

  est_col <- names(z)[
    grepl(
      "^Estimate$",
      names(z),
      ignore.case = TRUE
    )
  ][1]

  se_col <- names(z)[
    grepl(
      "Std\\.? ?Error",
      names(z),
      ignore.case = TRUE
    )
  ][1]

  p_cols <- names(z)[
    grepl(
      "Pr",
      names(z),
      fixed = TRUE
    )
  ]

  if (length(p_cols) == 0) {
    nums <- names(z)[
      vapply(
        z,
        is.numeric,
        logical(1)
      )
    ]
    p_col <- tail(nums, 1)
  } else {
    p_col <- p_cols[1]
  }

  data.table(
    exposure = c("BMI", "OSA")[seq_len(nrow(z))],
    b = as.numeric(z[[est_col]]),
    se = as.numeric(z[[se_col]]),
    pval = as.numeric(z[[p_col]])
  )
}

parse_Q <- function(x) {

  if (inherits(x, "try-error")) {
    return(
      list(
        Q = NA_real_,
        p = NA_real_
      )
    )
  }

  vals <- unlist(
    lapply(
      as.data.frame(x),
      function(v) {
        if (is.numeric(v)) v else numeric(0)
      }
    ),
    use.names = FALSE
  )

  vals <- vals[is.finite(vals)]

  if (length(vals) >= 2) {
    return(
      list(
        Q = vals[1],
        p = vals[length(vals)]
      )
    )
  }

  list(
    Q = NA_real_,
    p = NA_real_
  )
}


# ==============================================================================
# 15. RUN ALL 9 COMBINATIONS FOR BOTH OUTCOMES
# ==============================================================================

grid_qc <- list()
grid_ivw <- list()
rinput_store <- list()

for (i in seq_len(nrow(combo_grid))) {

  cg <- combo_grid[i]
  snps <- union_sets[[cg$combo_id]]$SNP

  for (oy in c("HFpEF", "AF")) {

    tag <- paste0(
      cg$combo_id,
      "__to__",
      oy
    )

    log_msg("MVMR grid: ", tag)

    pp <- prepare_mvmr(
      snps = snps,
      outcome_trait = oy,
      combo_id = cg$combo_id
    )

    if (pp$status != "OK") {

      grid_qc[[length(grid_qc) + 1]] <- data.table(
        combo_id = cg$combo_id,
        outcome = oy,
        BMI_tier = cg$BMI_tier,
        OSA_tier = cg$OSA_tier,
        BMI_p = cg$BMI_p,
        OSA_p = cg$OSA_p,
        BMI_rank = cg$BMI_rank,
        OSA_rank = cg$OSA_rank,
        relax_score = cg$relax_score,
        n_final = nrow(pp$data),
        F_BMI = NA_real_,
        F_OSA = NA_real_,
        min_conditional_F = NA_real_,
        both_F_ge10 = FALSE,
        Q = NA_real_,
        Q_p = NA_real_,
        status = pp$status
      )

      next
    }

    d <- pp$data

    BX <- as.matrix(
      d[, .(
        beta_BMI,
        beta_OSA
      )]
    )

    colnames(BX) <- c("BMI", "OSA")

    seBX <- as.matrix(
      d[, .(
        se_BMI,
        se_OSA
      )]
    )

    colnames(seBX) <- c("BMI", "OSA")

    rin <- MVMR::format_mvmr(
      BXGs = BX,
      BYG = d$beta_OUT,
      seBXGs = seBX,
      seBYG = d$se_OUT,
      RSID = d$SNP
    )

    rinput_store[[tag]] <- rin

    st <- try(
      MVMR::strength_mvmr(
        r_input = rin,
        gencov = 0
      ),
      silent = TRUE
    )

    Fv <- extract_two_F(st)

    ivw <- try(
      MVMR::ivw_mvmr(
        r_input = rin,
        gencov = 0
      ),
      silent = TRUE
    )

    pl <- try(
      MVMR::pleiotropy_mvmr(
        r_input = rin,
        gencov = 0
      ),
      silent = TRUE
    )

    # Additional Q-minimization based strength diagnostic
    st_het <- try(
      MVMR::strhet_mvmr(
        r_input = rin,
        gencov = 0
      ),
      silent = TRUE
    )

    Fhet <- extract_two_F(st_het)

    qv <- parse_Q(pl)

    ivw_dt <- parse_ivw(ivw)

    ivw_dt[
      ,
      `:=`(
        combo_id = cg$combo_id,
        outcome = oy,
        BMI_tier = cg$BMI_tier,
        OSA_tier = cg$OSA_tier,
        BMI_p = cg$BMI_p,
        OSA_p = cg$OSA_p,
        relax_score = cg$relax_score,
        n_final = nrow(d),
        F_BMI = Fv[1],
        F_OSA = Fv[2],
        min_conditional_F = min(
          Fv,
          na.rm = TRUE
        ),
        Q = qv$Q,
        Q_p = qv$p
      )
    ]

    grid_ivw[[length(grid_ivw) + 1]] <- ivw_dt

    grid_qc[[length(grid_qc) + 1]] <- data.table(
      combo_id = cg$combo_id,
      outcome = oy,
      BMI_tier = cg$BMI_tier,
      OSA_tier = cg$OSA_tier,
      BMI_p = cg$BMI_p,
      OSA_p = cg$OSA_p,
      BMI_rank = cg$BMI_rank,
      OSA_rank = cg$OSA_rank,
      relax_score = cg$relax_score,
      n_final = nrow(d),
      F_BMI = Fv[1],
      F_OSA = Fv[2],
      min_conditional_F = min(
        Fv,
        na.rm = TRUE
      ),
      both_F_ge10 =
        is.finite(Fv[1]) &&
        is.finite(Fv[2]) &&
        Fv[1] >= F_PASS &&
        Fv[2] >= F_PASS,
      F_BMI_strhet = Fhet[1],
      F_OSA_strhet = Fhet[2],
      Q = qv$Q,
      Q_p = qv$p,
      status = "OK"
    )

    capture.output(
      st,
      file = file.path(
        DIR_MVMR,
        paste0(
          tag,
          "__conditional_F.txt"
        )
      )
    )

    capture.output(
      st_het,
      file = file.path(
        DIR_MVMR,
        paste0(
          tag,
          "__strhet_F.txt"
        )
      )
    )

    capture.output(
      ivw,
      file = file.path(
        DIR_MVMR,
        paste0(
          tag,
          "__IVW.txt"
        )
      )
    )

    capture.output(
      pl,
      file = file.path(
        DIR_MVMR,
        paste0(
          tag,
          "__Q.txt"
        )
      )
    )
  }
}

grid_qc <- rbindlist(
  grid_qc,
  fill = TRUE
)

grid_ivw <- rbindlist(
  grid_ivw,
  fill = TRUE
)

fwrite(
  grid_qc,
  file.path(
    DIR_REP,
    "STEP3C_MVMR_grid_strength_QC.csv"
  )
)

fwrite(
  grid_ivw,
  file.path(
    DIR_MVMR,
    "STEP3C_MVMR_grid_IVW_results.csv"
  )
)


# ==============================================================================
# 16. SELECT LEAST-RELAXED STRENGTH-PASSING TIER — OUTCOME INDEPENDENT OF EFFECT
# ==============================================================================

selected <- list()

for (oy in c("HFpEF", "AF")) {

  d <- grid_qc[
    outcome == oy &
    status == "OK"
  ]

  pass <- d[
    both_F_ge10 == TRUE
  ]

  if (nrow(pass) > 0) {

    min_score <- min(
      pass$relax_score,
      na.rm = TRUE
    )

    cand <- pass[
      relax_score == min_score
    ]

    setorder(
      cand,
      -min_conditional_F,
      BMI_rank,
      OSA_rank
    )

    pickrow <- cand[1]

    pickrow[
      ,
      selection_status := "PASS_BOTH_CONDITIONAL_F_GE10"
    ]

  } else {

    # No valid conventional MVMR tier.
    # Report the strongest tier for diagnostics, but DO NOT call it valid.
    setorder(
      d,
      -min_conditional_F,
      relax_score
    )

    pickrow <- d[1]

    pickrow[
      ,
      selection_status := "NO_TIER_PASSES_BOTH_F_GE10"
    ]
  }

  selected[[oy]] <- pickrow
}

selected_qc <- rbindlist(
  selected,
  fill = TRUE
)

fwrite(
  selected_qc,
  file.path(
    DIR_REP,
    "STEP3C_SELECTED_strength_optimized_tiers.csv"
  )
)


# ==============================================================================
# 17. EXTRACT SELECTED MVMR EFFECTS
# ==============================================================================

selected_effects <- list()

for (oy in c("HFpEF", "AF")) {

  s <- selected_qc[
    outcome == oy
  ]

  eff <- grid_ivw[
    outcome == oy &
    combo_id == s$combo_id
  ]

  if (nrow(eff) == 0) next

  eff[
    ,
    `:=`(
      selection_status = s$selection_status,
      F_BMI_selected = s$F_BMI,
      F_OSA_selected = s$F_OSA,
      selected_min_F = s$min_conditional_F,
      selected_relax_score = s$relax_score
    )
  ]

  eff[
    ,
    analysis_class := fifelse(
      selection_status ==
        "PASS_BOTH_CONDITIONAL_F_GE10" &
        BMI_p <= 5e-8 &
        OSA_p <= 5e-8,
      "PRIMARY_STRONG",
      fifelse(
        selection_status ==
          "PASS_BOTH_CONDITIONAL_F_GE10",
        "EXPLORATORY_STRENGTH_OPTIMIZED",
        "INVALID_WEAK_INSTRUMENT_DIAGNOSTIC"
      )
    )
  ]

  selected_effects[[length(selected_effects) + 1]] <- eff
}

selected_effects <- rbindlist(
  selected_effects,
  fill = TRUE
)

# Multiplicity across four direct effects
selected_effects[
  ,
  bonferroni_threshold := 0.05 / 4
]

selected_effects[
  ,
  BH_FDR_q := p.adjust(
    pval,
    method = "BH",
    n = 4
  )
]

selected_effects[
  ,
  `:=`(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se,
    OR = exp(b),
    OR_low = exp(b - 1.96 * se),
    OR_high = exp(b + 1.96 * se)
  )
]

selected_effects[
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

fwrite(
  selected_effects,
  file.path(
    DIR_MVMR,
    "STEP3C_SELECTED_MVMR_effects.csv"
  )
)


# ==============================================================================
# 18. OPTIONAL Q-HETEROGENEITY MINIMIZATION
# ==============================================================================

qhet_results <- list()

if (
  is.finite(BMI_OSA_PHENO_COR) &&
  abs(BMI_OSA_PHENO_COR) < 1
) {

  log_msg(
    "Running qhet_mvmr with phenotypic rho = ",
    BMI_OSA_PHENO_COR
  )

  pcor <- matrix(
    c(
      1,
      BMI_OSA_PHENO_COR,
      BMI_OSA_PHENO_COR,
      1
    ),
    nrow = 2,
    byrow = TRUE
  )

  for (oy in c("HFpEF", "AF")) {

    ss <- selected_qc[outcome == oy]
    tag <- paste0(
      ss$combo_id,
      "__to__",
      oy
    )

    rin <- rinput_store[[tag]]

    if (is.null(rin)) next

    qh <- try(
      MVMR::qhet_mvmr(
        r_input = rin,
        pcor = pcor,
        CI = QHET_BOOTSTRAP_CI,
        iterations = QHET_ITERATIONS
      ),
      silent = TRUE
    )

    capture.output(
      qh,
      file = file.path(
        DIR_MVMR,
        paste0(
          "QHET__",
          tag,
          "__rho_",
          BMI_OSA_PHENO_COR,
          ".txt"
        )
      )
    )

    qhet_results[[oy]] <- qh
  }

} else {

  writeLines(
    c(
      "qhet_mvmr was NOT run.",
      "",
      "Reason:",
      "A defensible BMI–OSA phenotypic correlation was not supplied.",
      "",
      "Do not substitute LDSC genetic correlation for phenotypic correlation.",
      "If you later obtain an appropriate phenotypic correlation, set:",
      "BMI_OSA_PHENO_COR <- <value>",
      "and rerun STEP3C."
    ),
    file.path(
      DIR_REP,
      "STEP3C_QHET_NOT_RUN.txt"
    )
  )
}


# ==============================================================================
# 19. CONDITIONAL-F QC FIGURE
# ==============================================================================

plot_qc <- copy(grid_qc[status == "OK"])

plot_qc[
  ,
  combo_label := paste0(
    "BMI ", BMI_tier,
    "\nOSA ", OSA_tier
  )
]

# one figure per outcome, per user's single-panel preference
for (oy in c("HFpEF", "AF")) {

  d <- plot_qc[outcome == oy]

  d <- melt(
    d,
    id.vars = c(
      "combo_id",
      "combo_label",
      "relax_score"
    ),
    measure.vars = c(
      "F_BMI",
      "F_OSA"
    ),
    variable.name = "Exposure",
    value.name = "Conditional_F"
  )

  d[
    ,
    Exposure := fifelse(
      Exposure == "F_BMI",
      "BMI",
      "OSA"
    )
  ]

  # order by relaxation score then label
  ord <- unique(
    plot_qc[
      outcome == oy
    ][
      order(
        relax_score,
        BMI_rank,
        OSA_rank
      ),
      combo_label
    ]
  )

  d[
    ,
    combo_label := factor(
      combo_label,
      levels = ord
    )
  ]

  p <- ggplot(
    d,
    aes(
      x = combo_label,
      y = Conditional_F,
      group = Exposure,
      shape = Exposure
    )
  ) +
    geom_hline(
      yintercept = 10,
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_line(
      linewidth = 0.65
    ) +
    geom_point(
      size = 2.6
    ) +
    labs(
      x = NULL,
      y = "Conditional F-statistic"
    ) +
    theme_classic(
      base_size = 10.5,
      base_family = BASE_FAMILY
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none"
    )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(
        "STEP3C_ConditionalF_grid_",
        oy,
        ".tiff"
      )
    ),
    p,
    width = 8.2,
    height = 4.3,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  ggsave(
    file.path(
      DIR_FIG,
      paste0(
        "STEP3C_ConditionalF_grid_",
        oy,
        ".pdf"
      )
    ),
    p,
    width = 8.2,
    height = 4.3,
    units = "in",
    device = cairo_pdf
  )
}

writeLines(
  c(
    "Conditional F figure legend",
    "",
    "Conditional F-statistics for BMI and OSA across the pre-specified MVMR instrument-threshold grid.",
    "The horizontal dashed line denotes the conventional F=10 threshold.",
    "Threshold selection was based solely on instrument strength and not on outcome effect estimates or P-values."
  ),
  file.path(
    DIR_FIG,
    "STEP3C_ConditionalF_grid_legend.txt"
  )
)


# ==============================================================================
# 20. SELECTED MVMR FOREST FIGURE
# ==============================================================================

valid_plot <- selected_effects[
  analysis_class !=
    "INVALID_WEAK_INSTRUMENT_DIAGNOSTIC"
]

if (nrow(valid_plot) > 0) {

  COL_BONF <- "#C65353"
  COL_NOM <- "#3F78A8"
  COL_NS <- "#8B9098"

  d <- copy(valid_plot)

  d[
    ,
    pair := paste0(
      exposure,
      " → ",
      outcome
    )
  ]

  d[
    ,
    pair := factor(
      pair,
      levels = rev(c(
        "BMI → HFpEF",
        "OSA → HFpEF",
        "BMI → AF",
        "OSA → AF"
      ))
    )
  ]

  d[
    ,
    point_col := fifelse(
      significance_class ==
        "Bonferroni significant",
      COL_BONF,
      fifelse(
        significance_class ==
          "Nominal",
        COL_NOM,
        COL_NS
      )
    )
  ]

  p <- ggplot(
    d,
    aes(
      y = pair,
      x = OR
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_errorbarh(
      aes(
        xmin = OR_low,
        xmax = OR_high,
        colour = point_col
      ),
      height = 0,
      linewidth = 0.7
    ) +
    geom_point(
      aes(colour = point_col),
      shape = 15,
      size = 3.2
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

  ggsave(
    file.path(
      DIR_FIG,
      "STEP3C_StrengthOptimized_MVMR.tiff"
    ),
    p,
    width = 6.5,
    height = 4.0,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw"
  )

  ggsave(
    file.path(
      DIR_FIG,
      "STEP3C_StrengthOptimized_MVMR.pdf"
    ),
    p,
    width = 6.5,
    height = 4.0,
    units = "in",
    device = cairo_pdf
  )

  fwrite(
    d,
    file.path(
      DIR_FIG,
      "STEP3C_StrengthOptimized_MVMR_source_data.csv"
    )
  )

  writeLines(
    c(
      "STEP3C strength-optimized MVMR legend",
      "",
      "Squares indicate IVW multivariable MR direct-effect estimates; horizontal lines indicate 95% confidence intervals.",
      "BMI and OSA were modeled simultaneously.",
      "The displayed tier was selected using only conditional instrument strength: both conditional F-statistics had to be >=10, and the least relaxed qualifying threshold combination was selected.",
      "If any selected instruments used P>5×10^-8, this figure represents exploratory/sensitivity MVMR and does not replace the strict genome-wide-significant primary MVMR.",
      "Red/pink: Bonferroni significant; blue: nominal P<0.05; gray: non-significant."
    ),
    file.path(
      DIR_FIG,
      "STEP3C_StrengthOptimized_MVMR_legend.txt"
    )
  )
}


# ==============================================================================
# 21. FINAL DECISION REPORT
# ==============================================================================

decision_lines <- c(
  "STEP 3C FINAL DECISION",
  "",
  paste0(
    "Conventional conditional F pass threshold: ",
    F_PASS
  ),
  "",
  "Selected tier(s):"
)

for (oy in c("HFpEF", "AF")) {

  s <- selected_qc[outcome == oy]

  decision_lines <- c(
    decision_lines,
    paste0(
      oy,
      ": ",
      s$combo_id,
      " | F_BMI=",
      signif(s$F_BMI, 4),
      " | F_OSA=",
      signif(s$F_OSA, 4),
      " | minF=",
      signif(s$min_conditional_F, 4),
      " | ",
      s$selection_status
    )
  )
}

decision_lines <- c(
  decision_lines,
  "",
  "Interpretation rule:",
  "- PASS_BOTH_CONDITIONAL_F_GE10: conventional IVW MVMR instrument-strength QC passed.",
  "- If selected tier contains P>5e-8 IVs, classify it as exploratory/sensitivity MVMR.",
  "- NO_TIER_PASSES_BOTH_F_GE10: do not interpret conventional IVW MVMR as a valid direct causal estimate; obtain stronger exposure GWAS/instruments or use a justified weak-instrument robust method."
)

writeLines(
  decision_lines,
  file.path(
    DIR_REP,
    "STEP3C_FINAL_DECISION.txt"
  )
)


# ==============================================================================
# 22. SESSION INFO + CLOSE
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_REP,
    "STEP3C_sessionInfo.txt"
  )
)

dbDisconnect(
  con,
  shutdown = TRUE
)

log_msg("============================================================")
log_msg("STEP 3C COMPLETE")
log_msg("============================================================")
log_msg("Output: ", OUT_DIR)
log_msg("")
log_msg("Please send back:")
log_msg("1) 05_reports/STEP3C_MVMR_grid_strength_QC.csv")
log_msg("2) 05_reports/STEP3C_SELECTED_strength_optimized_tiers.csv")
log_msg("3) 03_mvmr_results/STEP3C_SELECTED_MVMR_effects.csv")
log_msg("4) 05_reports/STEP3C_FINAL_DECISION.txt")
log_msg("5) 05_reports/STEP3C_master_console.log")
log_msg("6) 04_figures/STEP3C_ConditionalF_grid_HFpEF.tiff")
log_msg("7) 04_figures/STEP3C_ConditionalF_grid_AF.tiff")
