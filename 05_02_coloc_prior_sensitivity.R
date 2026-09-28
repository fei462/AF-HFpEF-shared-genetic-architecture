# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP8B_COLOC_PRIOR_SENSITIVITY.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP8B — PRIOR SENSITIVITY FOR TARGETED coloc.abf RESULTS
# Project: AF / HFpEF / BMI / OSA
#
# Purpose:
#   Re-run the 11 pre-specified STEP8A pair-region coloc tests
#   across a grid of p12 priors while keeping p1=p2=1e-4.
#
# Why:
#   coloc is Bayesian and conclusions may depend on the prior
#   probability that a SNP is associated with BOTH traits (p12).
#
# Primary sensitivity range:
#   p12 = 1e-6 ... 1e-5
#   (more conservative through the default)
#
# Extended exploratory range:
#   p12 = 1e-7 ... 1e-4
#
# We DO NOT change regions or SNPs based on the observed results.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")
COLOC_ROOT <- file.path(DATA_ROOT, "STEP8_COLOC")

REGION_DIR <- file.path(
  COLOC_ROOT,
  "01_regional_sumstats"
)

TABLE_DIR <- file.path(
  COLOC_ROOT,
  "03_tables"
)

SENS_DIR <- file.path(
  COLOC_ROOT,
  "05_prior_sensitivity"
)

LOG_DIR <- file.path(
  COLOC_ROOT,
  "04_logs"
)

dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

MANIFEST_FILE <- file.path(
  TABLE_DIR,
  "STEP8A_test_manifest.csv"
)

INFO_FILE <- file.path(
  LAVA_ROOT,
  "02_input",
  "input.info.txt"
)

TRAIT_REGION_FILES <- c(
  AF = file.path(REGION_DIR, "STEP8A_AF_target_regions.tsv.gz"),
  HFpEF = file.path(REGION_DIR, "STEP8A_HFpEF_target_regions.tsv.gz"),
  BMI = file.path(REGION_DIR, "STEP8A_BMI_target_regions.tsv.gz"),
  OSA = file.path(REGION_DIR, "STEP8A_OSA_target_regions.tsv.gz")
)

required_files <- c(
  MANIFEST_FILE,
  INFO_FILE,
  TRAIT_REGION_FILES
)

if (any(!file.exists(required_files))) {
  stop(
    "Missing required STEP8A/STEP7 input file(s):\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    )
  )
}

# ============================================================
# 1. PACKAGES
# ============================================================

for (p in c("data.table", "coloc")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(data.table)
library(coloc)

cat(
  "\ncoloc version: ",
  as.character(packageVersion("coloc")),
  "\n",
  sep = ""
)

# ============================================================
# 2. LOAD MANIFEST + TRAIT METADATA
# ============================================================

TESTS <- fread(MANIFEST_FILE)

required_manifest_cols <- c(
  "region",
  "trait1",
  "trait2",
  "upstream_basis",
  "test_id"
)

if (!all(required_manifest_cols %in% names(TESTS))) {
  stop(
    "STEP8A manifest lacks columns: ",
    paste(
      setdiff(required_manifest_cols, names(TESTS)),
      collapse = ", "
    )
  )
}

INFO <- fread(
  INFO_FILE,
  na.strings = c("NA", "")
)

TRAIT_TYPE <- c(
  AF = "cc",
  HFpEF = "cc",
  BMI = "quant",
  OSA = "cc"
)

CASE_PROP <- setNames(
  rep(NA_real_, length(TRAIT_TYPE)),
  names(TRAIT_TYPE)
)

for (tr in names(TRAIT_TYPE)) {

  if (TRAIT_TYPE[[tr]] == "cc") {

    ca <- as.numeric(
      INFO[phenotype == tr, cases]
    )

    co <- as.numeric(
      INFO[phenotype == tr, controls]
    )

    if (
      length(ca) != 1 ||
      length(co) != 1 ||
      !is.finite(ca) ||
      !is.finite(co) ||
      ca <= 0 ||
      co <= 0
    ) {
      stop(
        "Invalid case/control metadata for ",
        tr
      )
    }

    CASE_PROP[[tr]] <- ca / (ca + co)
  }
}

# ============================================================
# 3. LOAD REGIONAL SUMSTATS
# ============================================================

REGIONAL <- lapply(
  TRAIT_REGION_FILES,
  fread,
  showProgress = FALSE
)

names(REGIONAL) <- names(TRAIT_REGION_FILES)

required_region_cols <- c(
  "region",
  "SNP",
  "POS",
  "MAF",
  "N",
  "Z",
  "pvalue"
)

for (tr in names(REGIONAL)) {

  x <- REGIONAL[[tr]]

  if (!all(required_region_cols %in% names(x))) {
    stop(
      tr,
      " regional file lacks: ",
      paste(
        setdiff(required_region_cols, names(x)),
        collapse = ", "
      )
    )
  }

  x[, SNP := as.character(SNP)]
  x[, region := as.character(region)]

  REGIONAL[[tr]] <- x
}

# ============================================================
# 4. PRIORS
# ============================================================

P1 <- 1e-4
P2 <- 1e-4

# Dense log-spaced grid for inspection
P12_GRID <- sort(unique(c(
  10^seq(-7, -4, length.out = 25),
  1e-6,
  3e-6,
  1e-5,
  3e-5,
  1e-4
)))

# Primary sensitivity interval:
# conservative p12 through default p12
PRIMARY_MIN <- 1e-6
PRIMARY_MAX <- 1e-5

# Descriptive decision rule, NOT a universal coloc standard.
# It mirrors the common idea that H4 should dominate H3.
passes_rule <- function(h3, h4) {
  is.finite(h3) &&
    is.finite(h4) &&
    h4 > 0.80 &&
    h4 > 3 * h3
}

# ============================================================
# 5. COLOC DATASET BUILDER
# ============================================================

build_dataset <- function(tr, x) {

  n_med <- round(
    median(
      x$N,
      na.rm = TRUE
    )
  )

  if (!is.finite(n_med) || n_med <= 0) {
    stop(
      "Invalid regional N for ",
      tr
    )
  }

  d <- list(
    pvalues = x$pvalue,
    MAF = x$MAF,
    N = n_med,
    type = TRAIT_TYPE[[tr]],
    snp = x$SNP,
    position = x$POS
  )

  if (TRAIT_TYPE[[tr]] == "cc") {
    d$s <- CASE_PROP[[tr]]
  }

  coloc::check_dataset(
    d,
    warn.minp = 1
  )

  d
}

# ============================================================
# 6. RUN SENSITIVITY GRID
# ============================================================

GRID_LIST <- list()

for (i in seq_len(nrow(TESTS))) {

  tt <- TESTS[i]

  test_id_i <- as.character(tt$test_id)
  region_i <- as.character(tt$region)
  trait1_i <- as.character(tt$trait1)
  trait2_i <- as.character(tt$trait2)
  basis_i <- as.character(tt$upstream_basis)

  cat(
    "\n[",
    i,
    "/",
    nrow(TESTS),
    "] ",
    test_id_i,
    "\n",
    sep = ""
  )

  d1 <- REGIONAL[[trait1_i]][
    region == region_i
  ]

  d2 <- REGIONAL[[trait2_i]][
    region == region_i
  ]

  # Pairwise intersection, preserving the same MAF and region
  m <- merge(
    d1[
      ,
      .(
        SNP,
        POS,
        MAF,
        N1 = N,
        P1v = pvalue
      )
    ],
    d2[
      ,
      .(
        SNP,
        N2 = N,
        P2v = pvalue
      )
    ],
    by = "SNP",
    all = FALSE
  )

  m <- m[
    is.finite(MAF) &
    MAF > 0 &
    MAF <= 0.5 &
    is.finite(P1v) &
    P1v > 0 &
    P1v <= 1 &
    is.finite(P2v) &
    P2v > 0 &
    P2v <= 1 &
    is.finite(N1) &
    N1 > 0 &
    is.finite(N2) &
    N2 > 0
  ]

  m <- unique(
    m,
    by = "SNP"
  )

  if (nrow(m) < 100) {
    stop(
      test_id_i,
      " has too few pairwise common SNPs: ",
      nrow(m)
    )
  }

  dat1 <- data.table(
    SNP = m$SNP,
    POS = m$POS,
    MAF = m$MAF,
    N = m$N1,
    pvalue = m$P1v
  )

  dat2 <- data.table(
    SNP = m$SNP,
    POS = m$POS,
    MAF = m$MAF,
    N = m$N2,
    pvalue = m$P2v
  )

  D1 <- build_dataset(
    trait1_i,
    dat1
  )

  D2 <- build_dataset(
    trait2_i,
    dat2
  )

  for (p12_i in P12_GRID) {

    res <- coloc::coloc.abf(
      dataset1 = D1,
      dataset2 = D2,
      p1 = P1,
      p2 = P2,
      p12 = p12_i
    )

    sm <- res$summary

    h0 <- as.numeric(sm["PP.H0.abf"])
    h1 <- as.numeric(sm["PP.H1.abf"])
    h2 <- as.numeric(sm["PP.H2.abf"])
    h3 <- as.numeric(sm["PP.H3.abf"])
    h4 <- as.numeric(sm["PP.H4.abf"])

    h4_cond <- if (
      is.finite(h3) &&
      is.finite(h4) &&
      (h3 + h4) > 0
    ) {
      h4 / (h3 + h4)
    } else {
      NA_real_
    }

    GRID_LIST[[length(GRID_LIST) + 1L]] <- data.table(
      test_id = test_id_i,
      region = region_i,
      trait1 = trait1_i,
      trait2 = trait2_i,
      upstream_basis = basis_i,
      nsnps = nrow(m),
      p1 = P1,
      p2 = P2,
      p12 = p12_i,
      PP.H0 = h0,
      PP.H1 = h1,
      PP.H2 = h2,
      PP.H3 = h3,
      PP.H4 = h4,
      H4_over_H3H4 = h4_cond,
      pass_rule = passes_rule(
        h3,
        h4
      )
    )
  }
}

GRID <- rbindlist(
  GRID_LIST,
  fill = TRUE
)

GRID_FILE <- file.path(
  SENS_DIR,
  "STEP8B_prior_sensitivity_grid.csv"
)

fwrite(
  GRID,
  GRID_FILE
)

# ============================================================
# 7. SUMMARISE ROBUSTNESS
# ============================================================

DEFAULT <- GRID[
  abs(log10(p12) - log10(1e-5)) < 1e-12
]

P12_1E6 <- GRID[
  abs(log10(p12) - log10(1e-6)) < 1e-12
]

PRIMARY <- GRID[
  p12 >= PRIMARY_MIN &
  p12 <= PRIMARY_MAX
]

SUMMARY <- PRIMARY[
  ,
  .(
    n_primary_grid = .N,
    min_PP_H4_primary = min(
      PP.H4,
      na.rm = TRUE
    ),
    max_PP_H4_primary = max(
      PP.H4,
      na.rm = TRUE
    ),
    min_H4_over_H3H4_primary = min(
      H4_over_H3H4,
      na.rm = TRUE
    ),
    max_H4_over_H3H4_primary = max(
      H4_over_H3H4,
      na.rm = TRUE
    ),
    all_primary_pass_rule = all(
      pass_rule
    ),
    proportion_primary_pass_rule = mean(
      pass_rule
    )
  ),
  by = .(
    test_id,
    region,
    trait1,
    trait2,
    upstream_basis
  )
]

SUMMARY <- merge(
  SUMMARY,
  DEFAULT[
    ,
    .(
      test_id,
      PP_H4_default = PP.H4,
      PP_H3_default = PP.H3,
      H4cond_default = H4_over_H3H4,
      pass_default_rule = pass_rule
    )
  ],
  by = "test_id",
  all.x = TRUE
)

SUMMARY <- merge(
  SUMMARY,
  P12_1E6[
    ,
    .(
      test_id,
      PP_H4_p12_1e6 = PP.H4,
      PP_H3_p12_1e6 = PP.H3,
      H4cond_p12_1e6 = H4_over_H3H4,
      pass_p12_1e6_rule = pass_rule
    )
  ],
  by = "test_id",
  all.x = TRUE
)

# Smallest p12 at which the descriptive rule is satisfied
P12_THRESHOLD <- GRID[
  pass_rule == TRUE,
  .(
    min_p12_passing_rule = min(
      p12
    )
  ),
  by = test_id
]

SUMMARY <- merge(
  SUMMARY,
  P12_THRESHOLD,
  by = "test_id",
  all.x = TRUE
)

SUMMARY[
  ,
  sensitivity_class :=
    fifelse(
      all_primary_pass_rule == TRUE,
      "Robust across p12=1e-6 to 1e-5",
      fifelse(
        pass_default_rule == TRUE,
        "Supports H4 at default prior but not across full conservative range",
        "Does not meet H4-dominance rule at default prior"
      )
    )
]

setorder(
  SUMMARY,
  -PP_H4_default
)

SUMMARY_FILE <- file.path(
  TABLE_DIR,
  "STEP8B_prior_sensitivity_summary.csv"
)

fwrite(
  SUMMARY,
  SUMMARY_FILE
)

# ============================================================
# 8. SIMPLE PRIOR-SENSITIVITY FIGURES
#
# One figure per test, no combined multi-panel output.
# ============================================================

FIG_DIR <- file.path(
  SENS_DIR,
  "figures"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

for (tid in unique(GRID$test_id)) {

  x <- GRID[
    test_id == tid
  ]

  png(
    filename = file.path(
      FIG_DIR,
      paste0(
        "STEP8B_",
        gsub(
          "[^A-Za-z0-9_]+",
          "_",
          tid
        ),
        "_prior_sensitivity.png"
      )
    ),
    width = 1800,
    height = 1400,
    res = 220
  )

  plot(
    x$p12,
    x$PP.H4,
    log = "x",
    type = "l",
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Prior probability p12",
    ylab = "Posterior probability",
    main = tid
  )

  lines(
    x$p12,
    x$PP.H3,
    lwd = 2,
    lty = 2
  )

  abline(
    v = 1e-5,
    lty = 3
  )

  abline(
    h = 0.8,
    lty = 3
  )

  dev.off()
}

# ============================================================
# 9. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP8B_sessionInfo.txt"
  )
)

# ============================================================
# 10. FINAL REPORT
# ============================================================

cat("\n====================================================\n")
cat("STEP8B COMPLETE — coloc PRIOR SENSITIVITY\n")
cat("====================================================\n")

cat("\nSensitivity summary:\n")
print(
  SUMMARY[
    ,
    .(
      test_id,
      PP_H4_default,
      PP_H4_p12_1e6,
      min_PP_H4_primary,
      H4cond_default,
      all_primary_pass_rule,
      min_p12_passing_rule,
      sensitivity_class
    )
  ]
)

cat("\nUPLOAD THESE 2 FILES:\n")
cat("1) ", SUMMARY_FILE, "\n", sep = "")
cat("2) ", GRID_FILE, "\n", sep = "")

cat(
  "\nAfter STEP8B is checked, the next step is coloc.susie ",
  "for strong/ambiguous regions using UKB v1.1 LD.\n",
  sep = ""
)

cat("====================================================\n")
