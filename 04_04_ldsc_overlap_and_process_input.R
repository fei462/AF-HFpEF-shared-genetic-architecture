# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP7B2_V4_STABLE_LDSC_OVERLAP_PROCESS_INPUT.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP7B2_V4 — STABLE LDSC SAMPLE OVERLAP + LAVA process.input QC
# Project: AF / HFpEF / BMI / OSA
#
# WHY THIS VERSION:
#   Previous versions failed because of log/CSV parsing and one
#   incorrect data.table ordering expression.
#
#   This version does NOT parse raw LDSC logs or result CSVs.
#   It uses the already completed and manually verified STEP1D
#   LDSC intercept values directly, then performs only:
#       1) sample-overlap matrix construction
#       2) input.info normalization
#       3) LAVA::process.input()
#       4) hard QC export
#
#   No missing values are set to zero.
#   No nearPD correction is used.
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PATHS
# ============================================================

DATA_ROOT <- "D:/A/data"
LAVA_ROOT <- file.path(DATA_ROOT, "STEP7_LAVA")

IN_DIR  <- file.path(LAVA_ROOT, "02_input")
TAB_DIR <- file.path(LAVA_ROOT, "07_tables")
LOG_DIR <- file.path(LAVA_ROOT, "08_logs")

dir.create(IN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

INFO_FILE <- file.path(IN_DIR, "input.info.txt")

REF_PREFIX <- file.path(
  LAVA_ROOT,
  "00_reference",
  "UKB_v1.1",
  "lava-ukb-v1.1"
)

PHEN <- c("AF", "HFpEF", "BMI", "OSA")

# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
library(data.table)

# ============================================================
# 2. FIXED, VERIFIED STEP1D LDSC VALUES
#
# Diagonal:
#   univariate LDSC intercepts
#
# Off-diagonal:
#   cross-trait LDSC Genetic Covariance intercepts (gcov_int)
# ============================================================

H2_INTERCEPT <- data.table(
  phenotype = c("AF", "HFpEF", "BMI", "OSA"),
  intercept = c(
    1.5097,
    0.9832,
    0.7339,
    1.1781
  )
)

GCOV_INTERCEPT <- data.table(
  trait1 = c(
    "AF",
    "AF",
    "AF",
    "HFpEF",
    "HFpEF",
    "BMI"
  ),
  trait2 = c(
    "HFpEF",
    "BMI",
    "OSA",
    "BMI",
    "OSA",
    "OSA"
  ),
  gcov_int = c(
    0.0636,
    0.0486,
    0.1072,
    0.0035,
    0.0031,
    0.0371
  )
)

# Save exactly what is being used
H2_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_univariate_LDSC_intercepts.csv"
)

GCOV_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_pairwise_gcov_intercepts.csv"
)

fwrite(H2_INTERCEPT, H2_FILE)
fwrite(GCOV_INTERCEPT, GCOV_FILE)

# ============================================================
# 3. BASIC VALIDATION OF FIXED LDSC VALUES
# ============================================================

stopifnot(
  nrow(H2_INTERCEPT) == 4,
  nrow(GCOV_INTERCEPT) == 6,
  setequal(H2_INTERCEPT$phenotype, PHEN),
  all(is.finite(H2_INTERCEPT$intercept)),
  all(H2_INTERCEPT$intercept > 0),
  all(is.finite(GCOV_INTERCEPT$gcov_int))
)

expected_pairs <- c(
  "AF__HFpEF",
  "AF__BMI",
  "AF__OSA",
  "HFpEF__BMI",
  "HFpEF__OSA",
  "BMI__OSA"
)

observed_pairs <- paste(
  GCOV_INTERCEPT$trait1,
  GCOV_INTERCEPT$trait2,
  sep = "__"
)

stopifnot(setequal(expected_pairs, observed_pairs))

# ============================================================
# 4. BUILD LDSC SAMPLING COVARIANCE MATRIX
# ============================================================

S_COV <- matrix(
  NA_real_,
  nrow = length(PHEN),
  ncol = length(PHEN),
  dimnames = list(PHEN, PHEN)
)

# diagonal
for (tr in PHEN) {
  S_COV[tr, tr] <- H2_INTERCEPT[
    phenotype == tr,
    intercept
  ]
}

# off-diagonal
for (i in seq_len(nrow(GCOV_INTERCEPT))) {

  a <- GCOV_INTERCEPT$trait1[i]
  b <- GCOV_INTERCEPT$trait2[i]
  v <- GCOV_INTERCEPT$gcov_int[i]

  S_COV[a, b] <- v
  S_COV[b, a] <- v
}

if (any(!is.finite(S_COV))) {
  stop("ERROR: Sampling covariance matrix contains NA/Inf.")
}

RAW_COV_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_LDSC_sampling_covariance_raw.txt"
)

write.table(
  S_COV,
  RAW_COV_FILE,
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

cat("\n================ RAW LDSC COVARIANCE ================\n")
print(round(S_COV, 5))

# ============================================================
# 5. CONVERT TO SAMPLING CORRELATION MATRIX
#
# This follows the official LAVA sample-overlap tutorial.
# ============================================================

S_COR <- round(cov2cor(S_COV), 5)

if (any(!is.finite(S_COR))) {
  stop("ERROR: Sampling correlation matrix contains NA/Inf.")
}

if (max(abs(S_COR - t(S_COR))) > 1e-12) {
  stop("ERROR: Sampling correlation matrix is not symmetric.")
}

if (any(abs(S_COR) > 1.00001)) {
  stop("ERROR: Sampling correlation outside [-1,1].")
}

EIG <- eigen(
  S_COR,
  symmetric = TRUE,
  only.values = TRUE
)$values

if (min(EIG) < -1e-8) {
  stop(
    "ERROR: Sampling correlation matrix is not positive semidefinite. ",
    "Minimum eigenvalue = ",
    min(EIG)
  )
}

OVERLAP_FILE <- file.path(
  IN_DIR,
  "sample.overlap.txt"
)

write.table(
  S_COR,
  OVERLAP_FILE,
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

MATRIX_QC <- data.table(
  metric = c(
    "min_eigenvalue",
    "max_eigenvalue",
    "max_abs_offdiag",
    "symmetric",
    "diag_all_1"
  ),
  value = c(
    min(EIG),
    max(EIG),
    max(abs(S_COR[row(S_COR) != col(S_COR)])),
    1,
    as.numeric(all(diag(S_COR) == 1))
  )
)

MATRIX_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_sample_overlap_QC.csv"
)

fwrite(MATRIX_QC, MATRIX_QC_FILE)

cat("\n================ SAMPLE OVERLAP MATRIX ================\n")
print(S_COR)

# Expected approximate matrix; hard guard against accidental edits
EXPECTED_COR <- matrix(
  c(
    1.00000, 0.05220, 0.04617, 0.08038,
    0.05220, 1.00000, 0.00412, 0.00288,
    0.04617, 0.00412, 1.00000, 0.03990,
    0.08038, 0.00288, 0.03990, 1.00000
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(PHEN, PHEN)
)

if (max(abs(S_COR - EXPECTED_COR)) > 0.00002) {
  stop(
    "ERROR: Constructed sample-overlap matrix differs from ",
    "the verified expected matrix."
  )
}

# ============================================================
# 6. VALIDATE / NORMALIZE input.info
#
# Current LAVA source treats a continuous phenotype as:
# cases=1, controls=0 -> prop_cases=1 -> binary=FALSE.
#
# We therefore normalize BMI to 1/0 here.
# ============================================================

if (!file.exists(INFO_FILE)) {
  stop(
    "ERROR: Missing input.info.txt:\n",
    INFO_FILE,
    "\nRun STEP7B_V2 first."
  )
}

INFO <- fread(
  INFO_FILE,
  na.strings = c("NA", "")
)

required_info_cols <- c(
  "phenotype",
  "cases",
  "controls",
  "filename"
)

if (!all(required_info_cols %in% names(INFO))) {
  stop(
    "ERROR: input.info.txt must contain columns: ",
    paste(required_info_cols, collapse = ", ")
  )
}

if (!setequal(INFO$phenotype, PHEN)) {
  stop(
    "ERROR: input.info phenotypes are not exactly AF, HFpEF, BMI, OSA."
  )
}

# Reorder safely using base matching
INFO <- INFO[match(PHEN, phenotype)]

# Normalize BMI to official continuous-trait sentinel
INFO[
  phenotype == "BMI",
  `:=`(
    cases = 1,
    controls = 0
  )
]

# Hard checks for binary phenotypes
expected_meta <- data.table(
  phenotype = c("AF", "OSA"),
  cases = c(228926, 38998),
  controls = c(1611415, 336659)
)

for (i in seq_len(nrow(expected_meta))) {

  tr <- expected_meta$phenotype[i]

  obs_case <- INFO[phenotype == tr, as.numeric(cases)]
  obs_ctrl <- INFO[phenotype == tr, as.numeric(controls)]

  if (
    length(obs_case) != 1 ||
    length(obs_ctrl) != 1 ||
    obs_case != expected_meta$cases[i] ||
    obs_ctrl != expected_meta$controls[i]
  ) {
    stop(
      "ERROR: input.info case/control metadata mismatch for ",
      tr
    )
  }
}

# HFpEF must be binary with positive case/control counts
hf_case <- INFO[phenotype == "HFpEF", as.numeric(cases)]
hf_ctrl <- INFO[phenotype == "HFpEF", as.numeric(controls)]

if (
  length(hf_case) != 1 ||
  length(hf_ctrl) != 1 ||
  !is.finite(hf_case) ||
  !is.finite(hf_ctrl) ||
  hf_case <= 0 ||
  hf_ctrl <= 0
) {
  stop("ERROR: Invalid HFpEF case/control counts in input.info.")
}

# All summary-stat files must exist
if (!all(file.exists(INFO$filename))) {
  print(INFO[!file.exists(filename)])
  stop("ERROR: One or more LAVA summary-stat files do not exist.")
}

# Overwrite normalized input.info
fwrite(
  INFO,
  INFO_FILE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

INFO_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_input_info_used.csv"
)

fwrite(INFO, INFO_QC_FILE)

cat("\n================ input.info USED ================\n")
print(INFO)

# ============================================================
# 7. UKB v1.1 REFERENCE HARD CHECK
# ============================================================

REF_QC <- rbindlist(
  lapply(1:22, function(chr) {

    f_info <- paste0(
      REF_PREFIX,
      "_chr",
      chr,
      ".info"
    )

    f_bcor <- paste0(
      REF_PREFIX,
      "_chr",
      chr,
      ".bcor"
    )

    data.table(
      chr = chr,
      info_exists = file.exists(f_info),
      bcor_exists = file.exists(f_bcor),
      info_file = f_info,
      bcor_file = f_bcor
    )
  })
)

REF_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_reference_QC.csv"
)

fwrite(REF_QC, REF_QC_FILE)

if (
  any(!REF_QC$info_exists) ||
  any(!REF_QC$bcor_exists)
) {
  print(REF_QC[!info_exists | !bcor_exists])
  stop(
    "ERROR: UKB v1.1 reference incomplete or misnamed."
  )
}

# ============================================================
# 8. LAVA PACKAGE CHECK
# ============================================================

if (!requireNamespace("LAVA", quietly = TRUE)) {

  stop(
    "ERROR: LAVA is not installed.\n",
    "Install it first with:\n",
    "if (!require('remotes')) install.packages('remotes')\n",
    "remotes::install_github('josefin-werme/LAVA')\n",
    "Then restart R and rerun this script."
  )
}

library(LAVA)

LAVA_VERSION <- as.character(
  packageVersion("LAVA")
)

cat("\nLAVA version:", LAVA_VERSION, "\n")

# ============================================================
# 9. RUN process.input()
# ============================================================

cat("\n====================================================\n")
cat("RUNNING LAVA::process.input()\n")
cat("====================================================\n")

input <- process.input(
  input.info.file = INFO_FILE,
  sample.overlap.file = OVERLAP_FILE,
  ref.prefix = REF_PREFIX,
  phenos = PHEN
)

# ============================================================
# 10. HARD CHECK process.input() RETURN OBJECT
# ============================================================

if (!is.environment(input)) {
  stop(
    "ERROR: process.input() did not return an environment."
  )
}

ENV_NAMES <- ls(input)

required_objects <- c(
  "info",
  "P",
  "sample.overlap",
  "sum.stats",
  "analysis.snps",
  "reference"
)

missing_objects <- setdiff(
  required_objects,
  ENV_NAMES
)

if (length(missing_objects)) {
  stop(
    "ERROR: process.input() missing expected objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

ENV_NAMES_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_process_input_environment_names.txt"
)

writeLines(
  ENV_NAMES,
  ENV_NAMES_FILE
)

# ============================================================
# 11. VERIFY PHENOTYPE CLASSIFICATION
# ============================================================

INFO_PROCESSED <- as.data.table(
  input$info
)

expected_binary <- c(
  AF = TRUE,
  HFpEF = TRUE,
  BMI = FALSE,
  OSA = TRUE
)

for (tr in PHEN) {

  obs <- INFO_PROCESSED[
    phenotype == tr,
    binary
  ]

  if (
    length(obs) != 1 ||
    is.na(obs) ||
    as.logical(obs) != expected_binary[[tr]]
  ) {
    stop(
      "ERROR: LAVA phenotype classification incorrect for ",
      tr
    )
  }
}

PROCESSED_INFO_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_process_input_info_QC.csv"
)

fwrite(
  INFO_PROCESSED,
  PROCESSED_INFO_FILE
)

# ============================================================
# 12. COMMON SNP / ALIGNMENT QC
# ============================================================

N_ANALYSIS <- length(
  input$analysis.snps
)

N_UNALIGNABLE <- if (
  "unalignable.snps" %in% ENV_NAMES
) {
  length(input$unalignable.snps)
} else {
  NA_integer_
}

if (N_ANALYSIS < 100000) {
  stop(
    "ERROR: Only ",
    N_ANALYSIS,
    " common analysis SNPs remain; likely an ID/alignment problem."
  )
}

COMMON_QC <- data.table(
  metric = c(
    "LAVA_version",
    "n_phenotypes",
    "n_analysis_snps",
    "n_unalignable_snps",
    "process_input_completed"
  ),
  value = c(
    LAVA_VERSION,
    input$P,
    N_ANALYSIS,
    N_UNALIGNABLE,
    1
  )
)

COMMON_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_process_input_common_QC.csv"
)

fwrite(
  COMMON_QC,
  COMMON_QC_FILE
)

# ============================================================
# 13. POST-HARMONIZATION SUMSTATS QC
# ============================================================

POST_QC <- rbindlist(
  lapply(PHEN, function(tr) {

    z <- input$sum.stats[[tr]]

    required <- c(
      "SNP",
      "A1",
      "A2",
      "STAT",
      "N"
    )

    if (!all(required %in% names(z))) {
      stop(
        "ERROR: Processed sumstats for ",
        tr,
        " lack required columns."
      )
    }

    data.table(
      phenotype = tr,
      rows_after_process_input = nrow(z),
      median_N = median(z$N, na.rm = TRUE),
      min_N = min(z$N, na.rm = TRUE),
      max_N = max(z$N, na.rm = TRUE),
      mean_STAT2 = mean(z$STAT^2, na.rm = TRUE),
      max_abs_STAT = max(abs(z$STAT), na.rm = TRUE),
      missing_STAT = sum(!is.finite(z$STAT)),
      missing_N = sum(!is.finite(z$N))
    )
  })
)

if (any(POST_QC$missing_STAT > 0) ||
    any(POST_QC$missing_N > 0)) {
  print(POST_QC)
  stop(
    "ERROR: Missing/non-finite processed STAT or N after process.input()."
  )
}

POST_QC_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_process_input_sumstats_QC.csv"
)

fwrite(
  POST_QC,
  POST_QC_FILE
)

# ============================================================
# 14. VERIFY SAMPLE OVERLAP AS SEEN BY LAVA
# ============================================================

LAVA_OVERLAP <- as.matrix(
  input$sample.overlap
)

if (
  !all(dim(LAVA_OVERLAP) == c(4, 4)) ||
  max(abs(LAVA_OVERLAP - S_COR)) > 0.00002
) {
  stop(
    "ERROR: LAVA-processed sample overlap differs from expected matrix."
  )
}

LAVA_OVERLAP_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_process_input_overlap_matrix.txt"
)

write.table(
  LAVA_OVERLAP,
  LAVA_OVERLAP_FILE,
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

# ============================================================
# 15. SAVE INPUT OBJECT FOR STEP7C
# ============================================================

INPUT_RDATA <- file.path(
  IN_DIR,
  "STEP7B2_V4_LAVA_input.RData"
)

save(
  input,
  file = INPUT_RDATA,
  compress = FALSE
)

# ============================================================
# 16. PROVENANCE / SESSION INFO
# ============================================================

PROVENANCE <- data.table(
  component = c(
    "LDSC diagonal",
    "LDSC off-diagonal",
    "sample overlap conversion",
    "LAVA reference",
    "LAVA version"
  ),
  source = c(
    "verified STEP1D univariate LDSC intercepts",
    "verified STEP1D cross-trait gcov_int values",
    "cov2cor()",
    REF_PREFIX,
    LAVA_VERSION
  )
)

PROVENANCE_FILE <- file.path(
  TAB_DIR,
  "STEP7B2_V4_provenance.csv"
)

fwrite(
  PROVENANCE,
  PROVENANCE_FILE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "STEP7B2_V4_sessionInfo.txt"
  )
)

# ============================================================
# 17. FINAL SUMMARY
# ============================================================

cat("\n====================================================\n")
cat("STEP7B2_V4 COMPLETE — process.input() PASS\n")
cat("====================================================\n")

cat(
  "\nCommon analysis SNPs across 4 traits + UKB reference: ",
  N_ANALYSIS,
  "\n",
  sep = ""
)

cat(
  "Unalignable SNPs removed: ",
  N_UNALIGNABLE,
  "\n",
  sep = ""
)

cat("\nProcessed phenotype classification:\n")
print(
  INFO_PROCESSED[
    ,
    .(
      phenotype,
      cases,
      controls,
      prop_cases,
      binary
    )
  ]
)

cat("\nPost-harmonization QC:\n")
print(POST_QC)

cat("\nUPLOAD THESE 6 FILES:\n")
cat("1) ", PROVENANCE_FILE, "\n", sep = "")
cat("2) ", H2_FILE, "\n", sep = "")
cat("3) ", GCOV_FILE, "\n", sep = "")
cat("4) ", MATRIX_QC_FILE, "\n", sep = "")
cat("5) ", COMMON_QC_FILE, "\n", sep = "")
cat("6) ", POST_QC_FILE, "\n", sep = "")

cat("\nOnly after these pass should STEP7C local h2 be started.\n")
cat("====================================================\n")
