# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP5D4_RUN_HFpEF_BMI_CONJFDR.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP5D4 — RUN HFpEF–BMI OFFICIAL conjFDR
# Reuse existing HFpEF/BMI .mat files and official reference
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- "D:/A/data/STEP5_PLEIOFDR"

SOFTWARE_DIR <- file.path(ROOT, "software")
REF_DIR      <- file.path(ROOT, "reference")
MAT_DIR      <- file.path(ROOT, "02_mat")
RUN_DIR      <- file.path(ROOT, "07_official_runs", "HFpEF_BMI")
RESULT_DIR   <- file.path(RUN_DIR, "results")
QC_DIR       <- file.path(ROOT, "00_qc")

dir.create(RUN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

CONJFDR_THRESHOLD <- 0.05
RANDPRUNE_N <- 500L
MAF_THRESHOLD <- 0.01
MHC_FROM <- 26000000L
MHC_TO   <- 34000000L

# ---------- Find pleioFDR repository ----------
runme_hits <- list.files(
  SOFTWARE_DIR,
  pattern = "^runme\\.m$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(runme_hits)) stop("Cannot find pleioFDR runme.m.")

repo_candidates <- unique(dirname(runme_hits))
repo_ok <- repo_candidates[
  file.exists(file.path(repo_candidates, "pleiotropy_analysis.m"))
]
if (!length(repo_ok)) stop("Cannot find pleiotropy_analysis.m.")
PLEIO_REPO <- repo_ok[1]

# ---------- Required inputs ----------
REF_MAT <- file.path(REF_DIR, "ref9545380_1kgPhase3eur_LDr2p1.mat")
REF_TXT <- file.path(REF_DIR, "9545380.ref")
HF_MAT  <- file.path(MAT_DIR, "HFpEF_pleioFDR.mat")
BMI_MAT <- file.path(MAT_DIR, "BMI_pleioFDR.mat")

required <- c(
  runme = file.path(PLEIO_REPO, "runme.m"),
  pleiotropy_analysis = file.path(PLEIO_REPO, "pleiotropy_analysis.m"),
  ref_mat = REF_MAT,
  ref_txt = REF_TXT,
  HFpEF_mat = HF_MAT,
  BMI_mat = BMI_MAT
)

check <- data.frame(
  item = names(required),
  exists = file.exists(required),
  size_MB = ifelse(
    file.exists(required),
    round(file.info(required)$size / 1024^2, 2),
    NA_real_
  ),
  path = unname(required),
  stringsAsFactors = FALSE
)

print(check, row.names = FALSE)
if (any(!check$exists)) stop("Required files are missing.")

if (file.info(HF_MAT)$size < 100 * 1024^2) stop("HFpEF .mat unexpectedly small.")
if (file.info(BMI_MAT)$size < 100 * 1024^2) stop("BMI .mat unexpectedly small.")

# ---------- Check FIG-only graphics patch ----------
PA_FILE <- file.path(PLEIO_REPO, "pleiotropy_analysis.m")
pa_txt <- readLines(PA_FILE, warn = FALSE)
fig_only_count <- sum(grepl("filetypes\\s*=\\s*\\{'fig'\\}", pa_txt))
cat("\nFIG-only graphics blocks detected:", fig_only_count, "\n")

if (fig_only_count < 1) {
  warning(
    "FIG-only patch not clearly detected. ",
    "If MATLAB fails during PNG/SVG export, reapply the previous patch."
  )
}

# ---------- Find MATLAB ----------
matlab_candidates <- unique(c(
  unname(Sys.which("matlab")),
  Sys.glob("C:/Program Files/MATLAB/R*/bin/matlab.exe"),
  Sys.glob("D:/Program Files/MATLAB/R*/bin/matlab.exe"),
  Sys.glob("E:/Program Files/MATLAB/R*/bin/matlab.exe"),
  Sys.glob("C:/MATLAB/R*/bin/matlab.exe"),
  Sys.glob("D:/MATLAB/R*/bin/matlab.exe"),
  Sys.glob("E:/MATLAB/R*/bin/matlab.exe")
))
matlab_candidates <- matlab_candidates[
  nzchar(matlab_candidates) & file.exists(matlab_candidates)
]
if (!length(matlab_candidates)) stop("MATLAB not detected.")
MATLAB_EXE <- sort(matlab_candidates, decreasing = TRUE)[1]

# ---------- Write config ----------
normp <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

CONFIG_FILE <- file.path(PLEIO_REPO, "config.txt")

cfg <- c(
  "# HFpEF-BMI OFFICIAL conjFDR",
  paste0("reffile=", normp(REF_MAT)),
  paste0("traitfolder=", normp(MAT_DIR)),
  "",
  "traitfile1=HFpEF_pleioFDR.mat",
  "traitname1=HFpEF",
  "traitfiles={'BMI_pleioFDR.mat'}",
  "traitnames={'BMI'}",
  "",
  paste0("outputdir=", normp(RESULT_DIR)),
  "stattype=conjfdr",
  paste0("fdrthresh=", CONJFDR_THRESHOLD),
  "",
  "randprune=true",
  paste0("randprune_n=", RANDPRUNE_N),
  "randprune_repeats=default",
  "reset_pruneidx=true",
  "",
  paste0("exclude_chr_pos=[6 ", MHC_FROM, " ", MHC_TO, "]"),
  "exclude_from_discovery=false",
  "",
  paste0("mafthresh=", MAF_THRESHOLD),
  "exclude_ambiguous_snps=true",
  "",
  "pthresh=1",
  "perform_gc=true",
  "use_standard_gc=false",
  "randprune_gc=true",
  "",
  paste0("refinfo=", normp(REF_TXT)),
  "",
  "onscreen=false",
  "dummy_zscore=false",
  "manh_plot=1",
  "exit_matlab_upon_completion=true"
)

writeLines(cfg, CONFIG_FILE)
writeLines(cfg, file.path(RUN_DIR, "config_HFpEF_BMI.txt"))

# ---------- Verify config ----------
cfg_check <- trimws(readLines(CONFIG_FILE, warn = FALSE))
must_have <- c(
  "traitfile1=HFpEF_pleioFDR.mat",
  "traitname1=HFpEF",
  "traitfiles={'BMI_pleioFDR.mat'}",
  "traitnames={'BMI'}",
  "stattype=conjfdr",
  "fdrthresh=0.05",
  "randprune=true",
  "randprune_n=500",
  "mafthresh=0.01",
  "exclude_ambiguous_snps=true",
  "exit_matlab_upon_completion=true"
)

missing_cfg <- must_have[!must_have %in% cfg_check]
if (length(missing_cfg)) {
  cat("\nMissing config lines:\n")
  cat(paste0("  ", missing_cfg), sep = "\n")
  stop("\nConfig verification failed.")
}

# ---------- Readiness ----------
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
library(data.table)

ready <- data.table(
  item = c(
    "pleioFDR runme",
    "pleiotropy_analysis",
    "official LD reference",
    "official refinfo",
    "HFpEF .mat",
    "BMI .mat",
    "MATLAB",
    "HFpEF-BMI config"
  ),
  found = c(
    file.exists(file.path(PLEIO_REPO, "runme.m")),
    file.exists(PA_FILE),
    file.exists(REF_MAT),
    file.exists(REF_TXT),
    file.exists(HF_MAT),
    file.exists(BMI_MAT),
    file.exists(MATLAB_EXE),
    file.exists(CONFIG_FILE)
  ),
  path = c(
    file.path(PLEIO_REPO, "runme.m"),
    PA_FILE,
    REF_MAT,
    REF_TXT,
    HF_MAT,
    BMI_MAT,
    MATLAB_EXE,
    CONFIG_FILE
  )
)

print(ready)
READINESS_FILE <- file.path(QC_DIR, "STEP5D4_HFpEF_BMI_readiness.csv")
fwrite(ready, READINESS_FILE)

if (any(!ready$found)) stop("HFpEF-BMI readiness check failed.")

# ---------- Launch MATLAB ----------
repo_matlab <- gsub("'", "''", normp(PLEIO_REPO))
matlab_win <- normalizePath(MATLAB_EXE, winslash = "\\", mustWork = TRUE)

matlab_code <- paste0(
  "cd('", repo_matlab, "'); run('runme.m');"
)

cmd <- paste0(
  'start "HFpEF-BMI pleioFDR" cmd /k ""',
  matlab_win,
  '" -batch "',
  matlab_code,
  '" "'
)

cat("\n====================================================\n")
cat("LAUNCHING HFpEF–BMI OFFICIAL conjFDR\n")
cat("====================================================\n")
cat("conjFDR threshold :", CONJFDR_THRESHOLD, "\n")
cat("randprune_n       :", RANDPRUNE_N, "\n")
cat("MAF threshold     :", MAF_THRESHOLD, "\n")
cat("Output folder     :", RESULT_DIR, "\n")
cat("====================================================\n")

shell(cmd, wait = FALSE)

cat("\nSTEP5D4 LAUNCHED SUCCESSFULLY\n")
cat("Do NOT start another pleioFDR pair simultaneously.\n")
cat("After completion, ZIP this folder:\n")
cat(RESULT_DIR, "\n")
