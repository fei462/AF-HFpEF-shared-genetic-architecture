# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP5D3_RUN_AF_OSA_CONJFDR.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP5D3 — RUN AF–OSA OFFICIAL conjFDR
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- "D:/A/data/STEP5_PLEIOFDR"

SOFTWARE_DIR <- file.path(ROOT, "software")
REF_DIR      <- file.path(ROOT, "reference")
MAT_DIR      <- file.path(ROOT, "02_mat")
RUN_DIR      <- file.path(ROOT, "07_official_runs", "AF_OSA")
RESULT_DIR   <- file.path(RUN_DIR, "results")
QC_DIR       <- file.path(ROOT, "00_qc")

dir.create(RUN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

T1 <- "AF"
T2 <- "OSA"

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

if (!length(runme_hits)) {
  stop("Cannot find pleioFDR runme.m under: ", SOFTWARE_DIR)
}

repo_candidates <- unique(dirname(runme_hits))
repo_ok <- repo_candidates[
  file.exists(file.path(repo_candidates, "pleiotropy_analysis.m"))
]

if (!length(repo_ok)) {
  stop("pleioFDR repository found, but pleiotropy_analysis.m is missing.")
}

PLEIO_REPO <- repo_ok[1]

# ---------- Required inputs ----------
REF_MAT <- file.path(
  REF_DIR,
  "ref9545380_1kgPhase3eur_LDr2p1.mat"
)
REF_TXT <- file.path(
  REF_DIR,
  "9545380.ref"
)
AF_MAT <- file.path(
  MAT_DIR,
  "AF_pleioFDR.mat"
)
OSA_MAT <- file.path(
  MAT_DIR,
  "OSA_pleioFDR.mat"
)

required <- c(
  runme = file.path(PLEIO_REPO, "runme.m"),
  pleiotropy_analysis = file.path(PLEIO_REPO, "pleiotropy_analysis.m"),
  ref_mat = REF_MAT,
  ref_txt = REF_TXT,
  AF_mat = AF_MAT,
  OSA_mat = OSA_MAT
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

cat("\n================ INPUT CHECK ================\n")
print(check, row.names = FALSE)

if (any(!check$exists)) {
  stop("Required AF/OSA pleioFDR files are missing.")
}

if (file.info(AF_MAT)$size < 100 * 1024^2) {
  stop("AF .mat is unexpectedly small.")
}
if (file.info(OSA_MAT)$size < 100 * 1024^2) {
  stop("OSA .mat is unexpectedly small.")
}

# ---------- Check previous MATLAB graphics patch ----------
PA_FILE <- file.path(PLEIO_REPO, "pleiotropy_analysis.m")
pa_txt <- readLines(PA_FILE, warn = FALSE)

fig_only_count <- sum(
  grepl("filetypes\\s*=\\s*\\{'fig'\\}", pa_txt)
)

cat("\nFIG-only graphics blocks detected:", fig_only_count, "\n")

if (fig_only_count < 1) {
  warning(
    "FIG-only graphics patch was not clearly detected. ",
    "If MATLAB again fails at PNG/SVG export, reapply the previous fix."
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

if (!length(matlab_candidates)) {
  stop("MATLAB not detected.")
}

MATLAB_EXE <- sort(matlab_candidates, decreasing = TRUE)[1]

# ---------- Write AF–OSA config ----------
normp <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

CONFIG_FILE <- file.path(PLEIO_REPO, "config.txt")

cfg <- c(
  "# AF-OSA OFFICIAL conjFDR",
  paste0("reffile=", normp(REF_MAT)),
  paste0("traitfolder=", normp(MAT_DIR)),
  "",
  "traitfile1=AF_pleioFDR.mat",
  "traitname1=AF",
  "traitfiles={'OSA_pleioFDR.mat'}",
  "traitnames={'OSA'}",
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
  paste0(
    "exclude_chr_pos=[6 ",
    MHC_FROM,
    " ",
    MHC_TO,
    "]"
  ),
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

CONFIG_ARCHIVE <- file.path(RUN_DIR, "config_AF_OSA.txt")
writeLines(cfg, CONFIG_ARCHIVE)

# ---------- Hard config check ----------
cfg_check <- trimws(readLines(CONFIG_FILE, warn = FALSE))

must_have <- c(
  "traitfile1=AF_pleioFDR.mat",
  "traitname1=AF",
  "traitfiles={'OSA_pleioFDR.mat'}",
  "traitnames={'OSA'}",
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
  cat("\n")
  stop("Config verification failed.")
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
    "AF .mat",
    "OSA .mat",
    "MATLAB",
    "AF-OSA config"
  ),
  found = c(
    file.exists(file.path(PLEIO_REPO, "runme.m")),
    file.exists(PA_FILE),
    file.exists(REF_MAT),
    file.exists(REF_TXT),
    file.exists(AF_MAT),
    file.exists(OSA_MAT),
    file.exists(MATLAB_EXE),
    file.exists(CONFIG_FILE)
  ),
  path = c(
    file.path(PLEIO_REPO, "runme.m"),
    PA_FILE,
    REF_MAT,
    REF_TXT,
    AF_MAT,
    OSA_MAT,
    MATLAB_EXE,
    CONFIG_FILE
  )
)

print(ready)

READINESS_FILE <- file.path(
  QC_DIR,
  "STEP5D3_AF_OSA_readiness.csv"
)
fwrite(ready, READINESS_FILE)

if (any(!ready$found)) {
  stop("AF–OSA readiness check failed.")
}

# ---------- Launch MATLAB ----------
repo_matlab <- gsub("'", "''", normp(PLEIO_REPO))
matlab_win <- normalizePath(
  MATLAB_EXE,
  winslash = "\\",
  mustWork = TRUE
)

matlab_code <- paste0(
  "cd('",
  repo_matlab,
  "'); run('runme.m');"
)

cmd <- paste0(
  'start "AF-OSA pleioFDR" cmd /k ""',
  matlab_win,
  '" -batch "',
  matlab_code,
  '" "'
)

cat("\n====================================================\n")
cat("LAUNCHING AF–OSA OFFICIAL conjFDR\n")
cat("====================================================\n")
cat("conjFDR threshold :", CONJFDR_THRESHOLD, "\n")
cat("randprune_n       :", RANDPRUNE_N, "\n")
cat("MAF threshold     :", MAF_THRESHOLD, "\n")
cat("MHC exclusion     : chr6:", MHC_FROM, "-", MHC_TO, "\n", sep = "")
cat("Output folder     :", RESULT_DIR, "\n")
cat("====================================================\n\n")

shell(cmd, wait = FALSE)

cat("\n====================================================\n")
cat("STEP5D3 LAUNCHED SUCCESSFULLY\n")
cat("====================================================\n")
cat("MATLAB should now run AF–OSA conjFDR in a separate window.\n")
cat("Do NOT start another pleioFDR pair at the same time.\n")
cat("MATLAB should close automatically when finished.\n\n")
cat("After completion, ZIP this folder:\n")
cat(RESULT_DIR, "\n")
cat("====================================================\n")
