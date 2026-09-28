# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP5B_OFFICIAL_pleioFDR_AF_HFpEF_V3.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP5B V3 — OFFICIAL pleioFDR/conjFDR setup + AF–HFpEF formal run
# Project: AF–HFpEF–BMI–OSA shared genetics
# V3 fix: GitHub source ZIPs are downloaded from scratch (NO Range resume).
#
# What this script does:
#   1) Download official precimed/pleiofdr source code
#   2) Download official EUR 1000G Phase 3 pleioFDR reference:
#        - ref9545380_1kgPhase3eur_LDr2p1.mat (~2.3 GB)
#        - 9545380.ref (~262 MB)
#   3) Download official precimed/python_convert source code
#   4) Create an isolated legacy-compatible conda environment for sumstats.py
#   5) Convert AF and HFpEF STEP5A inputs into official pleioFDR .mat files
#      using the official sumstats.py "mat" command
#   6) Generate an AF–HFpEF conjFDR config:
#        conjFDR < 0.05
#        random pruning n = 500
#        EUR 1000G reference
#        MAF >= 0.01
#        ambiguous SNP exclusion
#        MHC chr6:26–34 Mb (hg19), matching the reference paper
#   7) Detect MATLAB on Windows
#   8) Generate and, by default, launch a Windows .bat file for the formal run
#
# IMPORTANT:
#   - MiXeR is NOT required and may continue in the background.
#   - This script prioritizes AF–HFpEF first to minimize time/resources.
#   - It does NOT fall back to Octave, because official Octave support is
#     experimental and is not preferred for the publication-grade run.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = 1000000, error = NULL)

# ============================================================
# [0] USER SETTINGS — normally no changes are needed
# ============================================================

ROOT <- "D:/A/data"
STEP5_ROOT <- file.path(ROOT, "STEP5_PLEIOFDR")

# STEP5A prepared inputs
STEP5A_INPUT <- file.path(STEP5_ROOT, "01_inputs_prepared")

# STEP5B folders
SOFTWARE_DIR <- file.path(STEP5_ROOT, "software")
REF_DIR      <- file.path(STEP5_ROOT, "reference")
STD_DIR      <- file.path(STEP5_ROOT, "06_official_standardized")
MAT_DIR      <- file.path(STEP5_ROOT, "02_mat")
RUN_DIR      <- file.path(STEP5_ROOT, "07_official_runs", "AF_HFpEF")
QC_DIR       <- file.path(STEP5_ROOT, "00_qc")
TMP_DIR      <- file.path(STEP5_ROOT, "tmp_STEP5B")

# User's conda installation from earlier LDSC work.
CONDA_EXE_CANDIDATES <- c(
  "E:/app/python/Scripts/conda.exe",
  "E:/app/python/condabin/conda.bat",
  Sys.which("conda")
)

CONDA_ENV <- "pleiofdr310"

# Publication-grade official run
RANDPRUNE_N <- 500L
CONJFDR_THRESHOLD <- 0.05
MAF_THRESHOLD <- 0.01

# Reference-paper MHC definition: chr6:26–34 Mb, hg19
MHC_FROM <- 26000000L
MHC_TO   <- 34000000L

# TRUE = launch formal AF–HFpEF MATLAB run automatically when preparation passes.
# Set FALSE only if you want to inspect files first.
AUTO_START_MATLAB <- TRUE

# Keep preparation moderate while MiXeR is running.
DUCKDB_THREADS <- 2L
DUCKDB_MEMORY <- "5GB"

# ============================================================
# [1] FOLDERS + PACKAGES
# ============================================================

dirs <- c(STEP5_ROOT, SOFTWARE_DIR, REF_DIR, STD_DIR, MAT_DIR,
          RUN_DIR, QC_DIR, TMP_DIR)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

pkgs <- c("DBI", "duckdb", "data.table")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

library(DBI)
library(duckdb)
library(data.table)

# ============================================================
# [2] HELPER FUNCTIONS
# ============================================================

norm_win <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

sql_escape <- function(x) {
  gsub("'", "''", norm_win(x), fixed = TRUE)
}

first_existing_file <- function(x) {
  x <- x[nzchar(x)]
  hit <- x[file.exists(x)]
  if (length(hit)) hit[1] else NA_character_
}

file_size_num <- function(x) {
  if (!file.exists(x)) return(0)
  as.numeric(file.info(x)$size)
}

source_zip_ok <- function(path, expected_file) {
  if (!file.exists(path) || file_size_num(path) <= 0) return(FALSE)
  z <- tryCatch(
    utils::unzip(path, list = TRUE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(z) || nrow(z) == 0) return(FALSE)
  any(grepl(paste0("(^|/)", expected_file, "$"), z$Name, ignore.case = TRUE))
}

download_resumable <- function(urls, dest, min_bytes = 1, label = basename(dest)) {
  # V3 downloader:
  # - GitHub source ZIPs are NEVER resumed because codeload/GitHub archive
  #   endpoints may reject Range requests (curl error 33).
  # - Existing source ZIPs are validated by ZIP contents rather than file size.
  # - Large S3 reference files still use resume (-C -) when possible.
  # - Engines: curl -> PowerShell -> R/libcurl.

  urls <- unique(as.character(urls))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  is_source_zip <- grepl("\\.zip$", dest, ignore.case = TRUE) &&
    grepl("source", label, ignore.case = TRUE)

  expected_file <- if (grepl("python_convert", basename(dest), ignore.case = TRUE)) {
    "sumstats.py"
  } else if (grepl("pleiofdr", basename(dest), ignore.case = TRUE)) {
    "runme.m"
  } else {
    NA_character_
  }

  complete_now <- if (is_source_zip && !is.na(expected_file)) {
    source_zip_ok(dest, expected_file)
  } else {
    file.exists(dest) && file_size_num(dest) >= min_bytes
  }

  if (complete_now) {
    message("[OK already downloaded/validated] ", label,
            " (", round(file_size_num(dest) / 1024^2, 2), " MB)")
    return(invisible(TRUE))
  }

  # Critical fix for the user's current state:
  # remove the 0.09-MB/incomplete GitHub ZIP before retrying.
  if (is_source_zip && file.exists(dest)) {
    message("[CLEAN] Removing incomplete/non-valid source ZIP: ", dest)
    try(unlink(dest, force = TRUE), silent = TRUE)
  }

  curl <- Sys.which("curl")
  powershell <- Sys.which("powershell")

  success_check <- function() {
    if (is_source_zip && !is.na(expected_file)) {
      source_zip_ok(dest, expected_file)
    } else {
      file.exists(dest) && file_size_num(dest) >= min_bytes
    }
  }

  for (url in urls) {
    message("\n[DOWNLOAD TRY] ", label)
    message("URL: ", url)
    message("DEST: ", dest)

    is_github <- grepl(
      "(^https?://)?(www\\.)?(github\\.com|codeload\\.github\\.com|api\\.github\\.com)",
      url, ignore.case = TRUE
    )

    # --------------------------------------------------------
    # Engine 1: curl
    # --------------------------------------------------------
    if (nzchar(curl)) {
      # Never resume GitHub source archives.
      if (is_source_zip || is_github) {
        if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE)
        resume_args <- character(0)
      } else {
        resume_args <- if (file.exists(dest) && file_size_num(dest) > 0) {
          c("-C", "-")
        } else {
          character(0)
        }
      }

      args <- c(
        "--http1.1",
        "-L",
        "--fail",
        "--retry", "8",
        "--retry-all-errors",
        "--retry-delay", "5",
        "--connect-timeout", "60",
        "--speed-time", "120",
        "--speed-limit", "512",
        resume_args,
        "-o", shQuote(norm_win(dest)),
        shQuote(url)
      )

      status <- suppressWarnings(system2(curl, args = args))

      if (success_check()) {
        message("[DOWNLOAD PASS via curl] ", label,
                " (", round(file_size_num(dest) / 1024^2, 2), " MB)")
        return(invisible(TRUE))
      }

      message("curl attempt incomplete/invalid. Current size = ",
              round(file_size_num(dest) / 1024^2, 2), " MB")

      # Never carry an invalid GitHub ZIP into the next attempt.
      if (is_source_zip && file.exists(dest)) {
        try(unlink(dest, force = TRUE), silent = TRUE)
      }
    }

    # --------------------------------------------------------
    # Engine 2: Windows PowerShell
    # --------------------------------------------------------
    if (nzchar(powershell)) {
      if (is_source_zip && file.exists(dest)) {
        try(unlink(dest, force = TRUE), silent = TRUE)
      }

      ps_cmd <- paste0(
        "$ProgressPreference='SilentlyContinue'; ",
        "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; ",
        "Invoke-WebRequest -UseBasicParsing -Uri '", url,
        "' -OutFile '", gsub("'", "''", norm_win(dest), fixed = TRUE), "'"
      )

      suppressWarnings(system2(
        powershell,
        args = c("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", shQuote(ps_cmd))
      ))

      if (success_check()) {
        message("[DOWNLOAD PASS via PowerShell] ", label,
                " (", round(file_size_num(dest) / 1024^2, 2), " MB)")
        return(invisible(TRUE))
      }

      if (is_source_zip && file.exists(dest)) {
        try(unlink(dest, force = TRUE), silent = TRUE)
      }
    }

    # --------------------------------------------------------
    # Engine 3: R/libcurl (fresh download for source ZIPs)
    # --------------------------------------------------------
    if (is_source_zip) {
      if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE)

      ok_r <- tryCatch({
        utils::download.file(
          url = url,
          destfile = dest,
          method = "libcurl",
          mode = "wb",
          quiet = FALSE
        )
        TRUE
      }, error = function(e) {
        message("R/libcurl error: ", conditionMessage(e))
        FALSE
      })

      if (isTRUE(ok_r) && success_check()) {
        message("[DOWNLOAD PASS via R/libcurl] ", label,
                " (", round(file_size_num(dest) / 1024^2, 2), " MB)")
        return(invisible(TRUE))
      }

      if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE)
    }
  }

  stop(
    "All download routes failed for: ", label,
    "\nThis is a network/download problem, not a GWAS or pleioFDR analysis error.",
    if (is_source_zip) paste0(
      "\nManual fallback: download the ZIP in your browser and save it exactly as:\n",
      norm_win(dest),
      "\nThen rerun this V3 script; it will validate the ZIP and skip the download."
    ) else paste0(
      "\nCurrent size = ", round(file_size_num(dest) / 1024^2, 2), " MB",
      "\nKeep any partial S3 reference file and rerun; resume will be attempted."
    )
  )
}

run_conda <- function(args, stdout = "", stderr = "") {
  if (is.na(CONDA_EXE) || !file.exists(CONDA_EXE)) {
    stop("Conda executable not found.")
  }

  # Direct .exe is preferred. If a .bat is the only match, use cmd /c.
  if (grepl("\\.bat$", CONDA_EXE, ignore.case = TRUE)) {
    cmdargs <- c("/c", shQuote(norm_win(CONDA_EXE)), args)
    system2("cmd.exe", cmdargs, stdout = stdout, stderr = stderr)
  } else {
    system2(CONDA_EXE, args, stdout = stdout, stderr = stderr)
  }
}

detect_matlab <- function() {
  direct <- Sys.which("matlab")
  if (nzchar(direct) && file.exists(direct)) return(norm_win(direct))

  candidates <- unique(c(
    Sys.glob("C:/Program Files/MATLAB/R*/bin/matlab.exe"),
    Sys.glob("D:/Program Files/MATLAB/R*/bin/matlab.exe"),
    Sys.glob("E:/Program Files/MATLAB/R*/bin/matlab.exe"),
    Sys.glob("C:/MATLAB/R*/bin/matlab.exe"),
    Sys.glob("D:/MATLAB/R*/bin/matlab.exe"),
    Sys.glob("E:/MATLAB/R*/bin/matlab.exe")
  ))

  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)

  # Use the lexicographically latest release path.
  norm_win(sort(candidates, decreasing = TRUE)[1])
}

# ============================================================
# [3] INPUT CHECK
# ============================================================

prepared <- c(
  AF    = file.path(STEP5A_INPUT, "AF_pleioFDR_input.tsv.gz"),
  HFpEF = file.path(STEP5A_INPUT, "HFpEF_pleioFDR_input.tsv.gz")
)

if (any(!file.exists(prepared))) {
  stop(
    "STEP5A AF/HFpEF prepared input missing:\n",
    paste(names(prepared), prepared, sep = " = ", collapse = "\n")
  )
}

message("\nSTEP5A inputs found:")
print(prepared)

# ============================================================
# [4] DOWNLOAD OFFICIAL SOFTWARE
# ============================================================

PLEIO_ZIP <- file.path(SOFTWARE_DIR, "pleiofdr-master.zip")
PYCONV_ZIP <- file.path(SOFTWARE_DIR, "python_convert-master.zip")

download_resumable(
  c(
    "https://codeload.github.com/precimed/pleiofdr/zip/refs/heads/master",
    "https://github.com/precimed/pleiofdr/archive/refs/heads/master.zip",
    "https://api.github.com/repos/precimed/pleiofdr/zipball/master"
  ),
  PLEIO_ZIP,
  min_bytes = 100000,
  label = "precimed/pleiofdr source"
)

download_resumable(
  c(
    "https://codeload.github.com/precimed/python_convert/zip/refs/heads/master",
    "https://github.com/precimed/python_convert/archive/refs/heads/master.zip",
    "https://api.github.com/repos/precimed/python_convert/zipball/master"
  ),
  PYCONV_ZIP,
  min_bytes = 100000,
  label = "precimed/python_convert source"
)

PLEIO_REPO <- file.path(SOFTWARE_DIR, "pleiofdr-master")
PYCONV_REPO <- file.path(SOFTWARE_DIR, "python_convert-master")

if (!source_zip_ok(PLEIO_ZIP, "runme.m")) {
  stop("pleiofdr source ZIP is not a valid complete archive: ", PLEIO_ZIP)
}
if (!source_zip_ok(PYCONV_ZIP, "sumstats.py")) {
  stop("python_convert source ZIP is not a valid complete archive: ", PYCONV_ZIP)
}

if (!file.exists(file.path(PLEIO_REPO, "runme.m"))) {
  message("[EXTRACT] pleiofdr")
  unzip(PLEIO_ZIP, exdir = SOFTWARE_DIR)
}
if (!file.exists(file.path(PYCONV_REPO, "sumstats.py"))) {
  message("[EXTRACT] python_convert")
  unzip(PYCONV_ZIP, exdir = SOFTWARE_DIR)
}

# GitHub API zipball fallback may use a commit-hash folder name instead of
# pleiofdr-master/python_convert-master. Detect the extracted folder safely.
if (!file.exists(file.path(PLEIO_REPO, "runme.m"))) {
  cand <- list.dirs(SOFTWARE_DIR, recursive = FALSE, full.names = TRUE)
  hit <- cand[file.exists(file.path(cand, "runme.m"))]
  if (length(hit)) PLEIO_REPO <- hit[1]
}
if (!file.exists(file.path(PYCONV_REPO, "sumstats.py"))) {
  cand <- list.dirs(SOFTWARE_DIR, recursive = FALSE, full.names = TRUE)
  hit <- cand[file.exists(file.path(cand, "sumstats.py"))]
  if (length(hit)) PYCONV_REPO <- hit[1]
}

if (!file.exists(file.path(PLEIO_REPO, "runme.m"))) {
  stop("pleiofdr extraction failed: runme.m not found")
}
if (!file.exists(file.path(PYCONV_REPO, "sumstats.py"))) {
  stop("python_convert extraction failed: sumstats.py not found")
}

# ============================================================
# [5] DOWNLOAD OFFICIAL 1000G EUR pleioFDR REFERENCES
# ============================================================

REF_MAT <- file.path(REF_DIR, "ref9545380_1kgPhase3eur_LDr2p1.mat")
REF_TXT <- file.path(REF_DIR, "9545380.ref")
ABOUT_TXT <- file.path(REF_DIR, "about.txt")

download_resumable(
  "https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt",
  ABOUT_TXT,
  min_bytes = 1000,
  label = "pleioFDR about.txt"
)

# ~262 MB
download_resumable(
  "https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/9545380.ref",
  REF_TXT,
  min_bytes = 200000000,
  label = "9545380.ref"
)

# ~2.3 GB
download_resumable(
  "https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat",
  REF_MAT,
  min_bytes = 2000000000,
  label = "ref9545380_1kgPhase3eur_LDr2p1.mat"
)

# ============================================================
# [6] REFERENCE QC
# ============================================================

ref_head <- fread(REF_TXT, nrows = 5, showProgress = FALSE)
fwrite(ref_head, file.path(QC_DIR, "STEP5B_reference_head5.tsv"),
       sep = "\t", quote = FALSE)

required_ref_cols <- c("CHR", "SNP", "BP", "A1", "A2")
if (!all(required_ref_cols %in% names(ref_head))) {
  stop(
    "9545380.ref does not contain expected columns: ",
    paste(required_ref_cols, collapse = ", "),
    "\nObserved: ", paste(names(ref_head), collapse = ", ")
  )
}

# ============================================================
# [7] BUILD OFFICIAL-STANDARDIZED TSV.GZ FOR AF/HFpEF
#
# We assign CHR/BP FROM THE OFFICIAL 9545380.ref by rsID.
# GWAS A1/A2 are retained, so official sumstats.py can perform
# allele matching and Z-score orientation against the reference.
# ============================================================

std_files <- c(
  AF    = file.path(STD_DIR, "AF_pleioFDR_standard.tsv.gz"),
  HFpEF = file.path(STD_DIR, "HFpEF_pleioFDR_standard.tsv.gz")
)

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = file.path(TMP_DIR, "STEP5B_prepare.duckdb"),
  read_only = FALSE
)

on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(con, sprintf("SET threads=%d;", DUCKDB_THREADS))
dbExecute(con, sprintf("SET memory_limit='%s';", DUCKDB_MEMORY))
dbExecute(
  con,
  sprintf("SET temp_directory='%s';", sql_escape(TMP_DIR))
)

ref_sql <- sql_escape(REF_TXT)

prep_qc <- list()

for (tr in names(prepared)) {
  infile <- sql_escape(prepared[[tr]])
  outfile <- sql_escape(std_files[[tr]])

  message("\n[STANDARDIZE FOR OFFICIAL sumstats.py] ", tr)

  # Current prepared input has:
  # SNP A1 A2 P Z N LOGP
  #
  # Official sumstats.py "mat" expects a standardized table including:
  # SNP CHR BP A1 A2 PVAL Z N
  #
  # The official converter then aligns by SNP and allele to 9545380.ref.
  query <- sprintf("
    COPY (
      SELECT
        CAST(g.SNP AS VARCHAR) AS SNP,
        CAST(r.CHR AS INTEGER) AS CHR,
        CAST(r.BP AS BIGINT) AS BP,
        upper(CAST(g.A1 AS VARCHAR)) AS A1,
        upper(CAST(g.A2 AS VARCHAR)) AS A2,
        CAST(g.P AS DOUBLE) AS PVAL,
        CAST(g.Z AS DOUBLE) AS Z,
        CAST(g.N AS DOUBLE) AS N
      FROM read_csv_auto(
             '%s',
             delim='\\t',
             header=true,
             compression='gzip',
             sample_size=100000
           ) AS g
      INNER JOIN read_csv_auto(
             '%s',
             delim='\\t',
             header=true,
             sample_size=100000
           ) AS r
      ON CAST(g.SNP AS VARCHAR) = CAST(r.SNP AS VARCHAR)
      WHERE
        g.SNP IS NOT NULL
        AND g.P IS NOT NULL
        AND g.Z IS NOT NULL
        AND g.N IS NOT NULL
        AND CAST(g.P AS DOUBLE) > 0
        AND CAST(g.P AS DOUBLE) <= 1
        AND isfinite(CAST(g.Z AS DOUBLE))
        AND isfinite(CAST(g.N AS DOUBLE))
    )
    TO '%s'
    (FORMAT CSV, DELIMITER '\\t', HEADER TRUE, COMPRESSION GZIP);
  ", infile, ref_sql, outfile)

  if (!file.exists(std_files[[tr]])) {
    dbExecute(con, query)
  } else {
    message("[OK already standardized] ", std_files[[tr]])
  }

  # Count rows and check columns.
  qcount <- sprintf("
    SELECT
      count(*) AS n,
      count(DISTINCT SNP) AS n_unique_snp,
      min(PVAL) AS min_p,
      max(abs(Z)) AS max_abs_z,
      median(N) AS median_n
    FROM read_csv_auto(
      '%s', delim='\\t', header=true, compression='gzip', sample_size=100000
    );
  ", outfile)

  qc <- as.data.table(dbGetQuery(con, qcount))
  qc[, trait := tr]
  qc[, file := std_files[[tr]]]
  setcolorder(qc, c("trait", "n", "n_unique_snp", "min_p",
                    "max_abs_z", "median_n", "file"))
  prep_qc[[tr]] <- qc
}

prep_qc_dt <- rbindlist(prep_qc)
fwrite(prep_qc_dt,
       file.path(QC_DIR, "STEP5B_official_standardized_QC.csv"))

dbDisconnect(con, shutdown = TRUE)
con <- NULL

# ============================================================
# [8] FIND / CREATE ISOLATED CONDA ENV FOR OFFICIAL python_convert
#
# The historical converter uses APIs removed in modern pandas/numpy.
# Use an isolated compatibility environment; do NOT modify ldsc39.
# ============================================================

CONDA_EXE <- first_existing_file(CONDA_EXE_CANDIDATES)

if (is.na(CONDA_EXE)) {
  stop(
    "Conda executable not found.\n",
    "Expected e.g. E:/app/python/Scripts/conda.exe\n",
    "Do not modify your ldsc39 environment."
  )
}

message("\nConda: ", CONDA_EXE)

env_list <- tryCatch(
  run_conda(c("env", "list"), stdout = TRUE, stderr = TRUE),
  error = function(e) character(0)
)

env_exists <- any(grepl(paste0("(^|[[:space:]])", CONDA_ENV, "([[:space:]]|$)"),
                        env_list))

if (!env_exists) {
  message("\n[CONDA] Creating isolated env: ", CONDA_ENV)
  status <- run_conda(c(
    "create", "-y", "-n", CONDA_ENV,
    "python=3.10",
    "numpy=1.23.5",
    "scipy=1.10.1",
    "pandas=1.5.3",
    "six=1.16.0"
  ))

  if (!identical(as.integer(status), 0L)) {
    stop("Failed to create conda environment ", CONDA_ENV)
  }
} else {
  message("[OK conda env exists] ", CONDA_ENV)
}

# Test versions
pyver <- run_conda(
  c("run", "-n", CONDA_ENV, "python", "-c",
    shQuote("import sys,numpy,pandas,scipy,six; print(sys.version); print('numpy',numpy.__version__); print('pandas',pandas.__version__); print('scipy',scipy.__version__)")),
  stdout = TRUE, stderr = TRUE
)
writeLines(pyver, file.path(QC_DIR, "STEP5B_python_environment.txt"))
cat(paste(pyver, collapse = "\n"), "\n")

# ============================================================
# [9] OFFICIAL sumstats.py -> MATLAB .mat
# ============================================================

SUMSTATS_PY <- file.path(PYCONV_REPO, "sumstats.py")

mat_files <- c(
  AF    = file.path(MAT_DIR, "AF_pleioFDR.mat"),
  HFpEF = file.path(MAT_DIR, "HFpEF_pleioFDR.mat")
)

conversion_log <- file.path(QC_DIR, "STEP5B_sumstats_mat_conversion.log")
if (file.exists(conversion_log)) file.remove(conversion_log)

for (tr in names(std_files)) {
  message("\n[OFFICIAL MAT CONVERSION] ", tr)

  if (file.exists(mat_files[[tr]]) && file_size_num(mat_files[[tr]]) > 100000000) {
    message("[OK existing .mat] ", mat_files[[tr]],
            " (", round(file_size_num(mat_files[[tr]]) / 1024^2, 1), " MB)")
    next
  }

  args <- c(
    "run", "-n", CONDA_ENV,
    "python", norm_win(SUMSTATS_PY),
    "mat",
    "--sumstats", norm_win(std_files[[tr]]),
    "--ref", norm_win(REF_TXT),
    "--out", norm_win(mat_files[[tr]]),
    "--chunksize", "500000",
    "--force"
  )

  out <- run_conda(args, stdout = TRUE, stderr = TRUE)
  cat(
    paste0("\n========== ", tr, " ==========\n"),
    paste(out, collapse = "\n"),
    "\n",
    file = conversion_log,
    append = TRUE
  )

  if (!file.exists(mat_files[[tr]]) || file_size_num(mat_files[[tr]]) < 100000000) {
    stop(
      tr, " .mat conversion did not produce a plausible output file.\n",
      "Inspect: ", conversion_log
    )
  }
}

# ============================================================
# [10] MAT QC USING SCIPY
# ============================================================

CHECK_PY <- file.path(TMP_DIR, "check_pleiofdr_mat.py")

writeLines(c(
  "import sys, os, csv",
  "import numpy as np",
  "import scipy.io as sio",
  "",
  "outfile = sys.argv[1]",
  "pairs = []",
  "for x in sys.argv[2:]:",
  "    trait, path = x.split('=', 1)",
  "    pairs.append((trait, path))",
  "",
  "rows = []",
  "for trait, path in pairs:",
  "    d = sio.loadmat(path, variable_names=['logpvec','zvec','nvec'])",
  "    lp = np.ravel(d['logpvec'])",
  "    z = np.ravel(d['zvec'])",
  "    n = np.ravel(d['nvec'])",
  "    rows.append({",
  "      'trait': trait,",
  "      'n_reference_slots': len(lp),",
  "      'finite_logp': int(np.isfinite(lp).sum()),",
  "      'finite_z': int(np.isfinite(z).sum()),",
  "      'finite_n': int(np.isfinite(n).sum()),",
  "      'min_logp': float(np.nanmin(lp)),",
  "      'max_logp': float(np.nanmax(lp)),",
  "      'max_abs_z': float(np.nanmax(np.abs(z))),",
  "      'median_n': float(np.nanmedian(n)),",
  "      'mat_size_mb': os.path.getsize(path)/(1024**2),",
  "      'path': path",
  "    })",
  "",
  "with open(outfile, 'w', newline='') as f:",
  "    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))",
  "    w.writeheader()",
  "    w.writerows(rows)",
  "print('MAT QC written:', outfile)"
), CHECK_PY)

MAT_QC_CSV <- file.path(QC_DIR, "STEP5B_MAT_QC.csv")

qc_args <- c(
  "run", "-n", CONDA_ENV,
  "python", norm_win(CHECK_PY),
  norm_win(MAT_QC_CSV),
  paste0("AF=", norm_win(mat_files["AF"])),
  paste0("HFpEF=", norm_win(mat_files["HFpEF"]))
)

qc_out <- run_conda(qc_args, stdout = TRUE, stderr = TRUE)
cat(paste(qc_out, collapse = "\n"), "\n")

if (!file.exists(MAT_QC_CSV)) {
  stop("MAT QC failed.")
}

mat_qc <- fread(MAT_QC_CSV)
print(mat_qc)

# Hard gate: official reference has 9,545,380 slots.
EXPECTED_REF_N <- 9545380L

if (any(mat_qc$n_reference_slots != EXPECTED_REF_N)) {
  stop(
    "MAT vector length does not match official reference length 9,545,380.\n",
    "Do NOT run pleioFDR until fixed."
  )
}

if (any(mat_qc$finite_z < 100000)) {
  stop(
    "Too few finite Z scores after official reference/allele matching.\n",
    "Inspect allele/reference compatibility before running pleioFDR."
  )
}

# ============================================================
# [11] MATLAB DETECTION
# ============================================================

MATLAB_EXE <- detect_matlab()

matlab_status <- data.table(
  item = c(
    "pleioFDR repository",
    "official LD reference .mat",
    "official reference info",
    "AF official .mat",
    "HFpEF official .mat",
    "MATLAB"
  ),
  found = c(
    file.exists(file.path(PLEIO_REPO, "runme.m")),
    file.exists(REF_MAT) && file_size_num(REF_MAT) >= 2000000000,
    file.exists(REF_TXT) && file_size_num(REF_TXT) >= 200000000,
    file.exists(mat_files["AF"]),
    file.exists(mat_files["HFpEF"]),
    !is.na(MATLAB_EXE)
  ),
  path = c(
    PLEIO_REPO,
    REF_MAT,
    REF_TXT,
    mat_files["AF"],
    mat_files["HFpEF"],
    MATLAB_EXE
  )
)

fwrite(
  matlab_status,
  file.path(QC_DIR, "STEP5B_official_readiness.csv")
)

print(matlab_status)

# ============================================================
# [12] AF–HFpEF OFFICIAL conjFDR CONFIG
# ============================================================

RESULT_DIR <- file.path(RUN_DIR, "results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

CONFIG_FILE <- file.path(PLEIO_REPO, "config.txt")

cfg <- c(
  "# STEP5B AF-HFpEF publication-grade conjFDR",
  "# Generated automatically; based on official precimed/pleiofdr config.",
  paste0("reffile=", norm_win(REF_MAT)),
  paste0("traitfolder=", norm_win(MAT_DIR)),
  "",
  "traitfile1=AF_pleioFDR.mat",
  "traitname1=AF",
  "traitfiles={'HFpEF_pleioFDR.mat'}",
  "traitnames={'HFpEF'}",
  "",
  paste0("outputdir=", norm_win(RESULT_DIR)),
  "stattype=conjfdr",
  paste0("fdrthresh=", CONJFDR_THRESHOLD),
  "",
  "randprune=true",
  paste0("randprune_n=", RANDPRUNE_N),
  paste0("exclude_chr_pos=[6 ", MHC_FROM, " ", MHC_TO, "]"),
  "",
  "manh_fontsize_genenames=12",
  "manh_yspace=0.75",
  "manh_ymargin=0.25",
  "manh_colorlist=[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]",
  paste0("refinfo=", norm_win(REF_TXT)),
  "",
  "reset_pruneidx=true",
  "randprune_repeats=default",
  "pthresh=1",
  "perform_gc=true",
  "use_standard_gc=false",
  "randprune_gc=true",
  "exclude_from_discovery=false",
  paste0("mafthresh=", MAF_THRESHOLD),
  "exclude_ambiguous_snps=true",
  "onscreen=false",
  "dummy_zscore=false",
  "exit_matlab_upon_completion=true"
)

writeLines(cfg, CONFIG_FILE)
file.copy(CONFIG_FILE,
          file.path(RUN_DIR, "config_AF_HFpEF_conjFDR.txt"),
          overwrite = TRUE)

# ============================================================
# [13] WINDOWS LAUNCHER
# ============================================================

BAT_FILE <- file.path(RUN_DIR, "RUN_AF_HFpEF_CONJFDR.bat")
MATLAB_LOG <- file.path(RUN_DIR, "AF_HFpEF_MATLAB_console.log")
STATUS_TXT <- file.path(RUN_DIR, "AF_HFpEF_exit_status.txt")

if (!is.na(MATLAB_EXE)) {
  bat <- c(
    "@echo off",
    "setlocal",
    paste0("cd /d \"", gsub("/", "\\\\", norm_win(PLEIO_REPO)), "\""),
    paste0("\"", gsub("/", "\\\\", MATLAB_EXE),
           "\" -nodisplay -nosplash < runme.m > \"",
           gsub("/", "\\\\", norm_win(MATLAB_LOG)), "\" 2>&1"),
    paste0("echo %ERRORLEVEL% > \"",
           gsub("/", "\\\\", norm_win(STATUS_TXT)), "\""),
    "endlocal"
  )
  writeLines(bat, BAT_FILE)

  message("\nMATLAB detected: ", MATLAB_EXE)
  message("Launcher created: ", BAT_FILE)

  if (AUTO_START_MATLAB) {
    message("\n============================================================")
    message("STARTING OFFICIAL AF-HFpEF conjFDR IN BACKGROUND")
    message("randprune_n = ", RANDPRUNE_N)
    message("conjFDR threshold = ", CONJFDR_THRESHOLD)
    message("============================================================")

    # Background launch so RStudio returns promptly.
    shell(
      paste0('start "" /B "', gsub("/", "\\\\", norm_win(BAT_FILE)), '"'),
      wait = FALSE
    )

    writeLines(
      c(
        paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("MATLAB: ", MATLAB_EXE),
        paste0("Config: ", CONFIG_FILE),
        paste0("Result dir: ", RESULT_DIR),
        paste0("Console log: ", MATLAB_LOG)
      ),
      file.path(RUN_DIR, "STEP5B_LAUNCHED.txt")
    )
  }
} else {
  warning(
    "\nMATLAB was NOT detected.\n",
    "All official pleioFDR resources and AF/HFpEF .mat inputs are prepared, ",
    "but the formal conjFDR run has NOT been launched.\n",
    "Install/locate MATLAB and rerun from section [11] onward.\n",
    "Do NOT substitute Octave for the publication-grade primary analysis."
  )
}

# ============================================================
# [14] FINAL STATUS / FILES TO RETURN
# ============================================================

final_note <- c(
  "STEP5B official pleioFDR preparation completed.",
  "",
  "Primary pair: AF-HFpEF",
  paste0("conjFDR threshold: ", CONJFDR_THRESHOLD),
  paste0("randprune_n: ", RANDPRUNE_N),
  paste0("MAF threshold: ", MAF_THRESHOLD),
  paste0("MHC excluded from model fitting: chr6:", MHC_FROM, "-", MHC_TO, " (hg19)"),
  "",
  "Official input conversion:",
  "- SNP/reference alignment performed against 9545380.ref",
  "- GWAS A1/A2 retained for official allele matching",
  "- Z orientation handled by official precimed/python_convert sumstats.py mat",
  "",
  "Files to send back BEFORE interpreting formal conjFDR:",
  file.path(QC_DIR, "STEP5B_official_standardized_QC.csv"),
  file.path(QC_DIR, "STEP5B_MAT_QC.csv"),
  file.path(QC_DIR, "STEP5B_official_readiness.csv"),
  file.path(QC_DIR, "STEP5B_sumstats_mat_conversion.log"),
  "",
  "When MATLAB run finishes, also send:",
  MATLAB_LOG,
  RESULT_DIR,
  "",
  "Do not alter conjFDR<0.05 if shared loci are few."
)

writeLines(final_note, file.path(RUN_DIR, "README_STEP5B.txt"))

message("\n============================================================")
message("STEP5B PREPARATION COMPLETE")
message("QC: ", file.path(QC_DIR, "STEP5B_MAT_QC.csv"))
message("Readiness: ", file.path(QC_DIR, "STEP5B_official_readiness.csv"))
if (!is.na(MATLAB_EXE) && AUTO_START_MATLAB) {
  message("Formal AF-HFpEF conjFDR has been launched in background.")
  message("Monitor log: ", MATLAB_LOG)
} else {
  message("Formal MATLAB run has not been started.")
}
message("============================================================")
