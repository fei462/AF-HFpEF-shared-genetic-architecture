# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP5D2_PREP_BMI_OSA_RUN_NEXT_CONJFDR.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP5D2 — Prepare BMI/OSA official pleioFDR .mat + run next pair
# Default next formal pair: AF–BMI
#
# After AF–BMI finishes, change ONE line:
#   PAIR_TO_RUN <- "AF_OSA"
# and run this script again. Existing BMI/OSA .mat files will be reused.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999, timeout = 1000000)

# ================= USER SETTINGS =================
ROOT <- "D:/A/data/STEP5_PLEIOFDR"

PAIR_TO_RUN <- "AF_BMI"    # <-- NEXT TIME change to "AF_OSA"
AUTO_START_MATLAB <- TRUE

RANDPRUNE_N <- 500L
CONJFDR_THRESHOLD <- 0.05
MAF_THRESHOLD <- 0.01
MHC_FROM <- 26000000L
MHC_TO <- 34000000L

# Keep moderate while MiXeR may still be running.
DUCKDB_THREADS <- 2L
DUCKDB_MEMORY <- "5GB"

# ================= PATHS =================
STEP5A_INPUT <- file.path(ROOT, "01_inputs_prepared")
SOFTWARE_DIR <- file.path(ROOT, "software")
REF_DIR <- file.path(ROOT, "reference")
STD_DIR <- file.path(ROOT, "06_official_standardized")
MAT_DIR <- file.path(ROOT, "02_mat")
QC_DIR <- file.path(ROOT, "00_qc")
TMP_DIR <- file.path(ROOT, "tmp_STEP5D2")
RUN_BASE <- file.path(ROOT, "07_official_runs")

invisible(lapply(
  c(STD_DIR, MAT_DIR, QC_DIR, TMP_DIR, RUN_BASE),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# ================= PACKAGES =================
pkgs <- c("DBI","duckdb","data.table")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(to_install)) {
  install.packages(to_install, repos="https://cloud.r-project.org")
}

library(DBI)
library(duckdb)
library(data.table)

normp <- function(x) normalizePath(x, winslash="/", mustWork=FALSE)
sq <- function(x) gsub("'", "''", normp(x), fixed=TRUE)

# ================= FIND SOFTWARE =================
# Find pleioFDR repo from runme.m
runme_hits <- list.files(
  SOFTWARE_DIR,
  pattern="^runme\\.m$",
  recursive=TRUE,
  full.names=TRUE
)
if (!length(runme_hits)) stop("Cannot find pleioFDR runme.m under ", SOFTWARE_DIR)
PLEIO_REPO <- dirname(runme_hits[1])

# Find official converter
sumstats_hits <- list.files(
  SOFTWARE_DIR,
  pattern="^sumstats\\.py$",
  recursive=TRUE,
  full.names=TRUE
)
if (!length(sumstats_hits)) stop("Cannot find python_convert/sumstats.py under ", SOFTWARE_DIR)
SUMSTATS_PY <- sumstats_hits[1]

REF_TXT <- file.path(REF_DIR, "9545380.ref")
REF_MAT <- file.path(REF_DIR, "ref9545380_1kgPhase3eur_LDr2p1.mat")

if (!file.exists(REF_TXT)) stop("Missing ", REF_TXT)
if (!file.exists(REF_MAT)) stop("Missing ", REF_MAT)

# ================= FIND CONDA =================
conda_candidates <- c(
  "E:/app/python/Scripts/conda.exe",
  "E:/app/python/condabin/conda.bat",
  unname(Sys.which("conda"))
)
conda_candidates <- conda_candidates[nzchar(conda_candidates) & file.exists(conda_candidates)]
if (!length(conda_candidates)) stop("Conda not found.")
CONDA_EXE <- conda_candidates[1]
CONDA_ENV <- "pleiofdr310"

run_conda <- function(args, stdout=TRUE, stderr=TRUE) {
  if (grepl("\\.bat$", CONDA_EXE, ignore.case=TRUE)) {
    system2("cmd.exe", c("/c", shQuote(normp(CONDA_EXE)), args),
            stdout=stdout, stderr=stderr)
  } else {
    system2(CONDA_EXE, args, stdout=stdout, stderr=stderr)
  }
}

envs <- run_conda(c("env","list"), stdout=TRUE, stderr=TRUE)
if (!any(grepl(paste0("\\b", CONDA_ENV, "\\b"), envs))) {
  stop("Conda environment '", CONDA_ENV,
       "' not found. It should already exist from STEP5B.")
}

# ================= PREPARE BMI + OSA =================
prepared <- c(
  BMI = file.path(STEP5A_INPUT, "BMI_pleioFDR_input.tsv.gz"),
  OSA = file.path(STEP5A_INPUT, "OSA_pleioFDR_input.tsv.gz")
)
if (any(!file.exists(prepared))) {
  stop("Missing STEP5A BMI/OSA input:\n",
       paste(names(prepared), prepared, sep=" = ", collapse="\n"))
}

std_files <- c(
  BMI = file.path(STD_DIR, "BMI_pleioFDR_standard.tsv.gz"),
  OSA = file.path(STD_DIR, "OSA_pleioFDR_standard.tsv.gz")
)

mat_files <- c(
  AF = file.path(MAT_DIR, "AF_pleioFDR.mat"),
  HFpEF = file.path(MAT_DIR, "HFpEF_pleioFDR.mat"),
  BMI = file.path(MAT_DIR, "BMI_pleioFDR.mat"),
  OSA = file.path(MAT_DIR, "OSA_pleioFDR.mat")
)

if (!file.exists(mat_files["AF"])) stop("AF .mat missing.")
if (!file.exists(mat_files["HFpEF"])) stop("HFpEF .mat missing.")

# ---- standardize by official 9545380.ref ----
con <- dbConnect(
  duckdb::duckdb(),
  dbdir=file.path(TMP_DIR, "prepare_BMI_OSA.duckdb"),
  read_only=FALSE
)
dbExecute(con, sprintf("SET threads=%d;", DUCKDB_THREADS))
dbExecute(con, sprintf("SET memory_limit='%s';", DUCKDB_MEMORY))
dbExecute(con, sprintf("SET temp_directory='%s';", sq(TMP_DIR)))

qc_std <- list()

for (tr in names(prepared)) {
  message("\n[STANDARDIZE] ", tr)

  if (!file.exists(std_files[[tr]])) {
    q <- sprintf("
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
          '%s', delim='\\t', header=true, compression='gzip',
          sample_size=100000
        ) g
        INNER JOIN read_csv_auto(
          '%s', delim='\\t', header=true, sample_size=100000
        ) r
        ON CAST(g.SNP AS VARCHAR) = CAST(r.SNP AS VARCHAR)
        WHERE
          g.SNP IS NOT NULL
          AND CAST(g.P AS DOUBLE) > 0
          AND CAST(g.P AS DOUBLE) <= 1
          AND isfinite(CAST(g.Z AS DOUBLE))
          AND isfinite(CAST(g.N AS DOUBLE))
      )
      TO '%s'
      (FORMAT CSV, DELIMITER '\\t', HEADER TRUE, COMPRESSION GZIP);
    ",
    sq(prepared[[tr]]), sq(REF_TXT), sq(std_files[[tr]]))

    dbExecute(con, q)
  } else {
    message("[reuse existing] ", std_files[[tr]])
  }

  q2 <- sprintf("
    SELECT
      count(*) n,
      count(DISTINCT SNP) n_unique_snp,
      min(PVAL) min_p,
      max(abs(Z)) max_abs_z,
      median(N) median_n
    FROM read_csv_auto(
      '%s', delim='\\t', header=true, compression='gzip',
      sample_size=100000
    );
  ", sq(std_files[[tr]]))

  x <- as.data.table(dbGetQuery(con, q2))
  x[, trait := tr]
  x[, file := std_files[[tr]]]
  qc_std[[tr]] <- x
}

dbDisconnect(con, shutdown=TRUE)

std_qc <- rbindlist(qc_std, fill=TRUE)
setcolorder(std_qc, c("trait","n","n_unique_snp","min_p","max_abs_z","median_n","file"))
fwrite(std_qc, file.path(QC_DIR, "STEP5D2_BMI_OSA_standardized_QC.csv"))

# ---- official sumstats.py -> .mat ----
conv_log <- file.path(QC_DIR, "STEP5D2_BMI_OSA_mat_conversion.log")

for (tr in c("BMI","OSA")) {
  message("\n[MAT CONVERSION] ", tr)

  if (file.exists(mat_files[[tr]]) &&
      as.numeric(file.info(mat_files[[tr]])$size) > 100*1024^2) {
    message("[reuse existing .mat] ", mat_files[[tr]])
    next
  }

  out <- run_conda(c(
    "run","-n",CONDA_ENV,
    "python", normp(SUMSTATS_PY),
    "mat",
    "--sumstats", normp(std_files[[tr]]),
    "--ref", normp(REF_TXT),
    "--out", normp(mat_files[[tr]]),
    "--chunksize","500000",
    "--force"
  ), stdout=TRUE, stderr=TRUE)

  cat(
    paste0("\n========== ", tr, " ==========\n"),
    paste(out, collapse="\n"), "\n",
    file=conv_log, append=TRUE
  )

  if (!file.exists(mat_files[[tr]]) ||
      as.numeric(file.info(mat_files[[tr]])$size) < 100*1024^2) {
    stop(tr, " .mat conversion failed. Inspect ", conv_log)
  }
}

# ================= MAT QC =================
check_py <- file.path(TMP_DIR, "check_mat.py")
writeLines(c(
  "import sys,csv,os",
  "import numpy as np",
  "import scipy.io as sio",
  "out=sys.argv[1]",
  "rows=[]",
  "for arg in sys.argv[2:]:",
  "  trait,path=arg.split('=',1)",
  "  d=sio.loadmat(path, variable_names=['logpvec','zvec','nvec'])",
  "  lp=np.ravel(d['logpvec']); z=np.ravel(d['zvec']); n=np.ravel(d['nvec'])",
  "  rows.append({'trait':trait,'n_reference_slots':len(z),",
  "   'finite_z':int(np.isfinite(z).sum()),",
  "   'finite_logp':int(np.isfinite(lp).sum()),",
  "   'median_n':float(np.nanmedian(n)),",
  "   'max_abs_z':float(np.nanmax(np.abs(z))),",
  "   'mat_size_mb':os.path.getsize(path)/(1024**2),'path':path})",
  "with open(out,'w',newline='') as f:",
  "  w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)",
  "print('QC written:',out)"
), check_py)

mat_qc_csv <- file.path(QC_DIR, "STEP5D2_BMI_OSA_MAT_QC.csv")

qc_out <- run_conda(c(
  "run","-n",CONDA_ENV,
  "python", normp(check_py),
  normp(mat_qc_csv),
  paste0("BMI=",normp(mat_files["BMI"])),
  paste0("OSA=",normp(mat_files["OSA"]))
), stdout=TRUE, stderr=TRUE)

cat(paste(qc_out, collapse="\n"), "\n")
mat_qc <- fread(mat_qc_csv)
print(mat_qc)

if (any(mat_qc$n_reference_slots != 9545380L)) {
  stop("BMI/OSA MAT reference length mismatch. Do NOT run conjFDR.")
}
if (any(mat_qc$finite_z < 500000L)) {
  stop("Too few finite Z values after official reference matching.")
}

# ================= SELECT PAIR =================
pair_map <- list(
  AF_BMI = c("AF","BMI"),
  AF_OSA = c("AF","OSA"),
  HFpEF_BMI = c("HFpEF","BMI"),
  HFpEF_OSA = c("HFpEF","OSA"),
  BMI_OSA = c("BMI","OSA")
)

if (!PAIR_TO_RUN %in% names(pair_map)) {
  stop("PAIR_TO_RUN must be one of: ", paste(names(pair_map), collapse=", "))
}

tr <- pair_map[[PAIR_TO_RUN]]
T1 <- tr[1]
T2 <- tr[2]

if (!file.exists(mat_files[[T1]]) || !file.exists(mat_files[[T2]])) {
  stop("Required .mat missing for selected pair.")
}

RUN_DIR <- file.path(RUN_BASE, PAIR_TO_RUN)
RESULT_DIR <- file.path(RUN_DIR, "results")
dir.create(RESULT_DIR, recursive=TRUE, showWarnings=FALSE)

# ================= CONFIG =================
CONFIG_FILE <- file.path(PLEIO_REPO, "config.txt")

cfg <- c(
  paste0("reffile=", normp(REF_MAT)),
  paste0("traitfolder=", normp(MAT_DIR)),
  "",
  paste0("traitfile1=", T1, "_pleioFDR.mat"),
  paste0("traitname1=", T1),
  paste0("traitfiles={'", T2, "_pleioFDR.mat'}"),
  paste0("traitnames={'", T2, "'}"),
  "",
  paste0("outputdir=", normp(RESULT_DIR)),
  "stattype=conjfdr",
  paste0("fdrthresh=", CONJFDR_THRESHOLD),
  "",
  "randprune=true",
  paste0("randprune_n=", RANDPRUNE_N),
  paste0("exclude_chr_pos=[6 ", MHC_FROM, " ", MHC_TO, "]"),
  "",
  paste0("refinfo=", normp(REF_TXT)),
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
  "manh_plot=1",
  "exit_matlab_upon_completion=true"
)

writeLines(cfg, CONFIG_FILE)
file.copy(
  CONFIG_FILE,
  file.path(RUN_DIR, paste0("config_",PAIR_TO_RUN,".txt")),
  overwrite=TRUE
)

# ================= MATLAB =================
matlab_candidates <- unique(c(
  unname(Sys.which("matlab")),
  "E:/app/MATLAB/bin/matlab.exe",
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
MATLAB_EXE <- matlab_candidates[1]

# Readiness
ready <- data.table(
  item=c("pleioFDR runme","reference mat","reference txt",
         paste0(T1," mat"),paste0(T2," mat"),"MATLAB","config"),
  found=c(
    file.exists(file.path(PLEIO_REPO,"runme.m")),
    file.exists(REF_MAT),
    file.exists(REF_TXT),
    file.exists(mat_files[[T1]]),
    file.exists(mat_files[[T2]]),
    file.exists(MATLAB_EXE),
    file.exists(CONFIG_FILE)
  ),
  path=c(
    file.path(PLEIO_REPO,"runme.m"),
    REF_MAT,REF_TXT,
    mat_files[[T1]],mat_files[[T2]],
    MATLAB_EXE,CONFIG_FILE
  )
)
fwrite(ready, file.path(QC_DIR, paste0("STEP5D2_",PAIR_TO_RUN,"_readiness.csv")))
print(ready)

if (any(!ready$found)) stop("Readiness check failed.")

# ================= LAUNCH ONE PAIR ONLY =================
if (AUTO_START_MATLAB) {
  repo_matlab <- gsub("'", "''", normp(PLEIO_REPO))
  matlab_win <- normalizePath(MATLAB_EXE, winslash="\\", mustWork=TRUE)

  matlab_code <- paste0(
    "cd('", repo_matlab, "'); run('runme.m');"
  )

  cmd <- paste0(
    'start "', PAIR_TO_RUN, ' pleioFDR" cmd /k ""',
    matlab_win,
    '" -batch "',
    matlab_code,
    '" "'
  )

  cat("\nLaunching formal ", PAIR_TO_RUN, " conjFDR...\n", sep="")
  cat(cmd, "\n")
  shell(cmd, wait=FALSE)
}

cat("\n====================================================\n")
cat("STEP5D2 PREPARATION COMPLETE\n")
cat("Selected pair:", PAIR_TO_RUN, "\n")
cat("Results folder:", RESULT_DIR, "\n")
cat("When finished, ZIP this result folder and send it back.\n")
cat("====================================================\n")
