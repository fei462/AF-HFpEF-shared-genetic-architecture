# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11F_FIX_V6_PATCH_PANDAS_AND_RUN_AF_H2CTS(1).R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11F FIX V6 — patch LDSC pandas compatibility + rerun AF h2-cts
#
# CONFIRMED ERROR:
#   ldscore/sumstats.py, cell_type_specific()
#   pd.merge(...).loc[:,1:]
#
# FIX:
#   .loc[:,1:]  ->  .iloc[:,1:]
#
# RATIONALE:
#   After pd.merge(), columns have string labels.
#   We need positional selection "all columns except SNP",
#   therefore iloc is the correct pandas accessor.
#
# SAFETY:
#   - backs up sumstats.py before editing
#   - edits ONLY the line containing:
#       pd.merge(keep_snps, ref_ld_cts_allsnps
#   - checks Python syntax with py_compile
#   - does NOT recompute STEP11A-E
#   - does NOT recompute 264 cell-type LD scores
#   - does NOT recompute 22 AtrialUnion LD scores
#   - reruns ONLY AF --h2-cts
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

library(data.table)
library(ggplot2)

# ============================================================
# 0. PATHS — all confirmed by the V5 LDSC log
# ============================================================

DATA_ROOT <- "D:/A/data"

STEP11_ROOT <- file.path(
  DATA_ROOT,
  "STEP11_GSE238242"
)

RESULT_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "06_AF_H2CTS"
)

LOG_DIR <- file.path(
  STEP11_ROOT,
  "03_SLDSC",
  "07_AF_H2CTS_LOGS"
)

TABLE_DIR <- file.path(
  STEP11_ROOT,
  "05_FINAL",
  "SLDSC_AF"
)

FIG_DIR <- file.path(
  TABLE_DIR,
  "FIGURES"
)

for (d in c(
  RESULT_DIR,
  LOG_DIR,
  TABLE_DIR,
  FIG_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PYTHON_EXE <- "E:/app/python/envs/ldsc39/python.exe"
LDSC_PY <- "D:/A/ldsc39/ldsc.py"
SUMSTATS_PY <- "D:/A/ldsc39/ldscore/sumstats.py"

AF_SUMSTATS <- "D:/A/data/STEP1D_LDSC/01_munged/AF.sumstats.gz"

BASELINE_PREFIX <- paste0(
  "D:/A/data/CELLULAR/04_SLDSC/STEP10D1/02_LDSC_reference/",
  "baselineLD/baselineLD."
)

WEIGHTS_PREFIX <- paste0(
  "D:/A/data/CELLULAR/04_SLDSC/STEP10D1/02_LDSC_reference/",
  "weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
)

LDCTS_FILE <- file.path(
  RESULT_DIR,
  "GSE238242_LAA_ATAC_AF_primary_with_union_FIXED_V5.ldcts"
)

# ============================================================
# 1. HARD INPUT CHECK
# ============================================================

required_files <- c(
  python = PYTHON_EXE,
  ldsc = LDSC_PY,
  sumstats_py = SUMSTATS_PY,
  AF_sumstats = AF_SUMSTATS,
  ldcts = LDCTS_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files)) {
  stop(
    paste0(
      "Missing required file(s):\n",
      paste(
        names(missing_files),
        missing_files,
        sep = " = ",
        collapse = "\n"
      )
    )
  )
}

# Confirm the already-proven reference prefixes still exist.
check_ldscore_22 <- function(prefix) {
  all(
    file.exists(
      paste0(
        prefix,
        1:22,
        ".l2.ldscore.gz"
      )
    ) |
    file.exists(
      paste0(
        prefix,
        1:22,
        ".l2.ldscore"
      )
    )
  )
}

if (!check_ldscore_22(BASELINE_PREFIX)) {
  stop("baselineLD chr1-22 are not complete at the confirmed V5 prefix.")
}

if (!check_ldscore_22(WEIGHTS_PREFIX)) {
  stop("weights chr1-22 are not complete at the confirmed V5 prefix.")
}

# ============================================================
# 2. READ + AUDIT LDSC SOURCE BEFORE PATCH
# ============================================================

src <- readLines(
  SUMSTATS_PY,
  warn = FALSE,
  encoding = "UTF-8"
)

target_idx <- grep(
  "pd\\.merge\\(keep_snps,\\s*ref_ld_cts_allsnps",
  src,
  perl = TRUE
)

if (length(target_idx) != 1L) {

  audit <- data.table(
    line_number = seq_along(src),
    line = src
  )

  fwrite(
    audit[
      grepl(
        "ref_ld_cts|keep_snps|Performing regression",
        line
      )
    ],
    file.path(
      TABLE_DIR,
      "STEP11F_FIX_V6_sumstats_py_target_audit.csv"
    )
  )

  stop(
    paste0(
      "Expected exactly one cell_type_specific merge line, found ",
      length(target_idx),
      ".\nUpload STEP11F_FIX_V6_sumstats_py_target_audit.csv."
    )
  )
}

target_before <- src[target_idx]

cat(
  "\n============================================================\n",
  "LDSC SOURCE LINE BEFORE PATCH\n",
  "============================================================\n",
  "Line ",
  target_idx,
  ":\n",
  target_before,
  "\n",
  sep = ""
)

# ============================================================
# 3. BACK UP ORIGINAL SOURCE
# ============================================================

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

BACKUP_FILE <- paste0(
  SUMSTATS_PY,
  ".BACKUP_BEFORE_H2CTS_PANDAS_FIX_",
  timestamp
)

ok_backup <- file.copy(
  SUMSTATS_PY,
  BACKUP_FILE,
  overwrite = FALSE
)

if (!isTRUE(ok_backup) || !file.exists(BACKUP_FILE)) {
  stop(
    "Could not back up sumstats.py. No source edit was made."
  )
}

cat(
  "\nBackup created:\n",
  BACKUP_FILE,
  "\n",
  sep = ""
)

# ============================================================
# 4. PATCH ONLY THE CELL-TYPE-SPECIFIC MERGE LINE
# ============================================================

patched_line <- target_before

# Current CBIIT/Python3 code path seen in the user's traceback.
patched_line <- gsub(
  "\\.loc\\[:,\\s*1:\\]",
  ".iloc[:, 1:]",
  patched_line,
  perl = TRUE
)

# Defensive support if a legacy source uses .ix on this exact same line.
patched_line <- gsub(
  "\\.ix\\[:,\\s*1:\\]",
  ".iloc[:, 1:]",
  patched_line,
  perl = TRUE
)

if (
  identical(patched_line, target_before) &&
  !grepl("\\.iloc\\[:,\\s*1:\\]", target_before, perl = TRUE)
) {
  stop(
    paste0(
      "The target line was found, but it contains neither ",
      ".loc[:,1:], .ix[:,1:], nor the already-fixed .iloc[:,1:].\n",
      "No edit was made."
    )
  )
}

already_patched <- identical(
  patched_line,
  target_before
) &&
  grepl(
    "\\.iloc\\[:,\\s*1:\\]",
    target_before,
    perl = TRUE
  )

if (!already_patched) {

  src[target_idx] <- patched_line

  writeLines(
    src,
    SUMSTATS_PY,
    useBytes = TRUE
  )

  cat(
    "\nSource patched successfully.\n"
  )

} else {

  cat(
    "\nTarget source line was already patched; no additional edit required.\n"
  )
}

# Re-read from disk and verify.
src_after <- readLines(
  SUMSTATS_PY,
  warn = FALSE,
  encoding = "UTF-8"
)

target_after <- src_after[target_idx]

if (!grepl(
  "\\.iloc\\[:,\\s*1:\\]",
  target_after,
  perl = TRUE
)) {
  stop(
    paste0(
      "Patch verification failed.\n",
      "Restore backup if needed:\n",
      BACKUP_FILE
    )
  )
}

if (grepl(
  "\\.(loc|ix)\\[:,\\s*1:\\]",
  target_after,
  perl = TRUE
)) {
  stop(
    "Patch verification found an old .loc/.ix positional slice still present."
  )
}

cat(
  "\n============================================================\n",
  "LDSC SOURCE LINE AFTER PATCH\n",
  "============================================================\n",
  "Line ",
  target_idx,
  ":\n",
  target_after,
  "\n",
  sep = ""
)

PATCH_AUDIT <- data.table(
  field = c(
    "sumstats_py",
    "backup_file",
    "target_line_number",
    "before",
    "after",
    "already_patched_before_this_run"
  ),
  value = c(
    SUMSTATS_PY,
    BACKUP_FILE,
    as.character(target_idx),
    target_before,
    target_after,
    as.character(already_patched)
  )
)

fwrite(
  PATCH_AUDIT,
  file.path(
    TABLE_DIR,
    "STEP11F_FIX_V6_pandas_patch_audit.csv"
  )
)

# ============================================================
# 5. PYTHON SYNTAX CHECK
# ============================================================

PYCOMPILE_LOG <- file.path(
  LOG_DIR,
  "STEP11F_FIX_V6_py_compile.log"
)

compile_status <- suppressWarnings(
  system2(
    PYTHON_EXE,
    args = c(
      "-m",
      "py_compile",
      SUMSTATS_PY
    ),
    stdout = PYCOMPILE_LOG,
    stderr = PYCOMPILE_LOG
  )
)

if (!identical(
  as.integer(compile_status),
  0L
)) {

  cat(
    "\nPython syntax check failed.\n"
  )

  if (file.exists(PYCOMPILE_LOG)) {
    cat(
      paste(
        readLines(
          PYCOMPILE_LOG,
          warn = FALSE
        ),
        collapse = "\n"
      ),
      "\n"
    )
  }

  stop(
    paste0(
      "Patched sumstats.py failed py_compile.\n",
      "Original backup remains at:\n",
      BACKUP_FILE
    )
  )
}

cat(
  "\nPython py_compile: PASS\n"
)

# ============================================================
# 6. RUN ONLY AF --h2-cts
# ============================================================

AF_OUT <- file.path(
  RESULT_DIR,
  "AF_GSE238242_LAA_ATAC_FIXED_V6"
)

CONSOLE_LOG <- file.path(
  LOG_DIR,
  "STEP11F_FIX_V6_AF_h2_cts_console.log"
)

RESULT_FILE <- paste0(
  AF_OUT,
  ".cell_type_results.txt"
)

# Remove only stale V6 outputs.
stale <- c(
  CONSOLE_LOG,
  RESULT_FILE,
  paste0(AF_OUT, ".log")
)

unlink(
  stale[
    file.exists(stale)
  ]
)

args_h2cts <- c(
  LDSC_PY,
  "--h2-cts",
  AF_SUMSTATS,
  "--ref-ld-chr",
  BASELINE_PREFIX,
  "--ref-ld-chr-cts",
  LDCTS_FILE,
  "--w-ld-chr",
  WEIGHTS_PREFIX,
  "--out",
  AF_OUT
)

cat(
  "\n============================================================\n",
  "RUNNING ONLY AF --h2-cts — FIX V6\n",
  "============================================================\n\n",
  "No LD scores are being recomputed.\n\n",
  "AF: ",
  AF_SUMSTATS,
  "\nBaseline: ",
  BASELINE_PREFIX,
  "\nWeights: ",
  WEIGHTS_PREFIX,
  "\nCTS: ",
  LDCTS_FILE,
  "\n\n",
  sep = ""
)

exit_code <- suppressWarnings(
  system2(
    PYTHON_EXE,
    args = args_h2cts,
    stdout = CONSOLE_LOG,
    stderr = CONSOLE_LOG
  )
)

# ============================================================
# 7. FAILURE DIAGNOSTIC
# ============================================================

if (
  !identical(
    as.integer(exit_code),
    0L
  ) ||
  !file.exists(RESULT_FILE)
) {

  cat(
    "\nAF h2-cts V6 failed AFTER the pandas patch.\n"
  )

  if (file.exists(CONSOLE_LOG)) {

    z <- readLines(
      CONSOLE_LOG,
      warn = FALSE
    )

    cat(
      "\nLAST 120 LOG LINES:\n\n",
      paste(
        tail(z, 120),
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )
  }

  stop(
    paste0(
      "STEP11F FIX V6 stopped during AF h2-cts.\n",
      "Upload ONLY:\n",
      CONSOLE_LOG,
      "\n",
      file.path(
        TABLE_DIR,
        "STEP11F_FIX_V6_pandas_patch_audit.csv"
      )
    )
  )
}

# ============================================================
# 8. PARSE FORMAL RESULTS
# ============================================================

RES <- fread(
  RESULT_FILE,
  showProgress = FALSE
)

required_cols <- c(
  "Name",
  "Coefficient",
  "Coefficient_std_error",
  "Coefficient_P_value"
)

if (!all(
  required_cols %in% names(RES)
)) {
  stop(
    paste0(
      "Unexpected h2-cts output columns:\n",
      paste(
        names(RES),
        collapse = ", "
      )
    )
  )
}

if (nrow(RES) != 12L) {
  stop(
    "Expected 12 cell-type tests; observed ",
    nrow(RES)
  )
}

RES[
  ,
  Coefficient_Z :=
    Coefficient /
    Coefficient_std_error
]

RES[
  ,
  FDR_BH :=
    p.adjust(
      Coefficient_P_value,
      method = "BH"
    )
]

RES[
  ,
  Bonferroni_P :=
    p.adjust(
      Coefficient_P_value,
      method = "bonferroni"
    )
]

RES[
  ,
  FDR05 :=
    FDR_BH < 0.05
]

RES[
  ,
  Bonferroni05 :=
    Bonferroni_P < 0.05
]

RES[
  ,
  rank_by_P :=
    frank(
      Coefficient_P_value,
      ties.method = "min"
    )
]

setorder(
  RES,
  Coefficient_P_value
)

FINAL_RESULTS <- file.path(
  TABLE_DIR,
  "STEP11F_AF_GSE238242_celltype_S_LDSC_results_FIXED_V6.csv"
)

fwrite(
  RES,
  FINAL_RESULTS
)

fwrite(
  RES[FDR_BH < 0.05],
  file.path(
    TABLE_DIR,
    "STEP11F_AF_FDR05_celltypes_FIXED_V6.csv"
  )
)

fwrite(
  RES[Bonferroni_P < 0.05],
  file.path(
    TABLE_DIR,
    "STEP11F_AF_Bonferroni05_celltypes_FIXED_V6.csv"
  )
)

# ============================================================
# 9. FIGURE SOURCE + PUBLICATION FIGURE
# ============================================================

FIG_SOURCE <- copy(RES)

FIG_SOURCE[
  ,
  significance :=
    fifelse(
      Bonferroni05,
      "Bonferroni < 0.05",
      fifelse(
        FDR05,
        "FDR < 0.05",
        "Not significant"
      )
    )
]

fwrite(
  FIG_SOURCE,
  file.path(
    TABLE_DIR,
    "STEP11F_AF_S_LDSC_figure_source_FIXED_V6.csv"
  )
)

FIG_SOURCE[
  ,
  Name :=
    factor(
      Name,
      levels = rev(
        Name[
          order(
            Coefficient_Z,
            decreasing = TRUE
          )
        ]
      )
    )
]

p <- ggplot(
  FIG_SOURCE,
  aes(
    x = Coefficient_Z,
    y = Name
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.45,
    colour = "grey55"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = Coefficient_Z,
      yend = Name
    ),
    linewidth = 0.65,
    colour = "grey78"
  ) +
  geom_point(
    aes(
      fill = significance
    ),
    shape = 21,
    size = 4.2,
    stroke = 0.55,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Bonferroni < 0.05" = "#B2182B",
      "FDR < 0.05" = "#EF8A62",
      "Not significant" = "#4D648D"
    )
  ) +
  labs(
    x = expression(tau / SE),
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

TIFF_FILE <- file.path(
  FIG_DIR,
  "Figure_AF_LAA_celltype_S_LDSC_tauZ_FIXED_V6.tiff"
)

PDF_FILE <- file.path(
  FIG_DIR,
  "Figure_AF_LAA_celltype_S_LDSC_tauZ_FIXED_V6.pdf"
)

ggsave(
  TIFF_FILE,
  p,
  width = 5.8,
  height = 5.3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  PDF_FILE,
  p,
  width = 5.8,
  height = 5.3,
  units = "in"
)

# ============================================================
# 10. FINAL READINESS
# ============================================================

READINESS <- data.table(
  check = c(
    "sumstats_py_backup_created",
    "cell_type_specific_merge_uses_iloc",
    "python_py_compile_pass",
    "baselineLD_22_present",
    "weights_22_present",
    "AF_h2_cts_exit_zero",
    "cell_type_results_created",
    "12_celltype_tests_present",
    "publication_figure_created"
  ),
  pass = c(
    file.exists(BACKUP_FILE),
    grepl(
      "\\.iloc\\[:,\\s*1:\\]",
      target_after,
      perl = TRUE
    ),
    identical(
      as.integer(compile_status),
      0L
    ),
    check_ldscore_22(BASELINE_PREFIX),
    check_ldscore_22(WEIGHTS_PREFIX),
    identical(
      as.integer(exit_code),
      0L
    ),
    file.exists(RESULT_FILE),
    nrow(RES) == 12L,
    file.exists(TIFF_FILE)
  )
)

fwrite(
  READINESS,
  file.path(
    TABLE_DIR,
    "STEP11F_FIX_V6_readiness.csv"
  )
)

METHOD <- data.table(
  field = c(
    "root_cause",
    "source_fix",
    "source_backup",
    "upstream_LD_scores_recomputed",
    "formal_test",
    "control",
    "multiple_testing"
  ),
  value = c(
    "pandas label-based .loc was used for positional column slicing in LDSC cell_type_specific()",
    "replace target merge slice .loc[:,1:] / .ix[:,1:] with .iloc[:,1:]",
    BACKUP_FILE,
    "No",
    "LDSC --h2-cts",
    "AtrialUnion12",
    "BH-FDR + Bonferroni across 12 pre-specified cell types"
  )
)

fwrite(
  METHOD,
  file.path(
    TABLE_DIR,
    "STEP11F_FIX_V6_method_provenance.csv"
  )
)

# ============================================================
# 11. CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================================\n",
  "STEP11F FIX V6 COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

print(
  RES[
    ,
    .(
      Name,
      Coefficient,
      Coefficient_std_error,
      Coefficient_Z,
      Coefficient_P_value,
      FDR_BH,
      Bonferroni_P
    )
  ]
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

cat(
  "\nUPLOAD:\n",
  TABLE_DIR,
  "\n\nMost important files:\n",
  "STEP11F_FIX_V6_readiness.csv\n",
  "STEP11F_FIX_V6_pandas_patch_audit.csv\n",
  "STEP11F_AF_GSE238242_celltype_S_LDSC_results_FIXED_V6.csv\n",
  "STEP11F_AF_FDR05_celltypes_FIXED_V6.csv\n",
  "STEP11F_AF_Bonferroni05_celltypes_FIXED_V6.csv\n",
  "STEP11F_AF_S_LDSC_figure_source_FIXED_V6.csv\n",
  sep = ""
)
