# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP11E_GSE238242_BUILD_LDSC_ANNOT_AND_LDSCORES.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP11E — GSE238242 cell-type ATAC -> LDSC annotations + LD scores
# Project: AF-priority left atrial cell-state genetics
#
# INPUTS
#   STEP11D hg19 merged BED files (12 author cell types)
#   1000G Phase 3 EUR PLINK files: 1000G.EUR.QC.{1..22}.{bed,bim,fam}
#   HapMap3 SNP list: hm3_no_MHC.list.txt
#   Existing LDSC Python environment/source from the earlier LDSC analysis
#
# OUTPUTS
#   12 x 22 thin .annot.gz files
#   12 x 22 custom .l2.ldscore.gz files + M / M_5_50
#   one .ldcts file for STEP11F AF cell-type S-LDSC
#
# IMPORTANT METHOD RULES
#   - Primary annotation remains STEP11C >=10% accessible cells.
#   - Coordinates are the frozen STEP11D hg19 unique/autosomal intervals.
#   - Annotation vectors are generated in EXACT BIM SNP order.
#   - LD scores are computed from the same 1000G EUR Phase 3 PLINK panel.
#   - HapMap3 SNPs are used for --print-snps.
#   - No biological result is interpreted in this step.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. USER SETTINGS
# ============================================================

DATA_ROOT <- "D:/A/data"
PROJECT_ROOT <- file.path(DATA_ROOT, "STEP11_GSE238242")

# TRUE = after constructing annotations, immediately compute all 264 LD-score jobs.
# The script is checkpointed: rerunning will skip completed LD-score outputs.
RUN_LDSC_NOW <- TRUE

# ============================================================
# 1. DIRECTORIES
# ============================================================

LIFT_DIR <- file.path(
  PROJECT_ROOT,
  "03_SLDSC",
  "02_LIFTOVER_GRCh37"
)

LD_DIR <- file.path(
  PROJECT_ROOT,
  "03_SLDSC",
  "03_CELLTYPE_LD_SCORES"
)

HM3_DIR <- file.path(
  PROJECT_ROOT,
  "03_SLDSC",
  "02B_HM3_BY_CHR"
)

QC_DIR <- file.path(
  PROJECT_ROOT,
  "00_QC"
)

LOG_DIR <- file.path(
  PROJECT_ROOT,
  "03_SLDSC",
  "04_LDSC_LOGS"
)

for (d in c(
  LD_DIR,
  HM3_DIR,
  QC_DIR,
  LOG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 2. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "GenomicRanges",
  "IRanges"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  stop(
    "Missing R package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(data.table)
library(GenomicRanges)
library(IRanges)

# ============================================================
# 3. CELL TYPES + STEP11D HARD CHECK
# ============================================================

CELL_TYPES <- c(
  "Adipo",
  "CM",
  "EC",
  "Endocardial",
  "FB",
  "Lymphoid",
  "Mast",
  "Mesothelial",
  "Myeloid",
  "Neuron",
  "PC",
  "SMC"
)

READINESS_D <- file.path(
  LIFT_DIR,
  "STEP11D_readiness.csv"
)

if (!file.exists(READINESS_D)) {
  stop(
    "Missing STEP11D readiness file: ",
    READINESS_D
  )
}

D_READY <- fread(
  READINESS_D
)

if (
  !"pass" %in% names(D_READY) ||
  !all(
    as.logical(
      D_READY$pass
    )
  )
) {
  stop(
    "STEP11D readiness is not fully PASS."
  )
}

BED_FILES <- setNames(
  file.path(
    LIFT_DIR,
    paste0(
      "GSE238242_",
      CELL_TYPES,
      "_ATAC_hg19_merged.bed.gz"
    )
  ),
  CELL_TYPES
)

missing_beds <- BED_FILES[
  !file.exists(BED_FILES)
]

if (length(missing_beds)) {
  stop(
    paste0(
      "Missing STEP11D BED file(s):\n",
      paste(missing_beds, collapse = "\n")
    )
  )
}

# ============================================================
# 4. AUTO-DETECT 1000G EUR PLINK PREFIX
# ============================================================

prefix_triplet_ok <- function(prefix, chr) {

  all(
    file.exists(
      paste0(
        prefix,
        chr,
        c(
          ".bed",
          ".bim",
          ".fam"
        )
      )
    )
  )
}

candidate_chr1_bims <- unique(
  c(
    # Most likely locations from this project / prior LDSC work.
    Sys.glob(
      file.path(
        DATA_ROOT,
        "STEP15_HFPEF_PREP",
        "**",
        "1000G.EUR.QC.1.bim"
      )
    ),
    Sys.glob(
      file.path(
        DATA_ROOT,
        "STEP1D_LDSC",
        "**",
        "1000G.EUR.QC.1.bim"
      )
    )
  )
)

# Base R Sys.glob does not guarantee recursive ** behavior on every Windows setup.
# Add an explicit recursive scan if necessary.
if (!length(candidate_chr1_bims)) {

  scan_roots <- c(
    file.path(DATA_ROOT, "STEP15_HFPEF_PREP"),
    file.path(DATA_ROOT, "STEP1D_LDSC"),
    DATA_ROOT
  )

  for (rr in scan_roots[dir.exists(scan_roots)]) {

    z <- tryCatch(
      list.files(
        rr,
        pattern = "^1000G[.]EUR[.]QC[.]1[.]bim$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = FALSE
      ),
      error = function(e) character(0)
    )

    candidate_chr1_bims <- unique(
      c(
        candidate_chr1_bims,
        z
      )
    )

    if (length(candidate_chr1_bims)) {
      break
    }
  }
}

if (!length(candidate_chr1_bims)) {
  stop(
    paste0(
      "Could not locate 1000G.EUR.QC.1.bim under D:/A/data.\n",
      "The 1000G Phase 3 PLINK reference must be extracted before STEP11E."
    )
  )
}

candidate_prefixes <- unique(
  sub(
    "1[.]bim$",
    "",
    normalizePath(
      candidate_chr1_bims,
      winslash = "/",
      mustWork = TRUE
    )
  )
)

prefix_complete <- vapply(
  candidate_prefixes,
  function(pre) {
    all(
      vapply(
        1:22,
        function(chr) {
          prefix_triplet_ok(
            pre,
            chr
          )
        },
        logical(1)
      )
    )
  },
  logical(1)
)

candidate_prefixes <- candidate_prefixes[
  prefix_complete
]

if (!length(candidate_prefixes)) {
  stop(
    "1000G candidate(s) were found, but no prefix has complete BED/BIM/FAM files for chr1-22."
  )
}

# Prefer previous project reference paths, then shortest path.
score_prefix <- function(x) {

  1000L * grepl(
    "STEP15_HFPEF_PREP",
    x,
    ignore.case = TRUE
  ) +
  500L * grepl(
    "STEP1D_LDSC",
    x,
    ignore.case = TRUE
  ) +
  100L * grepl(
    "1000G",
    x,
    ignore.case = TRUE
  ) -
  nchar(x) / 1000
}

REF_PREFIX <- candidate_prefixes[
  which.max(
    vapply(
      candidate_prefixes,
      score_prefix,
      numeric(1)
    )
  )
]

cat(
  "\n1000G EUR PLINK prefix selected:\n",
  REF_PREFIX,
  "\n",
  sep = ""
)

# ============================================================
# 5. AUTO-DETECT HapMap3 SNP LIST
# ============================================================

HM3_CANDIDATES <- c(
  file.path(
    DATA_ROOT,
    "hm3_no_MHC.list.txt"
  ),
  file.path(
    DATA_ROOT,
    "hm3_no_MHC.list"
  ),
  file.path(
    DATA_ROOT,
    "w_hm3.snplist"
  )
)

HM3_CANDIDATES <- HM3_CANDIDATES[
  file.exists(
    HM3_CANDIDATES
  )
]

if (!length(HM3_CANDIDATES)) {

  hm_scan <- tryCatch(
    list.files(
      DATA_ROOT,
      pattern = "hm3.*(list|snplist|snp).*",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    ),
    error = function(e) character(0)
  )

  HM3_CANDIDATES <- hm_scan[
    file.exists(
      hm_scan
    )
  ]
}

if (!length(HM3_CANDIDATES)) {
  stop(
    "Could not locate a HapMap3 SNP list."
  )
}

HM3_FILE <- normalizePath(
  HM3_CANDIDATES[1],
  winslash = "/",
  mustWork = TRUE
)

HM3_RAW <- fread(
  HM3_FILE,
  header = FALSE,
  select = 1,
  showProgress = FALSE
)

HM3_SNPS <- unique(
  as.character(
    HM3_RAW[[1]]
  )
)

HM3_SNPS <- HM3_SNPS[
  nzchar(HM3_SNPS)
]

if (length(HM3_SNPS) < 500000L) {
  warning(
    "The detected HapMap3 list contains only ",
    length(HM3_SNPS),
    " entries. Inspect before STEP11F if unexpected."
  )
}

cat(
  "\nHapMap3 SNP list:\n",
  HM3_FILE,
  "\nEntries: ",
  format(length(HM3_SNPS), big.mark = ","),
  "\n",
  sep = ""
)

# ============================================================
# 6. AUTO-DETECT LDSC SOURCE + PYTHON ENVIRONMENT
# ============================================================

LDSC_CANDIDATES <- c(
  "D:/A/ldsc39/ldsc.py",
  file.path(
    DATA_ROOT,
    "ldsc39",
    "ldsc.py"
  )
)

LDSC_CANDIDATES <- LDSC_CANDIDATES[
  file.exists(
    LDSC_CANDIDATES
  )
]

if (!length(LDSC_CANDIDATES)) {

  ldsc_scan_roots <- c(
    "D:/A",
    DATA_ROOT
  )

  for (rr in ldsc_scan_roots[dir.exists(ldsc_scan_roots)]) {

    z <- tryCatch(
      list.files(
        rr,
        pattern = "^ldsc[.]py$",
        recursive = TRUE,
        full.names = TRUE
      ),
      error = function(e) character(0)
    )

    if (length(z)) {
      LDSC_CANDIDATES <- unique(
        c(
          LDSC_CANDIDATES,
          z
        )
      )
      break
    }
  }
}

if (!length(LDSC_CANDIDATES)) {
  stop(
    "Could not locate ldsc.py. Expected prior LDSC source under D:/A/ldsc39."
  )
}

LDSC_PY <- normalizePath(
  LDSC_CANDIDATES[1],
  winslash = "/",
  mustWork = TRUE
)

PYTHON_CANDIDATES <- unique(
  c(
    "E:/app/python/envs/ldsc39/python.exe",
    "D:/A/python/envs/ldsc39/python.exe",
    unname(
      Sys.which(
        "python"
      )
    )
  )
)

PYTHON_CANDIDATES <- PYTHON_CANDIDATES[
  nzchar(PYTHON_CANDIDATES) &
  file.exists(PYTHON_CANDIDATES)
]

if (!length(PYTHON_CANDIDATES)) {
  stop(
    "Could not locate the Python executable for the LDSC environment."
  )
}

# Prefer paths explicitly containing ldsc39.
PYTHON_EXE <- PYTHON_CANDIDATES[
  order(
    !grepl(
      "ldsc39",
      PYTHON_CANDIDATES,
      ignore.case = TRUE
    ),
    nchar(PYTHON_CANDIDATES)
  )
][1]

PYTHON_EXE <- normalizePath(
  PYTHON_EXE,
  winslash = "/",
  mustWork = TRUE
)

cat(
  "\nLDSC Python:\n",
  PYTHON_EXE,
  "\nLDSC script:\n",
  LDSC_PY,
  "\n",
  sep = ""
)

# ============================================================
# 7. LOAD CELL-TYPE BED FILES ONCE
#    BED is 0-based / end-exclusive.
#    Convert to 1-based GRanges: start = BED_start + 1, end = BED_end.
# ============================================================

CT_GR <- list()
BED_QC_LIST <- list()

for (ct in CELL_TYPES) {

  bed <- fread(
    BED_FILES[[ct]],
    header = FALSE,
    col.names = c(
      "chr",
      "start0",
      "end0",
      "name"
    ),
    showProgress = FALSE
  )

  if (!nrow(bed)) {
    stop(
      "Empty BED for cell type: ",
      ct
    )
  }

  bed[
    ,
    start1 :=
      as.integer(start0) + 1L
  ]

  bed[
    ,
    end1 :=
      as.integer(end0)
  ]

  valid <- (
    bed$start1 >= 1L &
    bed$end1 >= bed$start1 &
    grepl(
      "^chr([1-9]|1[0-9]|2[0-2])$",
      bed$chr
    )
  )

  if (!all(valid)) {
    stop(
      "Invalid hg19 BED coordinates detected for ",
      ct
    )
  }

  CT_GR[[ct]] <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(
      start = bed$start1,
      end = bed$end1
    )
  )

  BED_QC_LIST[[
    length(BED_QC_LIST) + 1L
  ]] <- data.table(
    cell_type = ct,
    n_intervals = nrow(bed),
    n_chromosomes = uniqueN(bed$chr),
    min_start = min(bed$start1),
    max_end = max(bed$end1)
  )
}

BED_QC <- rbindlist(
  BED_QC_LIST
)

fwrite(
  BED_QC,
  file.path(
    QC_DIR,
    "STEP11E_input_BED_QC.csv"
  )
)

# ============================================================
# 8. BUILD THIN .annot.gz FILES IN EXACT BIM ORDER
# ============================================================

ANNOT_QC_LIST <- list()
HM3_QC_LIST <- list()

for (chr in 1:22) {

  cat(
    "\n============================================================\n",
    "Building annotations for chromosome ",
    chr,
    "\n",
    "============================================================\n",
    sep = ""
  )

  BIM_FILE <- paste0(
    REF_PREFIX,
    chr,
    ".bim"
  )

  BIM <- fread(
    BIM_FILE,
    header = FALSE,
    col.names = c(
      "CHR",
      "SNP",
      "CM",
      "BP",
      "A1",
      "A2"
    ),
    showProgress = FALSE
  )

  if (!nrow(BIM)) {
    stop(
      "Empty BIM: ",
      BIM_FILE
    )
  }

  if (!all(
    as.integer(BIM$CHR) == chr
  )) {
    stop(
      "Unexpected chromosome values in ",
      BIM_FILE
    )
  }

  SNP_GR <- GRanges(
    seqnames = paste0(
      "chr",
      chr
    ),
    ranges = IRanges(
      start = as.integer(
        BIM$BP
      ),
      width = 1L
    )
  )

  # Per-chromosome HapMap3 list, preserving BIM order.
  HM3_CHR <- BIM[
    SNP %chin% HM3_SNPS,
    SNP
  ]

  if (length(HM3_CHR) < 1000L) {
    warning(
      "Only ",
      length(HM3_CHR),
      " HapMap3 SNPs on chr",
      chr
    )
  }

  HM3_CHR_FILE <- file.path(
    HM3_DIR,
    paste0(
      "hm3_no_MHC.",
      chr,
      ".snp"
    )
  )

  fwrite(
    data.table(
      SNP = HM3_CHR
    ),
    HM3_CHR_FILE,
    col.names = FALSE,
    sep = "\t"
  )

  HM3_QC_LIST[[
    length(HM3_QC_LIST) + 1L
  ]] <- data.table(
    chr = chr,
    n_bim_snps = nrow(BIM),
    n_hm3_snps = length(HM3_CHR),
    hm3_fraction = length(HM3_CHR) / nrow(BIM)
  )

  for (ct in CELL_TYPES) {

    ct_chr <- CT_GR[[ct]][
      as.character(
        seqnames(
          CT_GR[[ct]]
        )
      ) == paste0(
        "chr",
        chr
      )
    ]

    if (!length(ct_chr)) {

      annot <- integer(
        nrow(BIM)
      )

    } else {

      hit <- findOverlaps(
        SNP_GR,
        ct_chr,
        ignore.strand = TRUE,
        select = "first"
      )

      annot <- as.integer(
        !is.na(hit)
      )
    }

    if (length(annot) != nrow(BIM)) {
      stop(
        "Annotation length mismatch for ",
        ct,
        " chr",
        chr
      )
    }

    ANNOT_FILE <- file.path(
      LD_DIR,
      paste0(
        "GSE238242_",
        ct,
        ".",
        chr,
        ".annot.gz"
      )
    )

    fwrite(
      data.table(
        ANNOT = annot
      ),
      ANNOT_FILE,
      sep = "\t",
      compress = "gzip"
    )

    ANNOT_QC_LIST[[
      length(ANNOT_QC_LIST) + 1L
    ]] <- data.table(
      cell_type = ct,
      chr = chr,
      n_bim_snps = nrow(BIM),
      n_annotated_snps = sum(annot),
      annotation_fraction = mean(annot),
      annot_file = normalizePath(
        ANNOT_FILE,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }
}

ANNOT_QC <- rbindlist(
  ANNOT_QC_LIST,
  fill = TRUE
)

HM3_QC <- rbindlist(
  HM3_QC_LIST,
  fill = TRUE
)

fwrite(
  ANNOT_QC,
  file.path(
    QC_DIR,
    "STEP11E_annotation_QC_by_celltype_chr.csv"
  )
)

fwrite(
  HM3_QC,
  file.path(
    QC_DIR,
    "STEP11E_HM3_QC_by_chr.csv"
  )
)

# Cell-type global annotation summary.
ANNOT_SUMMARY <- ANNOT_QC[
  ,
  .(
    n_reference_snps = sum(
      n_bim_snps
    ),
    n_annotated_snp_chr_sum = sum(
      n_annotated_snps
    ),
    min_chr_annotation_fraction = min(
      annotation_fraction
    ),
    max_chr_annotation_fraction = max(
      annotation_fraction
    ),
    mean_chr_annotation_fraction = weighted.mean(
      annotation_fraction,
      n_bim_snps
    )
  ),
  by = cell_type
]

fwrite(
  ANNOT_SUMMARY,
  file.path(
    QC_DIR,
    "STEP11E_annotation_summary.csv"
  )
)

# ============================================================
# 9. COMPUTE CUSTOM LD SCORES
# ============================================================

L2_STATUS_LIST <- list()

if (RUN_LDSC_NOW) {

  total_jobs <- length(
    CELL_TYPES
  ) * 22L

  job_no <- 0L

  for (ct in CELL_TYPES) {

    for (chr in 1:22) {

      job_no <- job_no + 1L

      cat(
        "\n[",
        job_no,
        "/",
        total_jobs,
        "] ",
        ct,
        " chr",
        chr,
        "\n",
        sep = ""
      )

      BFILE_PREFIX <- paste0(
        REF_PREFIX,
        chr
      )

      ANNOT_FILE <- file.path(
        LD_DIR,
        paste0(
          "GSE238242_",
          ct,
          ".",
          chr,
          ".annot.gz"
        )
      )

      HM3_CHR_FILE <- file.path(
        HM3_DIR,
        paste0(
          "hm3_no_MHC.",
          chr,
          ".snp"
        )
      )

      OUT_PREFIX <- file.path(
        LD_DIR,
        paste0(
          "GSE238242_",
          ct,
          ".",
          chr
        )
      )

      EXPECTED_OUTPUTS <- c(
        paste0(
          OUT_PREFIX,
          ".l2.ldscore.gz"
        ),
        paste0(
          OUT_PREFIX,
          ".l2.M"
        ),
        paste0(
          OUT_PREFIX,
          ".l2.M_5_50"
        )
      )

      if (all(
        file.exists(
          EXPECTED_OUTPUTS
        )
      )) {

        cat(
          "  Existing LD-score outputs detected; skipping.\n"
        )

        L2_STATUS_LIST[[
          length(L2_STATUS_LIST) + 1L
        ]] <- data.table(
          cell_type = ct,
          chr = chr,
          status = "SKIP_COMPLETE",
          exit_code = 0L,
          output_prefix = normalizePath(
            OUT_PREFIX,
            winslash = "/",
            mustWork = FALSE
          )
        )

        next
      }

      LOG_FILE <- file.path(
        LOG_DIR,
        paste0(
          "STEP11E_",
          ct,
          "_chr",
          chr,
          ".log"
        )
      )

      args <- c(
        shQuote(
          LDSC_PY
        ),
        "--l2",
        "--bfile",
        shQuote(
          BFILE_PREFIX
        ),
        "--ld-wind-cm",
        "1",
        "--annot",
        shQuote(
          normalizePath(
            ANNOT_FILE,
            winslash = "/",
            mustWork = TRUE
          )
        ),
        "--thin-annot",
        "--out",
        shQuote(
          normalizePath(
            OUT_PREFIX,
            winslash = "/",
            mustWork = FALSE
          )
        ),
        "--print-snps",
        shQuote(
          normalizePath(
            HM3_CHR_FILE,
            winslash = "/",
            mustWork = TRUE
          )
        )
      )

      exit_code <- suppressWarnings(
        system2(
          PYTHON_EXE,
          args = args,
          stdout = LOG_FILE,
          stderr = LOG_FILE
        )
      )

      output_ok <- all(
        file.exists(
          EXPECTED_OUTPUTS
        )
      )

      status <- if (
        identical(
          as.integer(exit_code),
          0L
        ) &&
        output_ok
      ) {
        "PASS"
      } else {
        "FAIL"
      }

      L2_STATUS_LIST[[
        length(L2_STATUS_LIST) + 1L
      ]] <- data.table(
        cell_type = ct,
        chr = chr,
        status = status,
        exit_code = as.integer(
          exit_code
        ),
        output_prefix = normalizePath(
          OUT_PREFIX,
          winslash = "/",
          mustWork = FALSE
        ),
        log_file = normalizePath(
          LOG_FILE,
          winslash = "/",
          mustWork = FALSE
        )
      )

      if (status == "FAIL") {

        L2_STATUS_NOW <- rbindlist(
          L2_STATUS_LIST,
          fill = TRUE
        )

        fwrite(
          L2_STATUS_NOW,
          file.path(
            QC_DIR,
            "STEP11E_LDscore_job_status.csv"
          )
        )

        stop(
          paste0(
            "LDSC LD-score computation failed for ",
            ct,
            " chr",
            chr,
            ".\nInspect log:\n",
            LOG_FILE
          )
        )
      }
    }
  }
}

if (length(L2_STATUS_LIST)) {

  L2_STATUS <- rbindlist(
    L2_STATUS_LIST,
    fill = TRUE
  )

} else {

  L2_STATUS <- data.table(
    cell_type = character(),
    chr = integer(),
    status = character(),
    exit_code = integer()
  )
}

fwrite(
  L2_STATUS,
  file.path(
    QC_DIR,
    "STEP11E_LDscore_job_status.csv"
  )
)

# ============================================================
# 10. CREATE .ldcts FILE FOR STEP11F
# ============================================================

LD_PREFIX_ROOT <- normalizePath(
  LD_DIR,
  winslash = "/",
  mustWork = TRUE
)

LDCTS <- data.table(
  label = CELL_TYPES,
  ldscore_prefix = paste0(
    LD_PREFIX_ROOT,
    "/GSE238242_",
    CELL_TYPES,
    "."
  )
)

LDCTS_FILE <- file.path(
  LD_DIR,
  "GSE238242_LAA_ATAC.ldcts"
)

fwrite(
  LDCTS,
  LDCTS_FILE,
  sep = "\t",
  col.names = FALSE
)

# ============================================================
# 11. FINAL READINESS QC
# ============================================================

N_EXPECTED_ANNOT <- 12L * 22L

all_annot_exist <- all(
  file.exists(
    unlist(
      lapply(
        CELL_TYPES,
        function(ct) {
          file.path(
            LD_DIR,
            paste0(
              "GSE238242_",
              ct,
              ".",
              1:22,
              ".annot.gz"
            )
          )
        }
      )
    )
  )
)

all_l2_exist <- all(
  file.exists(
    unlist(
      lapply(
        CELL_TYPES,
        function(ct) {
          file.path(
            LD_DIR,
            paste0(
              "GSE238242_",
              ct,
              ".",
              1:22,
              ".l2.ldscore.gz"
            )
          )
        }
      )
    )
  )
)

all_m_exist <- all(
  file.exists(
    unlist(
      lapply(
        CELL_TYPES,
        function(ct) {
          c(
            file.path(
              LD_DIR,
              paste0(
                "GSE238242_",
                ct,
                ".",
                1:22,
                ".l2.M"
              )
            ),
            file.path(
              LD_DIR,
              paste0(
                "GSE238242_",
                ct,
                ".",
                1:22,
                ".l2.M_5_50"
              )
            )
          )
        }
      )
    )
  )
)

READINESS <- data.table(
  check = c(
    "STEP11D_all_pass",
    "12_celltype_BEDs_loaded",
    "1000G_EUR_chr1_22_complete",
    "HapMap3_list_detected",
    "LDSC_python_detected",
    "264_thin_annotation_files_created",
    "264_LDscore_files_created",
    "M_and_M_5_50_files_created",
    "ldcts_file_created"
  ),
  pass = c(
    all(as.logical(D_READY$pass)),
    length(CT_GR) == 12L,
    all(
      vapply(
        1:22,
        function(chr) {
          prefix_triplet_ok(
            REF_PREFIX,
            chr
          )
        },
        logical(1)
      )
    ),
    file.exists(HM3_FILE),
    file.exists(PYTHON_EXE) &&
      file.exists(LDSC_PY),
    all_annot_exist &&
      nrow(ANNOT_QC) == N_EXPECTED_ANNOT,
    if (RUN_LDSC_NOW) all_l2_exist else NA,
    if (RUN_LDSC_NOW) all_m_exist else NA,
    file.exists(LDCTS_FILE)
  )
)

fwrite(
  READINESS,
  file.path(
    QC_DIR,
    "STEP11E_readiness.csv"
  )
)

# ============================================================
# 12. PROVENANCE
# ============================================================

PROVENANCE <- data.table(
  field = c(
    "dataset",
    "tissue",
    "annotation_build",
    "primary_peak_threshold",
    "cell_types",
    "reference_panel",
    "reference_prefix",
    "hapmap3_file",
    "annotation_format",
    "annotation_order_rule",
    "ld_window",
    "print_snps",
    "ldsc_python",
    "ldsc_source",
    "run_ld_scores_now"
  ),
  value = c(
    "GSE238242",
    "Left atrial appendage",
    "GRCh37/hg19",
    ">=10% accessible cells from STEP11C",
    paste(
      CELL_TYPES,
      collapse = ","
    ),
    "1000 Genomes Phase 3 EUR",
    REF_PREFIX,
    HM3_FILE,
    "thin binary annot.gz",
    "exact same SNP count/order as each 1000G.EUR.QC chromosome BIM",
    "1 cM",
    "HapMap3 no-MHC SNPs",
    PYTHON_EXE,
    LDSC_PY,
    as.character(
      RUN_LDSC_NOW
    )
  )
)

fwrite(
  PROVENANCE,
  file.path(
    QC_DIR,
    "STEP11E_method_provenance.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    QC_DIR,
    "STEP11E_sessionInfo.txt"
  )
)

# ============================================================
# 13. CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================================\n",
  "STEP11E COMPLETE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Reference prefix:\n",
  REF_PREFIX,
  "\n\nHapMap3 list:\n",
  HM3_FILE,
  "\n\nCustom LD-score folder:\n",
  LD_DIR,
  "\n\n.ldcts file:\n",
  LDCTS_FILE,
  "\n\n",
  sep = ""
)

cat(
  "Annotation summary:\n"
)

print(
  ANNOT_SUMMARY
)

cat(
  "\nReadiness:\n"
)

print(
  READINESS
)

cat(
  "\nPlease upload these first:\n",
  "  STEP11E_readiness.csv\n",
  "  STEP11E_annotation_summary.csv\n",
  "  STEP11E_annotation_QC_by_celltype_chr.csv\n",
  "  STEP11E_HM3_QC_by_chr.csv\n",
  "  STEP11E_LDscore_job_status.csv\n",
  "  STEP11E_method_provenance.csv\n",
  "  GSE238242_LAA_ATAC.ldcts\n",
  "\nIf a job fails, also upload the corresponding STEP11E_*_chr*.log.\n",
  sep = ""
)
