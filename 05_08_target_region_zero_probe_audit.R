# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP9C1_V3_ZERO_PROBE_SAFE_4_REGION_AUDIT.R
# See repository README.md for execution order and external dependencies.


# ============================================================
# STEP9C1_V3
# ZERO-PROBE-SAFE 4-REGION SMR/HEIDI AUDIT
# AF / HFpEF / BMI / OSA
#
# IMPORTANT:
# - DOES NOT rerun SMR/HEIDI
# - reads the completed STEP9B master table
# - keeps all 4 pre-specified loci even if a locus has zero GTEx probes
# - a zero-probe locus is classified as NOT_EVALUABLE_NO_GTEX_PROBE
#   rather than treated as an error or a biological negative result
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
library(data.table)

# ============================================================
# 0. PATHS
# ============================================================

ROOT <- "D:/A/data/STEP9_SMR"

MASTER_FILE <- file.path(
  ROOT,
  "05_tables",
  "STEP9B",
  "STEP9B_all_32_SMR_HEIDI_results.csv.gz"
)

OUT_DIR <- file.path(
  ROOT,
  "05_tables",
  "STEP9C1_V3"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(MASTER_FILE)) {
  stop("Missing STEP9B master file: ", MASTER_FILE)
}

# ============================================================
# 1. LOAD STEP9B MASTER
# ============================================================

ALL <- fread(
  MASTER_FILE,
  showProgress = TRUE
)

required_cols <- c(
  "trait",
  "tissue",
  "ProbeID",
  "Probe_Chr",
  "Gene",
  "Probe_bp",
  "p_SMR",
  "p_HEIDI",
  "nsnp_HEIDI",
  "SMR_Bonferroni_pass",
  "HEIDI_status",
  "Primary_pass",
  "gene_key"
)

missing_cols <- setdiff(
  required_cols,
  names(ALL)
)

if (length(missing_cols)) {
  stop(
    "STEP9B master table lacks columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

ALL[
  ,
  Probe_Chr_num := suppressWarnings(
    as.integer(
      gsub(
        "^chr",
        "",
        as.character(Probe_Chr),
        ignore.case = TRUE
      )
    )
  )
]

# ============================================================
# 2. FOUR PRE-SPECIFIED TARGET LOCI
#    hg19 / GRCh37
# ============================================================

REGIONS <- data.table(
  region = c(
    "chr11_recurrent",
    "chr16_locus2135_FTO",
    "chr3_BMI_centered",
    "chr2_AF_HFpEF"
  ),

  chr = c(
    11L,
    16L,
    3L,
    2L
  ),

  start = c(
    16024798L,
    53393883L,
    185170000L,
    200401427L
  ),

  stop = c(
    17024798L,
    54866095L,
    186170000L,
    201401427L
  ),

  definition = c(
    "rs11529589-centered +/-500 kb recurrent chr11 region",
    "Exact LAVA locus 2135 / FTO-region",
    "Fixed 1 Mb recurrent chr3 BMI-centered region",
    "rs203768-centered +/-500 kb AF-HFpEF region"
  )
)

REGION_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_frozen_target_regions.csv"
)

fwrite(
  REGIONS,
  REGION_FILE
)

stopifnot(
  nrow(REGIONS) == 4L
)

# ============================================================
# 3. EXTRACT TARGET-REGION ROWS
# ============================================================

TARGET_LIST <- vector(
  "list",
  nrow(REGIONS)
)

for (i in seq_len(nrow(REGIONS))) {

  rr <- REGIONS[i]

  z <- ALL[
    Probe_Chr_num == rr$chr &
    is.finite(Probe_bp) &
    Probe_bp >= rr$start &
    Probe_bp <= rr$stop
  ]

  if (nrow(z)) {

    z[
      ,
      `:=`(
        target_region = rr$region,
        target_start = rr$start,
        target_stop = rr$stop,
        target_definition = rr$definition
      )
    ]
  }

  TARGET_LIST[[i]] <- z
}

TARGET <- rbindlist(
  TARGET_LIST,
  fill = TRUE
)

if (nrow(TARGET)) {
  setorder(
    TARGET,
    target_region,
    trait,
    p_SMR
  )
}

TARGET_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_all_target_region_SMR_results.csv"
)

fwrite(
  TARGET,
  TARGET_FILE
)

# ============================================================
# 4. REGION EXTRACTION QC
# ============================================================

if (nrow(TARGET)) {

  REGION_OBS <- TARGET[
    ,
    .(
      n_rows = .N,
      n_probes = uniqueN(ProbeID),
      n_traits = uniqueN(trait),
      traits = paste(
        sort(unique(trait)),
        collapse = ";"
      )
    ),
    by = target_region
  ]

} else {

  REGION_OBS <- data.table(
    target_region = character(),
    n_rows = integer(),
    n_probes = integer(),
    n_traits = integer(),
    traits = character()
  )
}

REGION_QC <- merge(
  REGIONS[
    ,
    .(
      target_region = region
    )
  ],
  REGION_OBS,
  by = "target_region",
  all.x = TRUE
)

for (cc in c("n_rows", "n_probes", "n_traits")) {
  REGION_QC[
    is.na(get(cc)),
    (cc) := 0L
  ]
}

REGION_QC[
  is.na(traits),
  traits := ""
]

REGION_QC[
  ,
  evaluation_status := ifelse(
    n_rows > 0L,
    "EVALUABLE",
    "NOT_EVALUABLE_NO_GTEX_PROBE"
  )
]

REGION_QC[
  ,
  region_order := match(
    target_region,
    REGIONS$region
  )
]

setorder(
  REGION_QC,
  region_order
)

REGION_QC[
  ,
  region_order := NULL
]

REGION_QC_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_target_region_extraction_QC.csv"
)

fwrite(
  REGION_QC,
  REGION_QC_FILE
)

# ============================================================
# 5. TARGET-REGION SUMMARY BY TRAIT
# ============================================================

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

if (nrow(TARGET)) {

  TARGET_SUMMARY_OBS <- TARGET[
    ,
    {
      idx <- which(
        is.finite(p_SMR)
      )

      best_idx <- if (length(idx)) {
        idx[
          which.min(
            p_SMR[idx]
          )
        ]
      } else {
        NA_integer_
      }

      list(
        n_tested_rows = .N,

        n_tested_probes =
          uniqueN(
            ProbeID
          ),

        min_p_SMR =
          safe_min(
            p_SMR
          ),

        best_gene =
          if (!is.na(best_idx)) {
            as.character(
              Gene[best_idx]
            )
          } else {
            NA_character_
          },

        best_probe =
          if (!is.na(best_idx)) {
            as.character(
              ProbeID[best_idx]
            )
          } else {
            NA_character_
          },

        best_tissue =
          if (!is.na(best_idx)) {
            as.character(
              tissue[best_idx]
            )
          } else {
            NA_character_
          },

        n_Bonferroni =
          sum(
            SMR_Bonferroni_pass,
            na.rm = TRUE
          ),

        n_primary_pass =
          sum(
            Primary_pass,
            na.rm = TRUE
          ),

        n_HEIDI_heterogeneity =
          sum(
            SMR_Bonferroni_pass &
            HEIDI_status == "FAIL_HETEROGENEITY",
            na.rm = TRUE
          ),

        n_HEIDI_not_evaluable =
          sum(
            SMR_Bonferroni_pass &
            HEIDI_status == "NOT_EVALUABLE",
            na.rm = TRUE
          ),

        primary_genes =
          paste(
            sort(
              unique(
                gene_key[
                  Primary_pass == TRUE
                ]
              )
            ),
            collapse = ";"
          )
      )
    },
    by = .(
      target_region,
      trait
    )
  ]

} else {

  TARGET_SUMMARY_OBS <- data.table(
    target_region = character(),
    trait = character()
  )
}

TRAITS <- c(
  "AF",
  "HFpEF",
  "BMI",
  "OSA"
)

GRID <- CJ(
  target_region = REGIONS$region,
  trait = TRAITS,
  unique = TRUE
)

TARGET_SUMMARY <- merge(
  GRID,
  TARGET_SUMMARY_OBS,
  by = c(
    "target_region",
    "trait"
  ),
  all.x = TRUE
)

count_cols <- c(
  "n_tested_rows",
  "n_tested_probes",
  "n_Bonferroni",
  "n_primary_pass",
  "n_HEIDI_heterogeneity",
  "n_HEIDI_not_evaluable"
)

for (cc in count_cols) {

  if (cc %in% names(TARGET_SUMMARY)) {

    TARGET_SUMMARY[
      is.na(get(cc)),
      (cc) := 0
    ]
  }
}

if ("primary_genes" %in% names(TARGET_SUMMARY)) {
  TARGET_SUMMARY[
    is.na(primary_genes),
    primary_genes := ""
  ]
}

TARGET_SUMMARY[
  ,
  evaluation_status := ifelse(
    n_tested_rows > 0,
    "EVALUABLE",
    "NOT_EVALUABLE_NO_GTEX_PROBE"
  )
]

TARGET_SUMMARY[
  ,
  region_order := match(
    target_region,
    REGIONS$region
  )
]

TARGET_SUMMARY[
  ,
  trait_order := match(
    trait,
    TRAITS
  )
]

setorder(
  TARGET_SUMMARY,
  region_order,
  trait_order
)

TARGET_SUMMARY[
  ,
  c(
    "region_order",
    "trait_order"
  ) := NULL
]

TARGET_SUMMARY_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_target_region_summary.csv"
)

fwrite(
  TARGET_SUMMARY,
  TARGET_SUMMARY_FILE
)

# ============================================================
# 6. TARGET-REGION STRICT PRIMARY HITS
# ============================================================

TARGET_PRIMARY <- TARGET[
  Primary_pass == TRUE
][
  order(
    target_region,
    trait,
    p_SMR
  )
]

TARGET_PRIMARY_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_target_region_primary_hits.csv"
)

fwrite(
  TARGET_PRIMARY,
  TARGET_PRIMARY_FILE
)

# ============================================================
# 7. STRICT CROSS-TRAIT SMR GENE OVERLAP
# ============================================================

PRIMARY <- ALL[
  Primary_pass == TRUE
]

PAIR_LIST <- list()
counter <- 0L

for (i in seq_len(length(TRAITS) - 1L)) {

  for (j in (i + 1L):length(TRAITS)) {

    a <- TRAITS[i]
    b <- TRAITS[j]

    ga <- unique(
      PRIMARY[
        trait == a,
        gene_key
      ]
    )

    gb <- unique(
      PRIMARY[
        trait == b,
        gene_key
      ]
    )

    shared <- intersect(
      ga,
      gb
    )

    if (!length(shared)) next

    for (g in shared) {

      aa <- PRIMARY[
        trait == a &
        gene_key == g
      ][
        order(p_SMR)
      ]

      bb <- PRIMARY[
        trait == b &
        gene_key == g
      ][
        order(p_SMR)
      ]

      counter <- counter + 1L

      PAIR_LIST[[counter]] <- data.table(
        trait1 = a,
        trait2 = b,
        gene_key = g,

        trait1_best_tissue =
          aa$tissue[1],

        trait1_min_p_SMR =
          aa$p_SMR[1],

        trait1_p_HEIDI =
          aa$p_HEIDI[1],

        trait2_best_tissue =
          bb$tissue[1],

        trait2_min_p_SMR =
          bb$p_SMR[1],

        trait2_p_HEIDI =
          bb$p_HEIDI[1],

        chromosome =
          aa$Probe_Chr_num[1],

        probe_bp =
          aa$Probe_bp[1]
      )
    }
  }
}

SHARED <- if (length(PAIR_LIST)) {

  rbindlist(
    PAIR_LIST,
    fill = TRUE
  )

} else {

  data.table(
    trait1 = character(),
    trait2 = character(),
    gene_key = character(),
    trait1_best_tissue = character(),
    trait1_min_p_SMR = numeric(),
    trait1_p_HEIDI = numeric(),
    trait2_best_tissue = character(),
    trait2_min_p_SMR = numeric(),
    trait2_p_HEIDI = numeric(),
    chromosome = integer(),
    probe_bp = numeric()
  )
}

SHARED_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_strict_cross_trait_SMR_gene_overlap.csv"
)

fwrite(
  SHARED,
  SHARED_FILE
)

# ============================================================
# 8. HEIDI-NOT-EVALUABLE AUDIT AGAINST ALL 4 LOCI
# ============================================================

NE <- ALL[
  SMR_Bonferroni_pass == TRUE &
  HEIDI_status == "NOT_EVALUABLE"
]

NE[
  ,
  `:=`(
    in_focal_region = FALSE,
    focal_region_name = NA_character_
  )
]

if (nrow(NE)) {

  for (i in seq_len(nrow(REGIONS))) {

    rr <- REGIONS[i]

    idx <- which(
      NE$Probe_Chr_num == rr$chr &
      is.finite(NE$Probe_bp) &
      NE$Probe_bp >= rr$start &
      NE$Probe_bp <= rr$stop
    )

    if (length(idx)) {

      NE$in_focal_region[idx] <- TRUE
      NE$focal_region_name[idx] <- rr$region
    }
  }
}

NE_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_HEIDI_not_evaluable_focal_region_audit.csv"
)

fwrite(
  NE,
  NE_FILE
)

# ============================================================
# 9. FINAL SUMMARY
# ============================================================

SUMMARY <- data.table(
  item = c(
    "master_rows",
    "primary_rows",
    "primary_AF_rows",
    "primary_HFpEF_rows",
    "primary_BMI_rows",
    "primary_OSA_rows",
    "strict_cross_trait_gene_pairs",
    "target_regions_frozen",
    "target_regions_evaluable",
    "target_regions_no_GTEx_probe",
    "target_region_primary_rows",
    "HEIDI_not_evaluable_rows",
    "HEIDI_not_evaluable_in_any_focal_region"
  ),

  value = c(
    nrow(ALL),
    nrow(PRIMARY),
    sum(PRIMARY$trait == "AF"),
    sum(PRIMARY$trait == "HFpEF"),
    sum(PRIMARY$trait == "BMI"),
    sum(PRIMARY$trait == "OSA"),
    nrow(SHARED),
    nrow(REGIONS),
    sum(REGION_QC$evaluation_status == "EVALUABLE"),
    sum(REGION_QC$evaluation_status == "NOT_EVALUABLE_NO_GTEX_PROBE"),
    nrow(TARGET_PRIMARY),
    nrow(NE),
    sum(
      NE$in_focal_region,
      na.rm = TRUE
    )
  )
)

SUMMARY_FILE <- file.path(
  OUT_DIR,
  "STEP9C1_V3_summary.csv"
)

fwrite(
  SUMMARY,
  SUMMARY_FILE
)

# ============================================================
# 10. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP9C1_V3 COMPLETE — ZERO-PROBE-SAFE 4-REGION AUDIT\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nRegion QC:\n"
)
print(
  REGION_QC
)

cat(
  "\nTarget-region summary:\n"
)
print(
  TARGET_SUMMARY
)

cat(
  "\nStrict cross-trait SMR genes:\n"
)
print(
  SHARED
)

cat(
  "\nOverall summary:\n"
)
print(
  SUMMARY
)

cat(
  "\nUPLOAD THESE 6 FILES:\n",
  "1) ", SUMMARY_FILE, "\n",
  "2) ", REGION_QC_FILE, "\n",
  "3) ", TARGET_SUMMARY_FILE, "\n",
  "4) ", TARGET_PRIMARY_FILE, "\n",
  "5) ", SHARED_FILE, "\n",
  "6) ", NE_FILE, "\n",
  sep = ""
)
