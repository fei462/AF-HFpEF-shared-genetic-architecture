# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP10B_V2_FINAL_gene_prioritization.R
# See repository README.md for execution order and external dependencies.


# ============================================================
# STEP10B_V2_FINAL
# Evidence-based candidate gene prioritization
# AF / HFpEF / BMI / OSA
#
# IMPORTANT PRINCIPLES
# --------------------
# 1) Each evidence layer contributes at most ONCE per gene.
#    A gene is NOT rewarded simply because it has more probes,
#    more tissues, or more repeated SMR rows.
#
# 2) SMR and HEIDI are treated as a COMPOSITE regulatory layer:
#       eQTL_SMR_HEIDI
#       sQTL_SMR_HEIDI
#    HEIDI is NOT scored independently to avoid double-counting.
#
# 3) n_primary_rows, n_traits, n_tissues, min_p_SMR, p_HEIDI
#    are retained as DESCRIPTIVE SUPPORT only.
#
# 4) This step does NOT rerun GWAS, conjFDR, LAVA, coloc,
#    SuSiE, FUMA, SMR, or HEIDI.
#    It only integrates FINAL confirmed evidence.
#
# 5) Main figures contain NO legend; legends can be written
#    separately in the manuscript/Word file.
#
# Input:
# D:/A/data/STEP10A/STEP10A_FINAL_V3_RESULTS/03_MAIN_TABLES/
#   Table1_priority_genes.csv
#   Table2_regulatory_evidence.csv
#
# Output:
# D:/A/data/STEP10B/STEP10B_V2_FINAL_RESULTS/
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# 0. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ============================================================
# 1. PATHS
# ============================================================

INPUT_DIR <- "D:/A/data/STEP10A/STEP10A_FINAL_V3_RESULTS/03_MAIN_TABLES"

GENE_SCORE_FILE <- file.path(
  INPUT_DIR,
  "Table1_priority_genes.csv"
)

REGULATORY_FILE <- file.path(
  INPUT_DIR,
  "Table2_regulatory_evidence.csv"
)

OUT_ROOT <- "D:/A/data/STEP10B/STEP10B_V2_FINAL_RESULTS"

DIRS <- c(
  "00_QC",
  "01_SOURCE_DATA",
  "02_MAIN_TABLES",
  "03_SUPPLEMENTARY_TABLES",
  "04_MAIN_FIGURES",
  "05_SUPPLEMENTARY_FIGURES",
  "06_LOGS"
)

for (d in DIRS) {
  dir.create(
    file.path(
      OUT_ROOT,
      d
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 2. HARD INPUT QC
# ============================================================

required_files <- c(
  GENE_SCORE_FILE,
  REGULATORY_FILE
)

if (any(!file.exists(required_files))) {
  stop(
    "Missing required STEP10B input(s):\n",
    paste(
      required_files[
        !file.exists(required_files)
      ],
      collapse = "\n"
    )
  )
}

gene_score_raw <- fread(
  GENE_SCORE_FILE,
  showProgress = FALSE
)

reg_raw <- fread(
  REGULATORY_FILE,
  showProgress = FALSE
)

if (!nrow(gene_score_raw)) {
  stop(
    "Table1_priority_genes.csv is empty."
  )
}

if (!nrow(reg_raw)) {
  stop(
    "Table2_regulatory_evidence.csv is empty."
  )
}

# Save exact source snapshots used by STEP10B.
fwrite(
  gene_score_raw,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SOURCE_Table1_priority_genes.csv"
  )
)

fwrite(
  reg_raw,
  file.path(
    OUT_ROOT,
    "01_SOURCE_DATA",
    "SOURCE_Table2_regulatory_evidence.csv"
  )
)

# ============================================================
# 3. HELPERS
# ============================================================

pick_col <- function(
  nms,
  candidates,
  required = TRUE,
  label = "column"
) {

  low <- tolower(nms)

  for (cand in candidates) {

    idx <- which(
      low == tolower(cand)
    )

    if (length(idx)) {
      return(
        nms[idx[1]]
      )
    }
  }

  if (required) {
    stop(
      "Cannot identify ",
      label,
      ". Available columns:\n",
      paste(
        nms,
        collapse = ", "
      )
    )
  }

  NA_character_
}

as_num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

as_logical_safe <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  y %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "PASS"
  )
}

clean_gene <- function(x) {

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  y[
    y %in% c(
      "",
      "NA",
      "N/A",
      "NAN",
      "NULL"
    )
  ] <- NA_character_

  y
}

collapse_unique <- function(x) {

  x <- sort(
    unique(
      x[
        !is.na(x) &
        nzchar(
          trimws(
            as.character(x)
          )
        )
      ]
    )
  )

  if (!length(x)) {
    return("")
  }

  paste(
    x,
    collapse = ";"
  )
}

# ============================================================
# 4. STANDARDIZE BASE GENE EVIDENCE TABLE
# ============================================================

base_gene_col <- pick_col(
  names(gene_score_raw),
  c(
    "gene",
    "Gene",
    "symbol",
    "Symbol"
  ),
  TRUE,
  "gene column in Table1"
)

BASE <- data.table(
  gene = clean_gene(
    gene_score_raw[[base_gene_col]]
  )
)

# Evidence columns previously assembled in STEP10A.
# Regardless of whether they are 0/1 or counts, here they are
# converted to binary PRESENT / ABSENT.
evidence_candidates <- list(
  conjFDR = c(
    "conjFDR",
    "conjfdr"
  ),

  LAVA = c(
    "LAVA",
    "lava"
  ),

  coloc = c(
    "coloc",
    "COLOC"
  ),

  SuSiE = c(
    "SuSiE",
    "susie",
    "SUSIE"
  ),

  FUMA = c(
    "FUMA",
    "fuma"
  )
)

for (nm in names(evidence_candidates)) {

  ccol <- pick_col(
    names(gene_score_raw),
    evidence_candidates[[nm]],
    FALSE
  )

  if (is.na(ccol)) {

    BASE[
      ,
      (nm) := 0L
    ]

  } else {

    vv <- as_num(
      gene_score_raw[[ccol]]
    )

    BASE[
      ,
      (nm) :=
        as.integer(
          is.finite(vv) &
          vv > 0
        )
    ]
  }
}

BASE <- unique(
  BASE[
    !is.na(gene)
  ],
  by = "gene"
)

# ============================================================
# 5. STANDARDIZE REGULATORY TABLE
# ============================================================

reg_gene_col <- pick_col(
  names(reg_raw),
  c(
    "gene",
    "Gene",
    "symbol",
    "Symbol"
  ),
  TRUE,
  "gene column in regulatory table"
)

trait_col <- pick_col(
  names(reg_raw),
  c(
    "trait",
    "Trait"
  ),
  FALSE
)

tissue_col <- pick_col(
  names(reg_raw),
  c(
    "tissue",
    "Tissue"
  ),
  FALSE
)

type_col <- pick_col(
  names(reg_raw),
  c(
    "Evidence_type",
    "evidence_type",
    "type",
    "Type",
    "layer",
    "Layer"
  ),
  FALSE
)

psmr_col <- pick_col(
  names(reg_raw),
  c(
    "p_SMR",
    "P_SMR",
    "p_smr"
  ),
  FALSE
)

pheidi_col <- pick_col(
  names(reg_raw),
  c(
    "p_HEIDI",
    "P_HEIDI",
    "p_heidi"
  ),
  FALSE
)

nheidi_col <- pick_col(
  names(reg_raw),
  c(
    "nsnp_HEIDI",
    "nSNP_HEIDI",
    "nsnp_heidi"
  ),
  FALSE
)

primary_col <- pick_col(
  names(reg_raw),
  c(
    "Primary_pass",
    "primary_pass"
  ),
  FALSE
)

bonf_col <- pick_col(
  names(reg_raw),
  c(
    "SMR_Bonferroni_pass",
    "smr_bonferroni_pass"
  ),
  FALSE
)

heidi_status_col <- pick_col(
  names(reg_raw),
  c(
    "HEIDI_status",
    "heidi_status"
  ),
  FALSE
)

REG <- data.table(
  gene = clean_gene(
    reg_raw[[reg_gene_col]]
  ),

  trait = if (
    !is.na(trait_col)
  ) {
    as.character(
      reg_raw[[trait_col]]
    )
  } else {
    NA_character_
  },

  tissue = if (
    !is.na(tissue_col)
  ) {
    as.character(
      reg_raw[[tissue_col]]
    )
  } else {
    NA_character_
  },

  evidence_type = if (
    !is.na(type_col)
  ) {
    as.character(
      reg_raw[[type_col]]
    )
  } else {
    NA_character_
  },

  p_SMR = if (
    !is.na(psmr_col)
  ) {
    as_num(
      reg_raw[[psmr_col]]
    )
  } else {
    NA_real_
  },

  p_HEIDI = if (
    !is.na(pheidi_col)
  ) {
    as_num(
      reg_raw[[pheidi_col]]
    )
  } else {
    NA_real_
  },

  nsnp_HEIDI = if (
    !is.na(nheidi_col)
  ) {
    as_num(
      reg_raw[[nheidi_col]]
    )
  } else {
    NA_real_
  }
)

# ------------------------------------------------------------
# Primary regulatory evidence definition
# ------------------------------------------------------------
# Preferred:
#   use frozen Primary_pass from STEP9 if available.
#
# Fallback:
#   use frozen SMR_Bonferroni_pass + HEIDI_status if available.
#
# Last fallback:
#   P_HEIDI > 0.01 AND nsnp_HEIDI >= 3
#   but WITHOUT reconstructing a new SMR significance threshold.
#   If the original Bonferroni status is unavailable, regulatory
#   evidence will NOT be promoted to primary solely by p_SMR < 0.05.

if (!is.na(primary_col)) {

  REG[
    ,
    Primary_pass :=
      as_logical_safe(
        reg_raw[[primary_col]]
      )
  ]

} else if (
  !is.na(bonf_col) &&
  !is.na(heidi_status_col)
) {

  bonf_v <- as_logical_safe(
    reg_raw[[bonf_col]]
  )

  hs <- as.character(
    reg_raw[[heidi_status_col]]
  )

  REG[
    ,
    Primary_pass :=
      bonf_v &
      hs ==
        "PASS_NO_HETEROGENEITY"
  ]

} else {

  REG[
    ,
    Primary_pass := FALSE
  ]

  warning(
    paste0(
      "Primary_pass / SMR_Bonferroni_pass was unavailable. ",
      "No new primary SMR significance threshold was invented; ",
      "regulatory rows are therefore kept as descriptive only."
    )
  )
}

# Normalize evidence type.
REG[
  ,
  regulatory_layer :=
    fifelse(
      grepl(
        "sqtl",
        evidence_type,
        ignore.case = TRUE
      ),
      "sQTL",
      fifelse(
        grepl(
          "eqtl",
          evidence_type,
          ignore.case = TRUE
        ),
        "eQTL",
        NA_character_
      )
    )
]

REG <- REG[
  !is.na(gene)
]

# ============================================================
# 6. COLLAPSE REGULATORY EVIDENCE TO ONE ROW PER GENE
# ============================================================

REG_SUM <- REG[
  ,
  {

    prim_idx <- which(
      Primary_pass == TRUE
    )

    eq_idx <- which(
      Primary_pass == TRUE &
      regulatory_layer == "eQTL"
    )

    sq_idx <- which(
      Primary_pass == TRUE &
      regulatory_layer == "sQTL"
    )

    best_idx <- if (
      length(prim_idx) &&
      any(
        is.finite(
          p_SMR[prim_idx]
        )
      )
    ) {

      valid <- prim_idx[
        is.finite(
          p_SMR[prim_idx]
        )
      ]

      valid[
        which.min(
          p_SMR[valid]
        )
      ]

    } else {
      NA_integer_
    }

    list(
      SMR_primary_any =
        as.integer(
          length(prim_idx) > 0
        ),

      eQTL_SMR_HEIDI =
        as.integer(
          length(eq_idx) > 0
        ),

      sQTL_SMR_HEIDI =
        as.integer(
          length(sq_idx) > 0
        ),

      min_primary_p_SMR =
        if (
          length(prim_idx) &&
          any(
            is.finite(
              p_SMR[prim_idx]
            )
          )
        ) {
          min(
            p_SMR[prim_idx],
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      best_primary_p_HEIDI =
        if (
          is.finite(best_idx)
        ) {
          p_HEIDI[best_idx]
        } else {
          NA_real_
        },

      best_primary_nsnp_HEIDI =
        if (
          is.finite(best_idx)
        ) {
          nsnp_HEIDI[best_idx]
        } else {
          NA_real_
        },

      n_primary_rows =
        length(prim_idx),

      n_primary_traits =
        uniqueN(
          trait[
            prim_idx
          ]
        ),

      n_primary_tissues =
        uniqueN(
          tissue[
            prim_idx
          ]
        ),

      primary_traits =
        collapse_unique(
          trait[
            prim_idx
          ]
        ),

      primary_tissues =
        collapse_unique(
          tissue[
            prim_idx
          ]
        )
    )
  },
  by = gene
]

# ============================================================
# 7. MERGE FINAL EVIDENCE
# ============================================================

ALL_GENES <- unique(
  data.table(
    gene = c(
      BASE$gene,
      REG_SUM$gene
    )
  )
)

FINAL <- merge(
  ALL_GENES,
  BASE,
  by = "gene",
  all.x = TRUE
)

FINAL <- merge(
  FINAL,
  REG_SUM,
  by = "gene",
  all.x = TRUE
)

binary_cols <- c(
  "conjFDR",
  "LAVA",
  "coloc",
  "SuSiE",
  "FUMA",
  "SMR_primary_any",
  "eQTL_SMR_HEIDI",
  "sQTL_SMR_HEIDI"
)

for (cc in binary_cols) {

  if (!cc %in% names(FINAL)) {
    FINAL[
      ,
      (cc) := 0L
    ]
  }

  FINAL[
    is.na(
      get(cc)
    ),
    (cc) := 0L
  ]

  FINAL[
    ,
    (cc) :=
      as.integer(
        get(cc) > 0
      )
  ]
}

for (cc in c(
  "n_primary_rows",
  "n_primary_traits",
  "n_primary_tissues"
)) {

  FINAL[
    is.na(
      get(cc)
    ),
    (cc) := 0L
  ]
}

FINAL[
  is.na(primary_traits),
  primary_traits := ""
]

FINAL[
  is.na(primary_tissues),
  primary_tissues := ""
]

# ============================================================
# 8. FINAL SCIENTIFICALLY CONSERVATIVE SCORING
# ============================================================
#
# Binary layer count:
#   conjFDR          1
#   LAVA             1
#   coloc            1
#   SuSiE            1
#   FUMA             1
#   eQTL-SMR+HEIDI   1
#   sQTL-SMR+HEIDI   1
#
# Evidence_count range = 0-7.
#
# Weighted Priority_score reflects inferential hierarchy:
#   conjFDR          2
#   LAVA             2
#   coloc            3
#   SuSiE            3
#   FUMA             1
#   eQTL-SMR+HEIDI   2
#   sQTL-SMR+HEIDI   2
#
# Maximum = 15.
#
# Importantly:
# - repeated tissues/probes DO NOT add points;
# - HEIDI is not scored separately;
# - n_primary_rows is descriptive only.
# ============================================================

FINAL[
  ,
  Evidence_count :=
    conjFDR +
    LAVA +
    coloc +
    SuSiE +
    FUMA +
    eQTL_SMR_HEIDI +
    sQTL_SMR_HEIDI
]

FINAL[
  ,
  Priority_score :=
    2 * conjFDR +
    2 * LAVA +
    3 * coloc +
    3 * SuSiE +
    1 * FUMA +
    2 * eQTL_SMR_HEIDI +
    2 * sQTL_SMR_HEIDI
]

FINAL[
  ,
  locus_level_any :=
    as.integer(
      conjFDR +
      LAVA +
      coloc +
      SuSiE >
        0
    )
]

FINAL[
  ,
  regulatory_any :=
    as.integer(
      eQTL_SMR_HEIDI +
      sQTL_SMR_HEIDI >
        0
    )
]

# Priority classes are evidence-structure labels, not causal claims.
FINAL[
  ,
  Priority_class :=
    fifelse(
      Evidence_count >= 4 &
      locus_level_any == 1 &
      regulatory_any == 1,
      "A_multi_layer",

      fifelse(
        Evidence_count >= 3 &
        locus_level_any == 1 &
        regulatory_any == 1,
        "B_convergent",

        fifelse(
          Evidence_count >= 2,
          "C_supportive",
          "D_single_layer"
        )
      )
    )
]

# Stable ranking.
setorder(
  FINAL,
  -Priority_score,
  -Evidence_count,
  min_primary_p_SMR,
  gene
)

FINAL[
  ,
  Rank :=
    seq_len(.N)
]

setcolorder(
  FINAL,
  c(
    "Rank",
    "gene",
    "Priority_class",
    "Priority_score",
    "Evidence_count",
    "conjFDR",
    "LAVA",
    "coloc",
    "SuSiE",
    "FUMA",
    "eQTL_SMR_HEIDI",
    "sQTL_SMR_HEIDI",
    "SMR_primary_any",
    "min_primary_p_SMR",
    "best_primary_p_HEIDI",
    "best_primary_nsnp_HEIDI",
    "n_primary_rows",
    "n_primary_traits",
    "n_primary_tissues",
    "primary_traits",
    "primary_tissues",
    "locus_level_any",
    "regulatory_any"
  )
)

# ============================================================
# 9. EXPORT MAIN + SUPPLEMENTARY TABLES
# ============================================================

ALL_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_STEP10B_all_gene_evidence_scores.csv"
)

fwrite(
  FINAL,
  ALL_FILE
)

TOP30 <- FINAL[
  Evidence_count >= 2
][
  1:min(
    30L,
    .N
  )
]

if (!nrow(TOP30)) {
  TOP30 <- FINAL[
    1:min(
      30L,
      .N
    )
  ]
}

TOP30_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table_STEP10B_TOP30_priority_genes.csv"
)

fwrite(
  TOP30,
  TOP30_FILE
)

# Higher-confidence table only.
HIGH_CONF <- FINAL[
  Priority_class %in%
    c(
      "A_multi_layer",
      "B_convergent"
    )
]

HIGH_CONF_FILE <- file.path(
  OUT_ROOT,
  "02_MAIN_TABLES",
  "Table_STEP10B_high_confidence_genes.csv"
)

fwrite(
  HIGH_CONF,
  HIGH_CONF_FILE
)

# Regulatory source details for top genes.
TOP_REG <- REG[
  gene %in%
    TOP30$gene &
  Primary_pass == TRUE
][
  order(
    gene,
    p_SMR
  )
]

TOP_REG_FILE <- file.path(
  OUT_ROOT,
  "03_SUPPLEMENTARY_TABLES",
  "TableS_STEP10B_TOP30_primary_SMR_details.csv"
)

fwrite(
  TOP_REG,
  TOP_REG_FILE
)

# ============================================================
# 10. MAIN FIGURE SOURCE DATA
# ============================================================

EVIDENCE_LONG <- TOP30[
  ,
  .(
    Rank,
    gene,
    Priority_score,
    Evidence_count,
    conjFDR,
    LAVA,
    coloc,
    SuSiE,
    FUMA,
    eQTL_SMR_HEIDI,
    sQTL_SMR_HEIDI
  )
] %>%
  pivot_longer(
    cols = c(
      "conjFDR",
      "LAVA",
      "coloc",
      "SuSiE",
      "FUMA",
      "eQTL_SMR_HEIDI",
      "sQTL_SMR_HEIDI"
    ),
    names_to = "Evidence_layer",
    values_to = "Present"
  ) %>%
  as.data.table()

evidence_labels <- c(
  conjFDR = "conjFDR",
  LAVA = "LAVA",
  coloc = "Coloc",
  SuSiE = "SuSiE",
  FUMA = "FUMA",
  eQTL_SMR_HEIDI = "eQTL-SMR",
  sQTL_SMR_HEIDI = "sQTL-SMR"
)

EVIDENCE_LONG[
  ,
  Evidence_label :=
    evidence_labels[
      Evidence_layer
    ]
]

EVIDENCE_LONG[
  ,
  Evidence_label :=
    factor(
      Evidence_label,
      levels = c(
        "conjFDR",
        "LAVA",
        "Coloc",
        "SuSiE",
        "FUMA",
        "eQTL-SMR",
        "sQTL-SMR"
      )
    )
]

gene_levels <- rev(
  TOP30$gene
)

EVIDENCE_LONG[
  ,
  gene :=
    factor(
      gene,
      levels = gene_levels
    )
]

FIG1_SOURCE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_Figure_STEP10B_evidence_matrix.csv"
)

fwrite(
  EVIDENCE_LONG,
  FIG1_SOURCE
)

# ============================================================
# 11. MAIN FIGURE — TOP30 MULTI-LAYER EVIDENCE MATRIX
# ============================================================
#
# Style:
# - white background
# - Arial
# - no internal legend
# - restrained high-impact-journal palette
# - single standalone figure
# ============================================================

p_matrix <- ggplot(
  EVIDENCE_LONG,
  aes(
    x = Evidence_label,
    y = gene
  )
) +
  geom_tile(
    fill = "#F4F4F4",
    colour = "#D9D9D9",
    linewidth = 0.35,
    width = 0.90,
    height = 0.90
  ) +
  geom_point(
    data = EVIDENCE_LONG[
      Present == 1
    ],
    aes(
      x = Evidence_label,
      y = gene
    ),
    shape = 21,
    size = 4.0,
    stroke = 0.45,
    fill = "#5B6C8F",
    colour = "#26344F"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0.5,
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      colour = "black",
      face = "italic"
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      12,
      16,
      12,
      12
    )
  )

FIG1_TIFF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_STEP10B_TOP30_evidence_matrix.tiff"
)

ggsave(
  FIG1_TIFF,
  p_matrix,
  width = 7.1,
  height = 8.4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

FIG1_PDF <- file.path(
  OUT_ROOT,
  "04_MAIN_FIGURES",
  "Figure_STEP10B_TOP30_evidence_matrix.pdf"
)

ggsave(
  FIG1_PDF,
  p_matrix,
  width = 7.1,
  height = 8.4,
  units = "in",
  device = cairo_pdf
)

# ============================================================
# 12. SUPPLEMENTARY FIGURE SOURCE DATA
# ============================================================

SCORE_SOURCE <- TOP30[
  ,
  .(
    Rank,
    gene,
    Priority_score,
    Evidence_count,
    Priority_class
  )
]

SCORE_SOURCE[
  ,
  gene :=
    factor(
      gene,
      levels = rev(
        TOP30$gene
      )
    )
]

FIGS1_SOURCE <- file.path(
  OUT_ROOT,
  "01_SOURCE_DATA",
  "SourceData_FigureS_STEP10B_priority_score.csv"
)

fwrite(
  SCORE_SOURCE,
  FIGS1_SOURCE
)

# ============================================================
# 13. SUPPLEMENTARY FIGURE — PRIORITY SCORE
# ============================================================

p_score <- ggplot(
  SCORE_SOURCE,
  aes(
    x = Priority_score,
    y = gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Priority_score,
      y = gene,
      yend = gene
    ),
    linewidth = 0.55,
    colour = "#B9BEC8"
  ) +
  geom_point(
    size = 3.2,
    shape = 21,
    stroke = 0.45,
    fill = "#7A8BA8",
    colour = "#26344F"
  ) +
  labs(
    x = "Multi-layer evidence score",
    y = NULL
  ) +
  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +
  theme(
    axis.text.y = element_text(
      face = "italic",
      colour = "black",
      size = 9
    ),
    axis.text.x = element_text(
      colour = "black"
    ),
    legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(
      10,
      14,
      10,
      10
    )
  )

FIGS1_TIFF <- file.path(
  OUT_ROOT,
  "05_SUPPLEMENTARY_FIGURES",
  "FigureS_STEP10B_TOP30_priority_score.tiff"
)

ggsave(
  FIGS1_TIFF,
  p_score,
  width = 6.3,
  height = 7.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# 14. QC REPORT
# ============================================================

QC <- data.table(
  item = c(
    "input_gene_score_rows",
    "input_regulatory_rows",
    "unique_base_genes",
    "unique_regulatory_genes",
    "final_unique_genes",
    "genes_with_any_primary_SMR",
    "genes_with_primary_eQTL_SMR_HEIDI",
    "genes_with_primary_sQTL_SMR_HEIDI",
    "genes_priority_A",
    "genes_priority_B",
    "genes_priority_C",
    "genes_priority_D",
    "top30_rows"
  ),

  value = c(
    nrow(gene_score_raw),
    nrow(reg_raw),
    uniqueN(BASE$gene),
    uniqueN(REG$gene),
    nrow(FINAL),
    sum(
      FINAL$SMR_primary_any == 1
    ),
    sum(
      FINAL$eQTL_SMR_HEIDI == 1
    ),
    sum(
      FINAL$sQTL_SMR_HEIDI == 1
    ),
    sum(
      FINAL$Priority_class ==
        "A_multi_layer"
    ),
    sum(
      FINAL$Priority_class ==
        "B_convergent"
    ),
    sum(
      FINAL$Priority_class ==
        "C_supportive"
    ),
    sum(
      FINAL$Priority_class ==
        "D_single_layer"
    ),
    nrow(TOP30)
  )
)

QC_FILE <- file.path(
  OUT_ROOT,
  "00_QC",
  "STEP10B_V2_QC_summary.csv"
)

fwrite(
  QC,
  QC_FILE
)

# Column provenance.
PROVENANCE <- data.table(
  field = c(
    "gene_score_input",
    "regulatory_input",
    "base_gene_column",
    "regulatory_gene_column",
    "primary_pass_column",
    "regulatory_type_column",
    "p_SMR_column",
    "p_HEIDI_column",
    "nsnp_HEIDI_column",
    "score_maximum"
  ),

  value = c(
    GENE_SCORE_FILE,
    REGULATORY_FILE,
    base_gene_col,
    reg_gene_col,
    ifelse(
      is.na(primary_col),
      "",
      primary_col
    ),
    ifelse(
      is.na(type_col),
      "",
      type_col
    ),
    ifelse(
      is.na(psmr_col),
      "",
      psmr_col
    ),
    ifelse(
      is.na(pheidi_col),
      "",
      pheidi_col
    ),
    ifelse(
      is.na(nheidi_col),
      "",
      nheidi_col
    ),
    "15"
  )
)

fwrite(
  PROVENANCE,
  file.path(
    OUT_ROOT,
    "00_QC",
    "STEP10B_V2_column_provenance.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_ROOT,
    "06_LOGS",
    "STEP10B_V2_sessionInfo.txt"
  )
)

# ============================================================
# 15. FINAL CONSOLE REPORT
# ============================================================

cat(
  "\n====================================================\n",
  "STEP10B_V2 FINAL COMPLETED\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nQC summary:\n"
)

print(
  QC
)

cat(
  "\nTop 30 priority genes:\n"
)

print(
  TOP30[
    ,
    .(
      Rank,
      gene,
      Priority_class,
      Priority_score,
      Evidence_count,
      conjFDR,
      LAVA,
      coloc,
      SuSiE,
      FUMA,
      eQTL_SMR_HEIDI,
      sQTL_SMR_HEIDI,
      min_primary_p_SMR,
      n_primary_traits,
      n_primary_tissues
    )
  ]
)

cat(
  "\nIMPORTANT:\n",
  "- SMR and HEIDI are NOT counted repeatedly by probe/tissue.\n",
  "- HEIDI is NOT an independent score layer.\n",
  "- Tissue count is descriptive only.\n",
  "- Priority classes are evidence-structure labels, not causal proof.\n",
  sep = ""
)

cat(
  "\nUPLOAD THESE FILES:\n",
  "1) ", QC_FILE, "\n",
  "2) ", TOP30_FILE, "\n",
  "3) ", HIGH_CONF_FILE, "\n",
  "4) ", ALL_FILE, "\n",
  "5) ", FIG1_SOURCE, "\n",
  "6) ", FIG1_TIFF, "\n",
  sep = ""
)

cat(
  "====================================================\n"
)
