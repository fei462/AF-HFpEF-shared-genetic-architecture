# CODE RELEASE v1.0
# Curated final script. Original working filename: STEP5A_pleioFDR_prep_condQQ_V2.R
# See repository README.md for execution order and external dependencies.

# ============================================================
# STEP5A — pleioFDR preparation + conditional Q-Q
# Project: AF–HFpEF–BMI–OSA shared genetics
# Purpose:
#   1) Prepare standardized pleioFDR inputs from the final full-GWAS standardized files
#   2) Generate bidirectional conditional Q-Q plots for all prespecified pairs
#   3) Create QC/readiness tables and official pleioFDR config templates
#
# IMPORTANT:
# - This step is independent of MiXeR; MiXeR is not part of the final manuscript.
# - The conditional Q-Q plots generated here are screening/pre-analysis plots.
# - Final manuscript-grade cond/conjFDR values and Q-Q plots should come from the
#   official pleioFDR pipeline after matching to the 1000G EUR reference and
#   repeated LD random-pruning.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# --------------------------- USER PATHS ---------------------------
ROOT <- "D:/A/data"
FULLGWAS_IN <- file.path(ROOT, "STEP5_PLEIOFDR", "00_fullGWAS_standardized", "01_inputs")
OUT <- file.path(ROOT, "STEP5_PLEIOFDR")

DIR_INPUT <- file.path(OUT, "01_inputs_prepared")
DIR_MAT <- file.path(OUT, "02_mat")
DIR_QQ <- file.path(OUT, "03_conditional_QQ")
DIR_CFG <- file.path(OUT, "04_config_templates")
DIR_RES <- file.path(OUT, "05_results")
DIR_QC <- file.path(OUT, "00_qc")
DIR_TMP <- file.path(OUT, "tmp")

dirs <- c(OUT, DIR_INPUT, DIR_MAT, DIR_QQ, DIR_CFG, DIR_RES, DIR_QC, DIR_TMP)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# --------------------------- PACKAGES -----------------------------
pkgs <- c("data.table", "ggplot2", "DBI", "duckdb")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(data.table)
library(ggplot2)
library(DBI)
library(duckdb)
library(grid)

# Memory-conscious settings for large full-GWAS files
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(MKL_NUM_THREADS = "2")

# --------------------------- INPUT FILES --------------------------
trait_files <- c(
  AF    = file.path(FULLGWAS_IN, "AF_fullGWAS.sumstats.gz"),
  HFpEF = file.path(FULLGWAS_IN, "HFpEF_fullGWAS.sumstats.gz"),
  BMI   = file.path(FULLGWAS_IN, "BMI_fullGWAS.sumstats.gz"),
  OSA   = file.path(FULLGWAS_IN, "OSA_fullGWAS.sumstats.gz")
)

missing_files <- trait_files[!file.exists(trait_files)]
if (length(missing_files)) {
  stop("Missing standardized full-GWAS input file(s):\n", paste(missing_files, collapse = "\n"))
}

pairs <- list(
  c("AF", "HFpEF"),
  c("AF", "BMI"),
  c("AF", "OSA"),
  c("HFpEF", "BMI"),
  c("HFpEF", "OSA"),
  c("BMI", "OSA")
)

# --------------------------- HELPERS ------------------------------

is_ambiguous <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

z_to_logp <- function(z) {
  lp_ln <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  -lp_ln / log(10)
}

logp_to_safe_p <- function(logp10) {
  max_logp <- -log10(.Machine$double.xmin)
  10^(-pmin(logp10, max_logp))
}

standardize_trait <- function(trait, infile, outfile) {
  message("\n[Prepare] ", trait, " <- ", infile)

  hdr <- names(fread(infile, nrows = 0))
  need <- c("SNP", "A1", "A2", "N", "Z")
  miss <- setdiff(need, hdr)
  if (length(miss)) stop(trait, ": missing columns: ", paste(miss, collapse = ", "))

  dt <- fread(infile, select = need, showProgress = TRUE)
  raw_n <- nrow(dt)

  dt[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    N = as.numeric(N),
    Z = as.numeric(Z)
  )]

  dt <- dt[
    !is.na(SNP) & nzchar(SNP) &
      A1 %chin% c("A", "C", "G", "T") &
      A2 %chin% c("A", "C", "G", "T") &
      A1 != A2 &
      is.finite(N) & N > 0 &
      is.finite(Z)
  ]
  n_valid <- nrow(dt)

  dt <- dt[!is_ambiguous(A1, A2)]
  n_nonamb <- nrow(dt)

  dt[, LOGP := z_to_logp(Z)]
  dt[, P := logp_to_safe_p(LOGP)]

  dup_n <- sum(duplicated(dt$SNP))
  if (dup_n > 0) {
    setorderv(dt, c("SNP", "LOGP"), c(1L, -1L))
    dt <- dt[!duplicated(SNP)]
  }

  setcolorder(dt, c("SNP", "A1", "A2", "P", "Z", "N", "LOGP"))

  fwrite(
    dt,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  fwrite(head(dt, 20),
         file.path(DIR_QC, paste0(trait, "_pleioFDR_input_head20.tsv")),
         sep = "\t", quote = FALSE, na = "NA")

  q <- data.table(
    trait = trait,
    raw_rows = raw_n,
    valid_ACGT_rows = n_valid,
    non_ambiguous_rows_before_dedup = n_nonamb,
    duplicate_SNP_rows_removed = dup_n,
    final_rows = nrow(dt),
    min_N = suppressWarnings(min(dt$N, na.rm = TRUE)),
    median_N = suppressWarnings(median(dt$N, na.rm = TRUE)),
    max_N = suppressWarnings(max(dt$N, na.rm = TRUE)),
    min_P_stored = suppressWarnings(min(dt$P, na.rm = TRUE)),
    max_abs_Z = suppressWarnings(max(abs(dt$Z), na.rm = TRUE)),
    output_file = outfile
  )

  rm(dt); gc()
  q
}

make_condqq_data <- function(lp_target, lp_cond,
                             target_name, cond_name,
                             cond_breaks = c(0, 1, 2, 3),
                             y_grid = seq(0, 7.3, by = 0.05)) {

  stopifnot(length(lp_target) == length(lp_cond))
  good <- is.finite(lp_target) & is.finite(lp_cond)
  lp_target <- lp_target[good]
  lp_cond <- lp_cond[good]

  out <- vector("list", length(cond_breaks))
  counts <- vector("list", length(cond_breaks))

  for (j in seq_along(cond_breaks)) {
    cb <- cond_breaks[j]
    keep <- lp_cond >= cb
    x <- lp_target[keep]
    n <- length(x)

    lab <- if (cb == 0) {
      "All SNPs"
    } else {
      paste0("P(", cond_name, ") <= 1e-", cb)
    }

    if (n < 100L) {
      warning(target_name, " | ", cond_name, ": only ", n,
              " variants in conditioning stratum ", lab)
    }

    if (n == 0L) {
      out[[j]] <- data.table(
        target = target_name, conditioned_on = cond_name,
        stratum = lab, cond_log10_threshold = cb,
        nominal_logp = y_grid, empirical_logq = NA_real_, n = 0L
      )
      counts[[j]] <- data.table(
        target = target_name, conditioned_on = cond_name,
        stratum = lab, cond_log10_threshold = cb, n = 0L
      )
      next
    }

    sx <- sort(x)
    idx_le <- findInterval(y_grid, sx)
    tail_n <- n - idx_le
    tail_n[y_grid == 0] <- n
    q_emp <- pmax(tail_n / n, 1 / n)
    empirical_logq <- -log10(q_emp)

    out[[j]] <- data.table(
      target = target_name,
      conditioned_on = cond_name,
      stratum = lab,
      cond_log10_threshold = cb,
      nominal_logp = y_grid,
      empirical_logq = empirical_logq,
      n = n
    )

    counts[[j]] <- data.table(
      target = target_name,
      conditioned_on = cond_name,
      stratum = lab,
      cond_log10_threshold = cb,
      n = n
    )
  }

  list(data = rbindlist(out), counts = rbindlist(counts))
}

save_condqq_plot <- function(qqdt, target_name, cond_name, stem) {
  labs_cond <- c(
    "All SNPs",
    paste0("P(", cond_name, ") <= 1e-1"),
    paste0("P(", cond_name, ") <= 1e-2"),
    paste0("P(", cond_name, ") <= 1e-3")
  )
  cols <- setNames(
    c("#5B5B5B", "#9ECAE1", "#4292C6", "#CB181D"),
    labs_cond
  )

  qqdt[, stratum := factor(stratum, levels = labs_cond)]

  p <- ggplot(qqdt[is.finite(empirical_logq)],
              aes(x = empirical_logq, y = nominal_logp,
                  color = stratum, group = stratum)) +
    geom_abline(intercept = 0, slope = 1,
                linewidth = 0.55, linetype = "dashed", color = "black") +
    geom_line(linewidth = 0.95, lineend = "round") +
    scale_color_manual(values = cols, drop = FALSE) +
    coord_cartesian(xlim = c(0, 7.3), ylim = c(0, 7.3), expand = FALSE) +
    scale_x_continuous(breaks = 0:7) +
    scale_y_continuous(breaks = 0:7) +
    labs(
      title = paste0(target_name, " | ", cond_name),
      x = bquote("Empirical " * -log[10](q)),
      y = bquote("Nominal " * -log[10](p))
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, face = "plain", hjust = 0),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      plot.margin = margin(8, 10, 8, 8)
    )

  ggsave(paste0(stem, ".pdf"), p, width = 4.8, height = 4.8,
         device = cairo_pdf)
  ggsave(paste0(stem, ".tiff"), p, width = 4.8, height = 4.8,
         dpi = 600, compression = "lzw", bg = "white")
}

save_legend <- function(cond_name, outfile) {
  labs <- c(
    "All SNPs",
    paste0("P(", cond_name, ") <= 0.1"),
    paste0("P(", cond_name, ") <= 0.01"),
    paste0("P(", cond_name, ") <= 0.001"),
    "Expected"
  )
  cols <- c("#5B5B5B", "#9ECAE1", "#4292C6", "#CB181D", "#000000")
  ltys <- c(1, 1, 1, 1, 2)

  tiff(outfile, width = 1600, height = 700, res = 300, compression = "lzw")
  grid.newpage()
  ys <- seq(0.83, 0.17, length.out = length(labs))
  for (i in seq_along(labs)) {
    grid.lines(x = unit(c(0.08, 0.24), "npc"),
               y = unit(c(ys[i], ys[i]), "npc"),
               gp = gpar(col = cols[i], lwd = 3, lty = ltys[i]))
    grid.text(labs[i],
              x = unit(0.29, "npc"), y = unit(ys[i], "npc"),
              just = "left", gp = gpar(fontfamily = "Arial", fontsize = 13))
  }
  dev.off()
}

# ---------------------- 1) PREPARE TRAITS ------------------------

prepared_files <- setNames(
  file.path(DIR_INPUT, paste0(names(trait_files), "_pleioFDR_input.tsv.gz")),
  names(trait_files)
)

qc_list <- vector("list", length(trait_files))
names(qc_list) <- names(trait_files)

for (tr in names(trait_files)) {
  qc_list[[tr]] <- standardize_trait(
    trait = tr,
    infile = trait_files[[tr]],
    outfile = prepared_files[[tr]]
  )
}

qc_traits <- rbindlist(qc_list, fill = TRUE)
fwrite(qc_traits,
       file.path(DIR_QC, "STEP5A_trait_preparation_QC.csv"))

# ---------------- 2) CONDITIONAL Q-Q, BOTH DIRECTIONS ------------

con <- dbConnect(duckdb::duckdb(), dbdir = file.path(DIR_TMP, "condqq.duckdb"))
dbExecute(con, "PRAGMA threads=2")
dbExecute(con, "PRAGMA memory_limit='4GB'")

qq_counts_all <- list()
qq_shift_all <- list()
pair_audit <- list()

sql_path <- function(x) gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = TRUE))

for (pr in pairs) {
  t1 <- pr[1]
  t2 <- pr[2]
  message("\n[Conditional Q-Q] ", t1, " <-> ", t2)

  f1 <- sql_path(prepared_files[[t1]])
  f2 <- sql_path(prepared_files[[t2]])

  sql <- sprintf("
    SELECT a.LOGP AS LP1, b.LOGP AS LP2
    FROM read_csv_auto('%s', delim='\\t', header=true, compression='gzip') a
    INNER JOIN read_csv_auto('%s', delim='\\t', header=true, compression='gzip') b
    ON a.SNP = b.SNP
    WHERE isfinite(a.LOGP) AND isfinite(b.LOGP)
  ", f1, f2)

  dat <- as.data.table(dbGetQuery(con, sql))
  n_overlap <- nrow(dat)

  pair_audit[[length(pair_audit) + 1L]] <- data.table(
    trait1 = t1, trait2 = t2, common_SNPs_for_screening_QQ = n_overlap
  )

  q1 <- make_condqq_data(dat$LP1, dat$LP2, t1, t2)
  stem1 <- file.path(DIR_QQ, paste0("CondQQ_", t1, "_given_", t2))
  fwrite(q1$data, paste0(stem1, "_source.csv"))
  save_condqq_plot(q1$data, t1, t2, stem1)
  save_legend(t2, paste0(stem1, "_legend.tiff"))

  q2 <- make_condqq_data(dat$LP2, dat$LP1, t2, t1)
  stem2 <- file.path(DIR_QQ, paste0("CondQQ_", t2, "_given_", t1))
  fwrite(q2$data, paste0(stem2, "_source.csv"))
  save_condqq_plot(q2$data, t2, t1, stem2)
  save_legend(t1, paste0(stem2, "_legend.tiff"))

  qq_counts_all[[length(qq_counts_all) + 1L]] <- q1$counts
  qq_counts_all[[length(qq_counts_all) + 1L]] <- q2$counts

  for (qqx in list(q1$data, q2$data)) {
    tmp <- qqx[abs(nominal_logp - 3) < 1e-12,
               .(target, conditioned_on, stratum,
                 cond_log10_threshold, n,
                 empirical_logq,
                 left_shift_at_nominal_1e3 = 3 - empirical_logq)]
    qq_shift_all[[length(qq_shift_all) + 1L]] <- tmp
  }

  rm(dat, q1, q2); gc()
}

fwrite(rbindlist(pair_audit),
       file.path(DIR_QC, "STEP5A_pair_overlap_QC.csv"))
fwrite(rbindlist(qq_counts_all, fill = TRUE),
       file.path(DIR_QC, "STEP5A_condQQ_stratum_counts.csv"))
fwrite(rbindlist(qq_shift_all, fill = TRUE),
       file.path(DIR_QC, "STEP5A_condQQ_left_shift_diagnostic.csv"))

dbDisconnect(con, shutdown = TRUE)
con <- NULL

# ---------------- 3) OFFICIAL pleioFDR RESOURCE PREFLIGHT --------

reference_candidates <- c(
  file.path(OUT, "reference", "ref9545380_1kgPhase3eur_LDr2p1.mat"),
  file.path(ROOT, "ref9545380_1kgPhase3eur_LDr2p1.mat"),
  file.path(ROOT, "pleiofdr", "ref9545380_1kgPhase3eur_LDr2p1.mat")
)

refinfo_candidates <- c(
  file.path(OUT, "reference", "9545380.ref"),
  file.path(ROOT, "9545380.ref"),
  file.path(ROOT, "pleiofdr", "9545380.ref")
)

repo_candidates <- c(
  file.path(OUT, "software", "pleiofdr"),
  file.path(ROOT, "pleiofdr")
)

first_existing <- function(x) {
  y <- x[file.exists(x) | dir.exists(x)]
  if (length(y)) y[1] else NA_character_
}

ref_mat <- first_existing(reference_candidates)
ref_txt <- first_existing(refinfo_candidates)
repo_dir <- first_existing(repo_candidates)
matlab_bin <- Sys.which("matlab")
git_bin <- Sys.which("git")

readiness <- data.table(
  item = c(
    "MATLAB executable",
    "git executable",
    "pleioFDR repository",
    "ref9545380_1kgPhase3eur_LDr2p1.mat",
    "9545380.ref",
    "AF prepared input",
    "HFpEF prepared input",
    "BMI prepared input",
    "OSA prepared input"
  ),
  found = c(
    nzchar(matlab_bin),
    nzchar(git_bin),
    !is.na(repo_dir),
    !is.na(ref_mat),
    !is.na(ref_txt),
    file.exists(prepared_files["AF"]),
    file.exists(prepared_files["HFpEF"]),
    file.exists(prepared_files["BMI"]),
    file.exists(prepared_files["OSA"])
  ),
  path = c(
    ifelse(nzchar(matlab_bin), matlab_bin, NA),
    ifelse(nzchar(git_bin), git_bin, NA),
    repo_dir,
    ref_mat,
    ref_txt,
    prepared_files["AF"],
    prepared_files["HFpEF"],
    prepared_files["BMI"],
    prepared_files["OSA"]
  )
)

fwrite(readiness,
       file.path(DIR_QC, "STEP5A_pleioFDR_readiness.csv"))

# ---------------- 4) CONFIG TEMPLATES FOR SIX PAIRS --------------

mat_names <- c(
  AF = "AF_pleioFDR.mat",
  HFpEF = "HFpEF_pleioFDR.mat",
  BMI = "BMI_pleioFDR.mat",
  OSA = "OSA_pleioFDR.mat"
)

ref_mat_cfg <- if (!is.na(ref_mat)) gsub("\\\\", "/", ref_mat) else
  gsub("\\\\", "/", file.path(OUT, "reference", "ref9545380_1kgPhase3eur_LDr2p1.mat"))

mat_dir_cfg <- gsub("\\\\", "/", DIR_MAT)
res_dir_cfg <- gsub("\\\\", "/", DIR_RES)

for (pr in pairs) {
  t1 <- pr[1]
  t2 <- pr[2]
  cfg_file <- file.path(DIR_CFG, paste0("config_", t1, "_", t2, ".txt"))

  cfg <- c(
    paste0("reffile=", ref_mat_cfg),
    paste0("traitfolder=", mat_dir_cfg),
    paste0("traitfile1=", mat_names[t1]),
    paste0("traitname1=", t1),
    paste0("traitfiles={'", mat_names[t2], "'}"),
    paste0("traitnames={'", t2, "'}"),
    paste0("outputdir=", res_dir_cfg, "/", t1, "_", t2),
    "stattype=conjfdr",
    "fdrthresh=0.05",
    "randprune=true",
    "randprune_n=500",
    "reset_pruneidx=true",
    "randprune_repeats=default",
    "pthresh=1",
    "exclude_chr_pos=[6 26000000 34000000]",
    "exclude_from_discovery=false",
    "mafthresh=0.01",
    "exclude_ambiguous_snps=true",
    "perform_gc=true",
    "use_standard_gc=false",
    "randprune_gc=true",
    "onscreen=false",
    "dummy_zscore=false",
    "exit_matlab_upon_completion=true"
  )

  writeLines(cfg, cfg_file)
}

next_txt <- c(
  "STEP5A completed.",
  "",
  "What is ready now:",
  "1) Standardized pleioFDR text inputs: 01_inputs_prepared/",
  "2) Bidirectional screening conditional Q-Q plots: 03_conditional_QQ/",
  "3) QC/readiness tables: 00_qc/",
  "4) Six official pleioFDR config templates: 04_config_templates/",
  "",
  "Important:",
  "- MiXeR is not used in the final manuscript or this step.",
  "- The screening Q-Q plots use common SNPs from the prepared files and reproduce",
  "  the pleioFDR visual grammar (empirical -log10(q) on x, nominal -log10(p) on y).",
  "- Final manuscript-grade Q-Q / condFDR / conjFDR must use the official pleioFDR",
  "  reference + random LD pruning.",
  "- The official reference itself is EUR 1000G Phase 3, MAF >1%, 503 individuals.",
  "- Your current GWAS files do not provide a usable INFO column; therefore INFO>0.9",
  "  cannot be re-applied here. This must be stated transparently in Methods.",
  "",
  "Files to send back after this script:",
  "00_qc/STEP5A_trait_preparation_QC.csv",
  "00_qc/STEP5A_pair_overlap_QC.csv",
  "00_qc/STEP5A_condQQ_stratum_counts.csv",
  "00_qc/STEP5A_condQQ_left_shift_diagnostic.csv",
  "00_qc/STEP5A_pleioFDR_readiness.csv",
  "and the primary conditional Q-Q plots."
)

writeLines(next_txt, file.path(OUT, "README_STEP5A.txt"))

message("\n============================================================")
message("STEP5A finished.")
message("Output root: ", OUT)
message("Please first inspect: ", file.path(DIR_QC, "STEP5A_pleioFDR_readiness.csv"))
message("============================================================")
