# Install R dependencies used by the released scripts.
# Run once in a clean R >= 4.2 environment.
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c(
  "data.table", "ggplot2", "ggrepel", "dplyr", "tidyr", "stringr",
  "Matrix", "DBI", "duckdb", "hdf5r", "zip", "remotes", "BiocManager", "patchwork"
)
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

# MR ecosystem
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
}
if (!requireNamespace("MVMR", quietly = TRUE)) remotes::install_github("WSpiller/MVMR", upgrade = "never")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) remotes::install_github("rondolab/MR-PRESSO", upgrade = "never")
if (!requireNamespace("LAVA", quietly = TRUE)) remotes::install_github("josefin-werme/LAVA", upgrade = "never")

# Colocalization / fine-mapping
for (p in c("coloc", "susieR")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# Bioconductor / single-cell regulatory analysis
bioc <- c("GenomicRanges", "IRanges", "S4Vectors", "SummarizedExperiment", "rtracklayer", "chromVAR")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("gchromVAR", quietly = TRUE)) remotes::install_github("caleblareau/gchromVAR", upgrade = "never")
if (!requireNamespace("SCAVENGE", quietly = TRUE)) {
  # SCAVENGE is distributed through its project repository; install from the source
  # specified in the README if not already available in your environment.
  message("SCAVENGE is not installed. See README.md > External software.")
}

cat("Dependency bootstrap complete.\n")
