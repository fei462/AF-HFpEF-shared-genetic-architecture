# Code-release validation report

Release: v1.0
Date: 2026-09-24

## Curation performed

- Removed superseded V1/V2/... scripts where a validated later version existed.
- Removed current-session-only diagnostic/recovery scripts from the public workflow.
- Removed MiXeR analysis code because MiXeR was explicitly removed from the final manuscript.
- Replaced the historical pleioFDR dependency on a `STEP2_MIXER` input folder with a dedicated final full-GWAS standardization step (`03_00_prepare_full_GWAS_standardized_inputs.R`). This preserves the standardized full-GWAS inputs while eliminating MiXeR as an analytical dependency.
- Repaired the GSE255612 grouped-median type issue directly in the released full script.
- Merged the validated STEP13C scope-safe downstream implementation into the released full AF/HFpEF cross-disease integration script.
- Retained only the bounded, mismatch-aware final AF–BMI 16p11.2 SuSiE implementation.
- Added a complete LDSC h2/rg wrapper because the original working archive contained post-processing/visualization scripts but not a clean final runner.
- Added figure-reproduction code for the final data-derived panels and conceptual-schematic scripts for Figure 1a / Figure 6e.
- Added README, run order, dependency bootstrap, license, Code Availability text, checksums, and excluded-script audit.

## Automated static checks

All released `.R` files were scanned for balanced parentheses/braces/brackets and unclosed quoted strings after stripping comments. No structural syntax imbalance was detected.

## Scope of validation

This packaging environment does not contain the multi-hundred-GB external GWAS/single-cell/reference data and therefore the entire end-to-end analysis was not re-executed inside the packaging container. The analysis scripts are curated from the final working versions used during the project, and the known runtime failures encountered in the working history were removed or directly patched as described above.

The archive is intended for code peer review/reproducibility together with the manuscript Source Data and the public/external datasets cited in the README.

## Repository-readiness audit added in v1.1 (2026-09-28)

- Confirmed 65 released R scripts and a complete ordered analytical workflow for the retained analysis branches.
- Confirmed MiXeR remains excluded from the final workflow.
- Added `.gitignore` to prevent accidental upload of large raw/reference datasets, generated objects, local configuration, and credentials.
- Added data-availability, reproducibility, and GitHub/Zenodo release documentation.
- Detected historical Windows absolute-path fallbacks in the released scripts. These are documented rather than globally rewritten because indiscriminate path substitution could alter the validated analysis logic; the repository supports relocation through documented configuration/environment variables where implemented.
- The uploaded archive did not include the complete final Supplementary Figure plotting-code/source-data set. No ungrounded replacement scripts were generated.
