# Reproducibility guide

## Scope

This release is intended to document and reproduce the final analytical workflow reported in the manuscript. MiXeR and superseded/debug-only analyses are intentionally excluded.

## Recommended execution

1. Clone/download the repository.
2. Download the external datasets and reference resources listed in `DATA_AVAILABILITY.md`.
3. Set `AFHFPEF_DATA_ROOT` to the local data directory.
4. Review `metadata/SOFTWARE_ENVIRONMENT.md` and install the required third-party tools.
5. Install R dependencies with `00_Setup/install_R_dependencies.R`.
6. Run scripts in `metadata/RUN_ORDER.csv` order.
7. Review each stage's QC/log outputs before continuing.
8. Reproduce the released main figure panels with `09_Figures/09_reproduce_main_figures.R`.

## Portability note

The released analysis scripts were curated from the actual Windows workflow used for the study. Many scripts retain `D:/A/data` as a historical fallback and some contain Windows-specific executable discovery logic. This is intentional provenance, not a claim of zero-configuration portability. Reviewers using another system should set the documented environment variables or edit the configuration block at the top of the relevant script.

## Validation status

Static structural checks were performed during packaging. The complete end-to-end workflow was not rerun in the packaging environment because the required external data exceed the size of this archive. See `metadata/VALIDATION_REPORT.md`.
