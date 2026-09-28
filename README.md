# Code release: AF–HFpEF–BMI–OSA shared genetic architecture study

Version: **1.1 (repository-ready curation, 2026-09-28)**

This repository contains the curated custom code used for the final manuscript analyses and figure generation. It was assembled from the complete working-code archive, with obsolete versions, failed debugging attempts, current-session hotfixes, and analyses not retained in the final manuscript removed.

## Important final-analysis decisions

- **MiXeR is not part of the final manuscript and is intentionally excluded from this release.**
- The primary genetic analyses retained are LDSC, conditional Q-Q/conjFDR, FUMA locus annotation, MR/MVMR, LAVA, coloc/SuSiE, SMR/HEIDI, S-LDSC, gchromVAR/SCAVENGE, and patient-level pseudobulk transcriptomic validation.
- Figure 4b from an earlier working draft was judged redundant with Figure 4a and is not reproduced as a final main-text panel.
- The AF–BMI 16p11.2 analysis is explicitly a **secondary/post-hoc** follow-up and is kept separate from the pre-specified coloc/SuSiE tests.

## Repository layout

- `00_Setup/` — dependency bootstrap and configuration example
- `01_GWAS_QC_LDSC/` — GWAS audit, HFpEF extraction, harmonization, LDSC h2/rg, LDSC post-processing
- `02_MR_MVMR/` — bidirectional UVMR, MVMR, weak-instrument sensitivity, final MR tables/QC
- `03_Pleiotropy_conjFDR_FUMA/` — conditional Q-Q, six pairwise conjFDR analyses, FUMA input preparation
- `04_LAVA/` — local h2, bivariate local rg, multivariate/conditional LAVA
- `05_Coloc_SuSiE_SMR_HEIDI/` — coloc ABF, prior sensitivity, mismatch-aware SuSiE, eQTL/sQTL SMR/HEIDI
- `06_Multilayer_Locus_Integration/` — multilayer locus/gene integration and secondary 16p11.2 follow-up
- `07_GSE238242_SLDSC_SCAVENGE/` — left-atrial snATAC workflow, S-LDSC, fine-mapping, gchromVAR, SCAVENGE, sensitivity analyses
- `08_Disease_Tissue_Validation/` — AF and HFpEF patient-level cardiomyocyte pseudobulk validation and cross-disease integration
- `09_Figures/` — panel-level Source Data and scripts to regenerate data-derived main figures plus conceptual schematics
- `metadata/` — code manifest, excluded-script audit, run order, validation report

## Quick start

1. Use R from the repository root.
2. Optionally set the project data root:

```r
Sys.setenv(AFHFPEF_DATA_ROOT = "D:/A/data")
```

3. Install dependencies:

```r
source("00_Setup/install_R_dependencies.R")
```

4. Execute scripts in the order given in `metadata/RUN_ORDER.csv`.

The released scripts preserve the original project fallback path `D:/A/data` because that is the path used for the reported analyses. Several scripts also support environment-variable overrides. Large public reference resources and raw GWAS/single-cell datasets are **not duplicated in this code archive**.

## External software and resources

The workflow uses third-party tools in addition to R:

- LDSC (Python 3 compatible implementation): https://github.com/bulik/ldsc
- pleioFDR/conjFDR: https://github.com/precimed/pleiofdr
- MATLAB with Statistics and Machine Learning Toolbox for pleioFDR
- PLINK 1.9: https://www.cog-genomics.org/plink/1.9/
- FUMA: https://fuma.ctglab.nl/
- LAVA: https://github.com/josefin-werme/LAVA
- SMR/HEIDI: https://yanglab.westlake.edu.cn/software/smr/
- SCAVENGE: https://github.com/sankaranlab/SCAVENGE
- gchromVAR: https://github.com/caleblareau/gchromVAR

Reference data include 1000 Genomes Phase 3 EUR LD resources, LDSC baseline/weights, and UK Biobank-based LAVA reference files as described in the corresponding scripts and manuscript Methods.

## Public datasets referenced by the code

- AF GWAS: GCST90624412 (European ancestry dataset used in this project)
- HFpEF GWAS: HERMES HFpEF resource used by the project
- BMI GWAS: GIANT/Locke 2015 summary statistics used by the project
- OSA GWAS: FinnGen R9 G6_SLEEPAPNO
- Human left-atrial multiome: GEO **GSE238242**
- Human AF left-atrial snRNA-seq: GEO **GSE255612**
- Human HFpEF snRNA-seq: Single Cell Portal **SCP3342**

Data access is governed by the original repositories. The code does not redistribute controlled or third-party data.

## Reproducibility notes

- Most analysis stages write their own QC tables, manifests, logs, and/or `sessionInfo()` to the project output folder.
- The final code release removes known failed versions and retains the highest validated working version for each analysis branch.
- Two historical runtime problems are repaired directly in the release rather than requiring reviewer-side hotfixes:
  - GSE255612 donor×cell-type summary now coerces grouped median cell counts to numeric to avoid `data.table` type instability.
  - AF/HFpEF cross-disease integration uses the validated scope-safe downstream implementation after permutation testing.
- The secondary AF–BMI 16p11.2 SuSiE workflow uses bounded mismatch-aware fitting and records non-convergence as unresolved rather than forcing convergence.

## Figure reproduction

The directory `09_Figures/source_data/` contains the panel-level CSVs used for the final figure scripts. Run:

```r
source("09_Figures/09_reproduce_main_figures.R")
```

Outputs are written to `09_Figures/reproduced_panels/` at 600 dpi. Figure 1a and Figure 6e are conceptual schematics and have separate scripts.

## Code availability statement for submission

See `CODE_AVAILABILITY_STATEMENT.txt`.

## License

The custom code in this release is provided under the MIT License (`LICENSE`). Third-party packages and software retain their own licenses.


## Repository-release audit (v1.1)

The uploaded code archive was reviewed for repository readiness. The analytical workflow is internally organized and contains 65 R scripts covering the retained final analyses. Additional repository files were added in v1.1: `.gitignore`, `DATA_AVAILABILITY.md`, `REPRODUCIBILITY.md`, and `GITHUB_ZENODO_RELEASE_CHECKLIST.md`.

A complete final Supplementary Figure plotting-code set was **not present in the uploaded archive**. No synthetic replacement scripts were created. This limitation is documented so that the public repository does not imply reproducibility beyond the supplied materials.
