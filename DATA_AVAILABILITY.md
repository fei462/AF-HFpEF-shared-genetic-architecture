# Data availability and expected inputs

This repository contains code and small figure Source Data files only. It does **not** redistribute the large third-party GWAS, LD-reference, QTL, or single-cell datasets.

## Public/external datasets used

- AF GWAS: GCST90624412 (European-ancestry dataset used in this project)
- HFpEF GWAS: HERMES HFpEF resource used in the project
- BMI GWAS: GIANT/Locke 2015 summary statistics
- OSA GWAS: FinnGen R9 G6_SLEEPAPNO
- Left-atrial multiome: GEO GSE238242
- AF left-atrial snRNA-seq: GEO GSE255612
- HFpEF snRNA-seq: Single Cell Portal SCP3342
- LD/reference resources: 1000 Genomes Phase 3 EUR; LDSC reference/baseline resources; UK Biobank-based LAVA reference
- Regulatory QTL resources: GTEx v8 cis-eQTL and cis-sQTL resources used by the SMR/HEIDI workflow

Access conditions remain those of the original repositories. Users should download these resources from their original providers and point the scripts to a local data root.

## Local configuration

The historical analysis root was `D:/A/data`. For a relocated installation, define the environment variable before running the workflow:

```r
Sys.setenv(AFHFPEF_DATA_ROOT = "/path/to/data")
```

See `00_Setup/config_example.R`, individual script headers, and `metadata/RUN_ORDER.csv`.

## Figure Source Data

Small CSV files required by the released main-figure reproduction script are included under `09_Figures/source_data/`.

**Important:** the submitted archive did not contain a complete set of final Supplementary Figure plotting scripts/source-data tables. They are therefore not fabricated in this release. If those final scripts are to be deposited, add the exact versions used for the submitted Supplementary Figures and update the manifest/checksums before creating the archival release.
