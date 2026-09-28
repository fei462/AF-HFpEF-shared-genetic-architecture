# Software environment and third-party dependencies

The custom code is primarily R. External programs are called where indicated by the individual scripts.

## Core environment

- R: analyses were developed and run under R 4.x; scripts write `sessionInfo()` or package-version tables at critical stages where possible.
- Python 3: used for LDSC. The project used a Python-3-compatible LDSC environment.
- MATLAB: used only for the official pleioFDR/conjFDR implementation; Statistics and Machine Learning Toolbox is required by the upstream package.

## External command-line/software dependencies

| Tool | Role | Upstream source |
|---|---|---|
| LDSC | SNP-h2, genetic correlation, cell-type S-LDSC | https://github.com/bulik/ldsc |
| pleioFDR | conditional FDR / conjFDR | https://github.com/precimed/pleiofdr |
| PLINK 1.9 | LD reference manipulation and LD matrices | https://www.cog-genomics.org/plink/1.9/ |
| FUMA | functional mapping and locus annotation | https://fuma.ctglab.nl/ |
| LAVA | local genetic correlation | https://github.com/josefin-werme/LAVA |
| SMR 1.3.x | SMR/HEIDI analyses | https://yanglab.westlake.edu.cn/software/smr/ |
| SCAVENGE | single-cell trait relevance | https://github.com/sankaranlab/SCAVENGE |
| gchromVAR | single-cell regulatory scoring | https://github.com/caleblareau/gchromVAR |

The code archive does not redistribute third-party software. Reviewers should install the official tools from the upstream sources.

## Reference resources

- 1000 Genomes Phase 3 European LD references
- LDSC baselineLD / HapMap3 weights and regression-score resources
- LAVA UK Biobank reference resources
- GTEx v8 eQTL/sQTL summary resources used by the SMR/HEIDI layer

Large reference files are intentionally excluded from this code archive.
