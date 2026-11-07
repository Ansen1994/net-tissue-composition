# Analysis code for: NET-related gene signatures in ER-positive breast cancer reflect tissue composition rather than neutrophil biology

Portable, path-independent copy of the analysis code. Every machine-specific absolute path in the original scripts has been replaced by an environment variable; the replacements are listed one by one in `PATHS.md`. The original working copies were not modified.

## Contents

| Path | What it is |
|---|---|
| `scripts/01_primary/` | Differential expression, LASSO and random-forest feature selection, Cox models, Kaplan-Meier and the original figure scripts |
| `scripts/02_drug/` | Exploratory drug-sensitivity analyses (GDSC2, ridge regression) |
| `scripts/03_work/` | Cohort assembly, module scores, single-cell attribution, table generation (`s9_supp_tables.R`, `make_table_s5_ci.R`, `make_table_s5c.R`, `make_table_s13f.R`, `make_table_s15b.R`, `j3_singlecell.R`) |
| `scripts/04_work_scripts/` | Validation-cohort and cross-platform scaling scripts |
| `scripts/05_figures_v2/` | Final figure scripts (`fig1.R`–`fig7.R`, `style.R`) |
| `scripts/06_python/` | Packaging helpers only — see "What the Python scripts do not do" below |
| `scripts/07_sc_umap/` | Standalone asset that recomputes the display-only UMAP of Supplementary Figure S3 (config, pipeline, plots, read-only loader); enters no statistical model |
| `00_config.R` | Central configuration: `WB_DATA_ROOT`, `WB_OUTPUT_DIR`, fixed seeds |
| `sessionInfo.txt` | Output of `sessionInfo()` on the machine that produced the results (R 4.4.2) |
| `PATHS.md` | Audit of every absolute path that was rewritten |
| `MANIFEST.md` | Script-by-script listing with file counts |
| `图表脚本对照表.md` | Crosswalk from each figure and table to its script |

Totals: **87** R scripts, **13** Python scripts, **22** path replacements across **15** files.

## Requirements

- R 4.4.2 with the packages recorded in `sessionInfo.txt` (DESeq2, glmnet, randomForest, survival, survminer, timeROC, clusterProfiler, ggplot2, patchwork, rms).
- Python 3.13 with openpyxl, python-docx, PyMuPDF and Pillow — needed only for packaging, not for any analysis.

## How to run

1. Download the public input data listed under "Input data" below into one directory and set `WB_DATA_ROOT` to it.
2. Set `WB_OUTPUT_DIR` to a writable directory for derived tables and figures.
3. Run the scripts in directory order (`01_primary` → `02_drug` → `03_work` → `04_work_scripts` → `05_figures_v2`), each from within its own directory after `source("../../00_config.R")`.

Both variables can be set as environment variables or edited in `00_config.R`; no script contains a machine path.

## Fixed seeds

- `2026` — LASSO (`cv.glmnet`) and random-forest feature selection.
- `42` — five-fold cross-validation, module scores and drug-sensitivity analyses.
- `20261106` — display-only UMAP embedding of Supplementary Figure S3.

## Input data

All inputs are public and were obtained at no cost:

- TCGA-BRCA — Genomic Data Commons (https://portal.gdc.cancer.gov/)
- GSE20685, GSE96058, GSE176078 — NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/)
- METABRIC — cBioPortal (https://www.cbioportal.org/)
- GDSC2 — Genomics of Drug Sensitivity in Cancer (https://www.cancerrxgene.org/)
- NET-related gene annotations — GeneCards (https://www.genecards.org/, relevance score > 10)

## What the Python scripts do not do

The Python scripts perform no statistical computation and produce no number that appears in the manuscript. They consolidate the supplementary workbook, convert figure PDFs to 600 dpi submission TIFF, audit reference records and assemble the manuscript package. Every numerical result was produced in R.

## Licence

Released under the MIT License; the full text is in `LICENSE`. A permissive licence was chosen so that the code can be reused and adapted; the manuscript's text and figures are not covered by it.
